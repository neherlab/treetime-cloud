from __future__ import print_function, division
import treetime
import pandas
import numpy as np
from Bio import Phylo, AlignIO, Align, Seq, SeqRecord
import zipfile
import time
import os, sys, json
from treetime import utils
from treetime.treetime import TreeTime
import traceback

dirname = (os.path.dirname(__file__))

treename = "in_tree.nwk"
alnname = "in_aln.fasta"
metaname = "in_meta.csv"
session_state_name = "session_state.txt"


# file names, uniform across the sessions
session_state_file = 'session_state.json'
in_tree = "in_tree.nwk"
in_aln = "in_aln.fasta"
in_meta = "in_meta.csv"
in_cfg = "config.json"

log_filename = "log.txt"
out_tree_json = "out_tree.json"
out_likelihoods_json = "out_likelihoods.json"
out_tree_nwk = "out_tree.nwk"
out_aln_fasta = "out_aln.fasta"
out_metadata_csv = "out_metadata.csv"
out_mol_clock_csv = 'molecular_clock.csv'
out_gtr = "out_GTR.txt"
zipname = 'treetime_results.zip'


def get_filepaths(root):
    """
    Get file locations to save the user-uploaded files and to run treetime
    """
    return {
        'tree' : os.path.join(root, treename),
        'aln' : os.path.join(root, alnname),
        'meta' : os.path.join(root, metaname),
        "session_state":os.path.join(root, session_state_name)
    }


def read_metadata_from_file(infile):
    """
    @brief      Reads a metadata from file or handle.
    @param      self    The object
    @param      infile  The input file or file handle
    """
    try:
        # read the metadata file into pandas dataframe.
        df = pandas.read_csv(infile, index_col=0, sep=r'\s*,\s*')
        # check the metadata has strain names in the first column
        if 'name' not in df.index.name.lower():
            print ("Cannot read metadata: first column should contain the names of the strains")
            return
        # look for the column containing sampling dates
        # We assume that the dates might be given in eihter human-readable format
        # (e.g. ISO dates), or be already converted to the numeric format.
        potential_date_columns = []
        potential_numdate_columns = []
        # Scan the dataframe columns and find ones which likely to store the
        # dates
        for ci,col in enumerate(df.columns):
            if 'date' in col.lower():
                try: #  avoid date parsing when can be parsed as float
                    tmp = float(df.iloc[0,ci])
                    potential_numdate_columns.append((ci, col))
                except: #  otherwise add as potential date column
                    potential_date_columns.append((ci, col))
        # if a potential numeric date column was found, use it
        # (use the first, if there are more than one)
        if len(potential_numdate_columns)>=1:
            name = potential_numdate_columns[0][1]
            # Use this column as numdate_given
            dates = df[name].to_dict()
        elif len(potential_date_columns)>=1:
            #try to parse the csv file with dates in the idx column:
            idx = potential_date_columns[0][0]
            name = potential_date_columns[0][1]
            # NOTE as the 0th column is the index, we should parse the dates
            # for the column idx + 1
            df = pandas.read_csv(infile, index_col=0, sep=r'\s*,\s*', parse_dates=[1+idx])
            dates = {k: utils.numeric_date(df.loc[k, name]) for k in df.index}
            df.loc[:, name] = map(lambda x: str(x.date()), df.loc[:, name])
        else:
            print ("Metadata file has no column which looks like a sampling date!")
        metadata = df.to_dict(orient='index')
        return dates, metadata
    except:
        print ("Cannot read the metadata file. Exception caught")
        raise
        return {}, {}

class SessionState(object):
    saving = "saving results"
    error = "error"
    running = "running"
    reading = "reading config"
    done = "done"

def _write_session_state(root, state, desc=""):

        dic = {"state":state,
                "desc": desc}
        with open (os.path.join(root, "session_state.txt"), 'w') as of:
            json.dump(dic, of, indent=True)

class TreeTimeWeb(treetime.TreeTime):

    def __init__(self, root, webconfig, metadata=True, *args, **kwargs):

        self._webconfig = webconfig

        self._root_dir = root
        self._log_file = os.path.join(self._root_dir, log_filename)

        if webconfig['build_tree'] is True or webconfig['build_tree'] == 'True':
            self.build_tree(root)

        # run treetime with the specified parameters:
        tree = Phylo.read(get_filepaths(root)['tree'], 'newick')
        aln = AlignIO.read(get_filepaths(root)['aln'], 'fasta')

        if metadata:
            dates, metadata = read_metadata_from_file(get_filepaths(root)['meta'])
            self._metadata = metadata
        else:
            dates = {}

        gtr = 'jc' if webconfig['gtr'] == 'infer' else webconfig['gtr']
        super(TreeTimeWeb, self).__init__(dates=dates, tree=tree, aln=aln,
                gtr=str(gtr), *args, **kwargs)

    def run(self, **kwargs):
        _write_session_state(self._root_dir, SessionState.reading)
        # get the run parameters
        infer_gtr  = self._webconfig['gtr'] == 'infer'
        root = self._webconfig['root']
        do_marginal = False if self._webconfig['do_marginal'] == 'False' or not self._webconfig['do_marginal'] else True
        resolve_polytomies = False if self._webconfig['polytomies'] == 'False' or not self._webconfig['polytomies'] else True
        slope = None if self._webconfig['slope'] == 'False' or not self._webconfig['slope'] else float(self._webconfig['slope_value'])
        Tc = None if self._webconfig['use_coalescent_prior'] == 'False' or not self._webconfig['use_coalescent_prior'] \
            else float(self._webconfig['coalescent_prior_value'])

        relax_clock = False if self._webconfig['use_relaxed_clock'] == 'False' or not self._webconfig['use_relaxed_clock'] else True
        # add kwargs only if the relax clock is set to true
        if relax_clock:
            kwargs['slack'] = float(self._webconfig['relaxed_clock']['slack'])
            kwargs['coupling'] = float(self._webconfig['relaxed_clock']['coupling'])

        # run treetime
        try:
            _write_session_state(self._root_dir, SessionState.running)
            super(TreeTimeWeb, self).run(root=root, infer_gtr=infer_gtr, relaxed_clock=relax_clock,
                resolve_polytomies=resolve_polytomies, max_iter=5, Tc=Tc, fixed_slope=slope,
                do_marginal=do_marginal, **kwargs)
        except:
            tb = traceback.format_exc()
            _write_session_state(self._root_dir, SessionState.error, desc="TreeTime crashed. {}".format(tb))
            return

        # save results
        try:
            self.logger("###TreeTimeWeb.run: Done treetime computations, saving the results",0)
            _write_session_state(self._root_dir, SessionState.saving)
            self.save_treetime_results()
            _write_session_state(self._root_dir, SessionState.done)
            self.logger("###TreeTimeWeb.run: All tasks completed successfully, exiting...",0)
        except:
            tb = traceback.format_exc()
            _write_session_state(self._root_dir, SessionState.error, desc="TreeTime crashed while saving results. {}".format(tb))
            return


    def save_treetime_results(self):

        from Bio import Align
        #  files to be displayed in the web interface
        self._tree_to_json()
        self._likelihoods_to_json()

        # files to be downloaded as .zip archive
        Phylo.write(self.tree, os.path.join(self._root_dir, out_tree_nwk), 'newick')
        self._save_alignment()
        self._save_metadata_to_csv()
        self._save_molecular_clock_to_csv()
        self._save_gtr()
        # zip all results to one file
        with zipfile.ZipFile(os.path.join(self._root_dir, zipname), 'w') as out_zip:
            out_zip.write(os.path.join(self._root_dir, out_tree_nwk), arcname=out_tree_nwk)
            out_zip.write(os.path.join(self._root_dir, out_aln_fasta), arcname=out_aln_fasta)
            out_zip.write(os.path.join(self._root_dir, out_metadata_csv), arcname=out_metadata_csv)
            out_zip.write(os.path.join(self._root_dir, out_tree_json), arcname=out_tree_json)
            #out_zip.write(os.path.join(self._root_dir, in_cfg), arcname=in_cfg)
            out_zip.write(os.path.join(self._root_dir, out_mol_clock_csv), arcname=out_mol_clock_csv)
            out_zip.write(os.path.join(self._root_dir, out_likelihoods_json), arcname=out_likelihoods_json)
            out_zip.write(os.path.join(self._root_dir, out_gtr), arcname=out_gtr)

    def _save_alignment(self):
        aln = Align.MultipleSeqAlignment([SeqRecord.SeqRecord (Seq.Seq(''.join(n.sequence)), id=n.name, name=n.name, description="")
            for n in self.tree.find_clades ()])
        AlignIO.write(aln, os.path.join(self._root_dir, out_aln_fasta), "fasta")

    def _save_metadata_to_csv(self):
        meta = {node: self._node_metadata(node) for node in self.tree.find_clades()}
        rows_dic = {key.name: {e['name'] : e['value']for e in meta[key]} for key in meta.keys()}
        df = pandas.DataFrame(rows_dic).T
        outf = os.path.join(self._root_dir, out_metadata_csv)
        df.to_csv(outf)

    def _save_molecular_clock_to_csv(self):
        #save molecular clock in normal format
        mclock = np.array([(tip.dist2root, tip.numdate_given)
            for tip in self.tree.get_terminals()
            if hasattr(tip, 'dist2root') and hasattr(tip, 'numdate_given')])
        np.savetxt(os.path.join(self._root_dir, out_mol_clock_csv), mclock,
            delimiter=',',
            header='Distance_to_root,Sampling_date')

    def _save_gtr(self):
        with open (os.path.join(self._root_dir, out_gtr), 'w') as outf:
            outf.write(str(self.gtr))

    def _tree_to_json(self):
        """
        Add extra attributes to the tree, and save the tree in JSON format
        """

        def _node_to_json(node):
            """
            convert node data to dictionary
            """

            tree_json = {}
            str_attr = ['clade','strain', 'date', 'muts', 'strseq']
            num_attr = ['xvalue', 'yvalue', 'tvalue', 'numdate']

            if hasattr(node, 'name'):
                tree_json['strain'] = node.name
                tree_json['name'] = node.name

            for prop in str_attr:
                if hasattr(node, prop):
                    tree_json[prop] = node.__getattribute__(prop)
            for prop in num_attr:
                if hasattr(node, prop):
                    try:
                        tree_json[prop] = round(node.__getattribute__(prop),5)
                    except:
                        self.logger("cannot round:", node.__getattribute__(prop), "assigned as is", 1, warn=True)
                        tree_json[prop] = node.__getattribute__(prop)

            if node.clades: # node is internal
                tree_json["children"] = []
                for ch in node.clades:
                    tree_json["children"].append(_node_to_json(ch))

            tree_json["metadata"] = self._node_metadata(node)

            return tree_json

        self._layout()
        tree_json = _node_to_json(self.tree.root)

        # save the result in the json file:
        outf = os.path.join(self._root_dir, out_tree_json)
        with open (outf,'w') as of:
            json.dump(tree_json, of, indent=False)

    def _layout(self):
        """Add clade, xvalue, yvalue, mutation and trunk attributes to all nodes in tree"""
        self.tree.root.branch_length = 0.001
        clade = 0
        yvalue = 0
        for node in self.tree.find_clades(order="preorder"):
            # set mutations
            if node.up is not None:
                node.muts = ', '.join([node.up.sequence[p] + str(p) + node.sequence[p]
                    for p in np.where(node.up.sequence != node.sequence)[0]])

            # set sequences
            node.strseq = "".join(node.sequence)

            # set clade No
            node.clade = clade
            clade += 1
            if node.up is not None: #try:
                # Set xValue, tValue, yValue
                node.xvalue = node.up.xvalue + node.mutation_length
                node.tvalue = node.numdate - self.tree.root.numdate
            else:
                node.xvalue = 0.0
                node.tvalue = 0.0
            if node.is_terminal():
                node.yvalue = yvalue
                yvalue += 1
            # check numdate
            if not hasattr(node, 'numdate'):
                node.numdate = 0.0
        for node in self.tree.get_nonterminals(order="postorder"):
            node.yvalue = np.mean([x.yvalue for x in node.clades])

    def _node_metadata(self, node):


        if hasattr(node, 'branch_length_interpolator') and node.branch_length_interpolator is not None:
            gamma = node.branch_length_interpolator.gamma
        else:
            gamma = 1.0

        meta = []
        # if node has user-provided metadata, append it
        if node.name in self._metadata:
            node_meta = self._metadata[node.name]
            meta += [{'name':k, 'value': node_meta[k]} for k in node_meta]

        # append numdates to the metadata
        if hasattr(node, 'numdate'):
            meta.append({"name":"numdate", "value":node.numdate})

        # append deviation of the branch length from average
        stretch = np.min([node.branch_length / (node.mutation_length / gamma + 0.000001), 2.])
        meta.append({
                "name": "Branch length stretch",
                "value": stretch
            })

        relax_clock = False if self._webconfig['use_relaxed_clock'] == 'False' or not self._webconfig['use_relaxed_clock'] else True
        if relax_clock:
            # append mutation rate deviation from average
            meta.append({
                "name": "Local mutation rate",
                "value": gamma
                })
            # else:
            #     meta.append({"name": "Relaxed mutation rate", "value":1.0})

        return meta

    def _distribution_to_human_readable(self, dist):

        if dist.is_delta:
            date = self.date2dist.to_numdate(dist.peak_pos)
            return [date-0.5, date-1e-10, date, date + 1e-10, date + 0.5],[0,0,1.0,0,0]

        peak_pos = dist.peak_pos
        fwhm = dist.fwhm
        raw_x = dist.x [(dist.x > peak_pos - 3*fwhm) & (dist.x < peak_pos + 3*fwhm)]
        dates_x = np.array(map(self.date2dist.to_numdate, raw_x))
        y = dist.prob_relative(raw_x)
        return dates_x, y

    def _likelihoods_to_json(self):
        """
        Save the likelihoods fro the node positions in JSON format
        """
        out_dic = {}
        for node in self.tree.find_clades():
            x,y = self._distribution_to_human_readable(node.marginal_pos_LH)
            arr = [{"x":f, "y":b} for f, b in zip(x, y)]
            out_dic[node.name] = arr

        outfile = os.path.join(self._root_dir, out_likelihoods_json)
        with open(outfile, 'w') as outf:
            json.dump(out_dic, outf, indent=False)

    def _tips_data_to_json(self):

        if not hasattr(self.tree.get_terminals()[0], 'xvalue'):
            self._layout()

        arr = [
        {
            'name': k.name,
            'strain':k.name,
            'numdate_given': k.numdate_given if hasattr(k, 'numdate_given') else 0.0,
            'numdate': k.numdate if hasattr(k, 'numdate') else 0.0,
            'xValue': k.xvalue if hasattr(k, 'xvalue') else 0.0,

        } for k in tt.tree.get_terminals()]

        with open (outf,'w') as of:
            json.dump(arr, of, indent=True)

    def logger(self, msg, level, warn=False):
        """
        @brief      Overriding the basic logger functionality to enable possibility to
                    redirect the output to the file

        @param      self   The object
        @param      msg    The message
        @param      level  The verbosity level (1-highest, 10-lowest)
        @param      warn   Warning flag. If set to True, the message will be printed
                            regardless the verbosity with Warn mark

        """
        if level<self.verbose or warn:
            dt = time.time() - self.t_start
            outstr = '\n' if level<2 else ''
            outstr+= format(dt, '4.2f')+'\t'
            outstr+= level*'-'
            outstr+=msg
            with open(self._log_file, 'a') as logf:
                print(outstr,file=logf)
                print (outstr)

    def build_tree(self, root):

        aln_filename = get_filepaths(root)['aln']
        tree_filename = get_filepaths(root)['tree']

        fast_tree_bin = os.path.join(dirname, "fasttree")
        if not os.path.exists(fast_tree_bin):
            raise (RuntimeError("Cannot find FastTree binary."))

        call = [fast_tree_bin, '-nt','-quiet', aln_filename, ' > ', tree_filename]

        res = os.system(' '.join(call))
        if res != 0:
            raise RuntimeError("FastTree: Exception caught while building the phylogenetic tree")

        return

    def run_treeanc(self, **kwargs):
        _write_session_state(self._root_dir, SessionState.reading)
        # get the run parameters
        infer_gtr  = self._webconfig['gtr'] == 'infer'
        do_marginal = False if self._webconfig['do_marginal'] == 'False' or not self._webconfig['do_marginal'] else True

        # run treeanc
        try:
            _write_session_state(self._root_dir, SessionState.running)
            super(TreeTimeWeb, self).optimize_seq_and_branch_len(reuse_branch_len=True, prune_short=True, max_iter=5, infer_gtr=infer_gtr, marginal=do_marginal, **kwargs)
        except:
            tb = traceback.format_exc()
            _write_session_state(self._root_dir, SessionState.error, desc="TreeTime crashed. {}".format(tb))
            return

        # save results
        try:
            self.logger("###TreeTimeWeb.run: Done treetime computations, saving the results",0)
            _write_session_state(self._root_dir, SessionState.saving)
            self.save_treeanc_results()
            _write_session_state(self._root_dir, SessionState.done)
            self.logger("###TreeTimeWeb.run: All tasks completed successfully, exiting...",0)
        except:
            tb = traceback.format_exc()
            _write_session_state(self._root_dir, SessionState.error, desc="TreeTime crashed while saving results. {}".format(tb))
            return

    def save_treeanc_results(self):
        from Bio import Align
        #  files to be displayed in the web interface
        Phylo.write(self.tree, os.path.join(self._root_dir, out_tree_nwk), 'newick')
        self._save_alignment()
        self._save_gtr()
        with zipfile.ZipFile(os.path.join(self._root_dir, zipname), 'w') as out_zip:
            out_zip.write(os.path.join(self._root_dir, out_tree_nwk), arcname=out_tree_nwk)
            out_zip.write(os.path.join(self._root_dir, out_aln_fasta), arcname=out_aln_fasta)
            out_zip.write(os.path.join(self._root_dir, out_gtr), arcname=out_gtr)

def run_treetime(root, webconfig):
    try:
        ttw = TreeTimeWeb(root, webconfig)
    except:
        tb = traceback.format_exc()
        _write_session_state(root,SessionState.error, desc="TreeTime crashed at the initialization phase. {}".format(tb))
        return
    ttw.run()

def run_treeanc(root, webconfig):
    try:
        ttw = TreeTimeWeb(root, webconfig, metadata=False)
    except:
        tb = traceback.format_exc()
        _write_session_state(root, SessionState.error, desc="TreeTime crashed at the initialization phase. {}".format(tb))
        return

    ttw.run_treeanc()

if __name__=="__main__":

    root = '../../sessions/AGJPHWRULXGO/'
    import  tree_time_config
    cfg = tree_time_config.treetime_webconfig

    ttw = TreeTimeWeb(root, cfg)
    ttw.run()



