from treetime.treetime import TreeTime
import os, json
from Bio import Phylo, AlignIO, Align, Seq, SeqRecord
import numpy as np
from treetime import utils
import pandas
import zipfile

myDir = os.path.dirname(os.path.abspath(__file__))

available_gtrs = ['Jukes-Cantor']

# file names, uniform across the sessions
session_state_file = 'session_state.json'
in_tree = "in_tree.nwk"
in_aln = "in_aln.fasta"
in_meta = "in_meta.csv"
in_cfg = "config.json"

out_tree_json = "out_tree.json"
out_likelihoods_json = "out_likelihoods.json"
out_tree_nwk = "out_tree.nwk"
out_aln_fasta = "out_aln.fasta"
out_metadata_csv = "out_metadata.csv"
out_mol_clock_csv = 'molecular_clock.csv'
out_gtr = "out_GTR.txt"
zipname = 'treetime_results.zip'


def build_tree(root):


    aln_filename = os.path.join(root, "in_aln.fasta")
    tree_filename = os.path.join(root, "in_tree.nwk")
    call = [os.path.join(myDir, 'fasttree'),
            '-nt','-quiet', aln_filename, ' > ', tree_filename]
    os.system(' '.join(call))

class TreeTimeWeb(TreeTime):
    """
    TreeTime class extension to be used with web interface.
    This class defines functionality to communicate with the server, provide an
    arbitrary data to display

    """

    def __init__(self, root_dir, *args,**kwargs):

        # set the default configurations
        self._build_tree = False
        self._infer_gtr = True
        self._root = None
        self._mutation_rate = None
        self._relaxed_clock = None
        self._Tc = None

        # set the directory, containing files for the session
        self._root_dir = root_dir
        self._nwk =  os.path.join(self._root_dir, in_tree)
        self._aln =  os.path.join(self._root_dir, in_aln)
        self._meta = os.path.join(self._root_dir, in_meta)
        self._cfg  = os.path.join(self._root_dir, in_cfg)

        # read the JSON configuration file
        with open(self._cfg) as ff:
            config_dic = json.load(ff)

        # Compose the list of steps to perform,
        # # save this as dictionary to session_state.json file (used to update
        # web browser state)
        self._init_session_state(config_dic)

        #  build tree if necessary
        if self._build_tree:
            print ("Building phylogenetic tree...")
            build_tree(self._root_dir)
            self._advance_session_progress()

        tree = Phylo.read(self._nwk, 'newick')
        aln = AlignIO.read(self._aln, 'fasta')
        #  read the metadata
        dates, metadata = self._read_metadata_from_file(self._meta)
        super(TreeTime, self).__init__(dates=dates, tree=tree, aln=aln,
            gtr=self._gtr, *args, **kwargs)

        self._metadata = metadata

    def _init_session_state(self, config_dic):

        session_state = {
            "state" : "running",
            "error" : "", #  no errors so far
            "todo" : [],
            "done" : [],
            "progress":""
            }

        # Do we need to build the tree?
        if 'build_tree' in config_dic and config_dic['build_tree']:
            session_state['todo'].append('Build phylogenetic tree')
            self._build_tree = True


        # set the GTR model
        if 'gtr' in config_dic:
            if config_dic['gtr'] == 'infer':
                self._infer_gtr = True
                self._gtr = "Jukes-Cantor"
                session_state['todo'].append('Reconstruct ancestral states, infer evolutionary model')

            elif config_dic['gtr'] in available_gtrs:
                self._infer_gtr = False
                self._gtr = config_dic['gtr']
                session_state['todo'].append('Reconstruct ancestral states')

            else:
                # NOTE cannot use the logger class before the super.__init__ is called
                # therefore, just use plain print
                print ("Configuration contains unknown GTR type. Ignoring, will use Jukes-Cantor.")
                self._infer_gtr = False
                self._gtr = "Jukes-Cantor"
                session_state['todo'].append('Reconstruct ancestral states')

        else:
            print ("Configuration has no gtr model set. Will use default (Jukes-Cantor)")
            self._infer_gtr = False
            self._gtr = "Jukes-Cantor"
            session_state['todo'].append('Reconstruct ancestral states')

        # Find best root
        if 'reroot' in config_dic and 'reroot' is not None:
            session_state['todo'].append('Find better root for the input tree')
            self._root = config_dic['reroot']

        if 'use_mu' in config_dic and config_dic['use_mu'] == True:
            try:
                self._mutation_rate=float(config_dic['mu'])
            except:
                self.logger("Cannot parse mutation rate from the configuration dictionary!", 1, warn=True)
                self._mutation_rate = None

            if self._mutation_rate == 0:
                self.logger("TreeTimeWeb got configuration 'use_mu'= True , nut the mutation rate is set to ZERO", 1, warn=True)
                self._mutation_rate = None

            if self._mutation_rate is not None and self._root == 'best' :
                self.logger("Warning! the mutation rate cannot be used together "
                    "with the tree rooting to optimize molecular clock! Setting mutation rate to None", 1, warn=True)
                self._mutation_rate = None


        # in any case, we set the make time tree as todo
        session_state['todo'].append('Make time tree')

        # Resolve polytomies
        if 'resolve_poly' in config_dic and config_dic['resolve_poly']:
            self._resolve_polytomies = True
            session_state['todo'].append("Resolve multiple mergers")
        else:
            self._resolve_polytomies = False

        # Model coalescent process
        if 'model_coalescent' in config_dic and config_dic['model_coalescent']==True:
            try:
                self._Tc = float(config_dic['coalescent'])
                if self._Tc == 0:
                    self.logger("TreeTimeWeb got configuration 'model_coalescence'= True , "
                        "but the coalescence timescale is set to ZERO", 1, warn=True)
                    self._Tc = None
                else:
                    session_state['todo'].append("Model coalescent process")
            except:
                self.logger("Cannot parse coalescent timescale from the config dictionary!", 1, warn=True)
                self._Tc = None

        # Relax molecular clock
        if 'do_relax_clock' in config_dic and config_dic['do_relax_clock']==True:
            try:
                rc = config_dic['relax_clock']
                slack, coupling = float(rc['slack']), float(rc['coupling'])
                if slack == 0 and coupling == 0:
                    self._relaxed_clock = None
                else:
                    self._relaxed_clock = (slack, coupling)
                    session_state['todo'].append('Relax molecular clock')
            except:
                self.logger("Cannot parse relaxed clock parameters from the config dict!", 1, warn=True)
                self._relaxed_clock = None


        session_state['todo'].append('Save results')

        self._session_state = session_state
        self._save_session_state()

    def _advance_session_progress(self):

        # current step is done:
        if self._session_state['progress']:
            self._session_state['done'].append(self._session_state['progress'])

        # if there are steps to do - pop first of them and show as running
        if len(self._session_state['todo']) > 0:
            self._session_state['progress'] = self._session_state['todo'].pop(0)

        else: #  otherwise, the session state is complete:
            self._session_state['state'] = 'complete'

        # save the state to the json file:
        self._save_session_state()

    def _save_session_state(self):
        """
        Save session state dictionary to the file in the session folder
        """
        with open(os.path.join(self._root_dir, session_state_file), 'w') as outf:
            json.dump(self._session_state, outf)

    def run(self, max_iter=1):
        self._advance_session_progress()

        self.optimize_seq_and_branch_len(infer_gtr=self._infer_gtr, prune_short=True)
        self._advance_session_progress()

        if self._root is not None:
            self.reroot(root=self._root)
            self._advance_session_progress()

        self.make_time_tree(slope=self._mutation_rate)
        self._advance_session_progress()

        if self._resolve_polytomies:
            # if polytomies are found, rerun the entire procedure
            if self.resolve_polytomies():
                self.prepare_tree()
                self.optimize_seq_and_branch_len(prune_short=False)
                if self._root == 'best':
                    self.reroot(root=self._root)
                self.make_time_tree()
            self._advance_session_progress()

        # set root.gamma bc root doesn't have a branch_length_interpolator but gamma is needed
        if not hasattr(self.tree.root, 'gamma'):
            self.tree.root.gamma = 1.0

        # add coalescent prior
        if self._Tc is not None:
            from merger_models import coalescent
            self.logger('TreeTime.run: adding coalescent prior',2)
            coalescent(self.tree, Tc=self._Tc)
            self.make_time_tree()
            self._advance_session_progress()

        if self._relaxed_clock  and len(self._relaxed_clock)==2:
            import ipdb; ipdb.set_trace()
            slack, coupling = self._relaxed_clock
            # iteratively reconstruct ancestral sequences and re-infer
            # time tree to ensure convergence.
            niter = 0
            while niter<max_iter:
                # estimate a relaxed molecular clock
                self.relaxed_clock(slack=slack, coupling=coupling)

                ndiff = self.reconstruct_anc('ml')
                if ndiff==0:
                    break
                self.make_time_tree()
                niter+=1

            self._advance_session_progress()

        self._save_results()
        self._advance_session_progress()

    def _read_metadata_from_file(self, infile):
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

            else:
                print ("Metadata file has no column which looks like a sampling date!")

            metadata = df.to_dict(orient='index')
            return dates, metadata

        except:

            print ("Cannot read the metadata file. Exception caught")

    def _save_results(self):

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
            out_zip.write(os.path.join(self._root_dir, in_cfg), arcname=in_cfg)
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
                        print "cannot round:", node.__getattribute__(prop), "assigned as is"
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

        meta = []
        # if node has user-provided metadata, append it
        if node.name in self._metadata:
            node_meta = self._metadata[node.name]
            meta += [{'name':k, 'value': node_meta[k]} for k in node_meta]

        # append numdates to the metadata
        if hasattr(node, 'numdate'):
            meta.append({"name":"numdate", "value":node.numdate})

        # append deviation of the branch length from average
        meta.append({
                "name": "Branch length stretch",
                "value": (node.branch_length) / (node.mutation_length + 0.001)
            })

        if self._relaxed_clock:
            # append mutation rate deviation from average
            if hasattr(node, 'branch_length_interplator'):
                meta.append({
                    "name": "Relaxed mutation rate",
                    "value": node.branch_length_interpolator.gamma
                    })
            else:
                meta.append({"name": "Relaxed mutation rate", "value":1.0})

        return meta

    def _distribution_to_human_readable(self, dist):

        if dist.is_delta:
            date = utils.numeric_date() -  self.date2dist.get_date(dist.peak_pos)

            return [date-0.5, date-1e-10, date, date + 1e-10, date + 0.5],[0,0,1.0,0,0]

        peak_pos = dist.peak_pos
        fwhm = dist.fwhm
        raw_x = dist.x [(dist.x > peak_pos - 3*fwhm) & (dist.x < peak_pos + 3*fwhm)]
        dates_x = utils.numeric_date() -  np.array(map(self.date2dist.get_date, raw_x))
        y = dist.prob_relative(raw_x)
        return dates_x, y

    def _likelihoods_to_json(self):
        """
        Save the likelihoods fro the node positions in JSON format
        """
        out_dic = {}
        for node in self.tree.find_clades():
            x,y = self._distribution_to_human_readable(node.marginal_lh)
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

if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt

    root_dir = '../../sessions/YVGFQBWREHHH'
    myTree = TreeTimeWeb(root_dir, verbose=4)
    myTree.run()

