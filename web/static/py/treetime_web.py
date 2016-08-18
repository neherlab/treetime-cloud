from treetime.treetime import TreeTime
import os, json

available_gtrs = ['Jukes-Cantor']

# file names, uniform across the sessions
session_state_file = 'session_state.json'
in_tree = "in_tree.nwk"
in_aln = "in_aln.fasta"
in_meta = "in_meta.csv"
in_cfg = "config.json"


class TreeTimeWeb(TreeTime):
    """
    TreeTime class extension to be used with web interface.
    This class defines functionality to communicate with the server, provide an
    arbitrary data to display

    """

    def __init__(self, root_dir, tree_builder, *args,**kwargs):

        # set the default configurations
        self._build_tree = False
        self._infer_gtr = True
        self._root = None
        self._mutation_rate = None
        self._Tc = None

        # set the directory, containing files for the session
        self._root_dir = root_dir
        self._nwk =  os.path.join(self._root_dir,in_tree)
        self._aln =  os.path.join(self._root_dir,in_aln)
        self._meta = os.path.join(self._root_dir,in_meta)
        self._cfg  = os.path.join(self._root_dir,in_cfg)
        #  read the metadata
        dates, metadata = self._read_metadata_from_file(self._meta)

        # read the JSON configuration file
        with open(self._cfg) as ff:
            config_dic = json.load(ff)

        if 'gtr' in config_dic:
            if config_dic['gtr'] == 'infer':
                self._infer_gtr = True
                gtr = "Jukes-Cantor"
            elif config_dic['gtr'] in available_gtrs:
                self._infer_gtr = False
                gtr = config_dic['gtr']
            else:
                # NOTE cannot use the logger class before the super.__init__ is called
                # therefore, just use plain print
                print ("Configuration contains unknown GTR type. Ignoring, will use Jukes-Cantor.")
                self._infer_gtr = False
                gtr = "Jukes-Cantor"
        else:
            print ("Configuration has no gtr model set. Will use default (Jukes-Cantor)")
            self._infer_gtr = False
            gtr = "Jukes-Cantor"

        super(TreeTime, self).__init__(dates=dates, tree=self._nwk, aln=self._aln, gtr=gtr, *args, **kwargs)
        self._init_session_state(config_dic)
        self._metadata = metadata

    def _init_session_state(self, config_dic):

        session_state = {
            "state" : "running",
            "error" : "", #  no errors so far
            "todo" : [],
            "done" : [],
            "progress":""
            }

        if 'build_tree' in config_dic and config_dic['build_tree']:
            session_state['todo'].append('Build phylogenetic tree')
            self._build_tree = True


        if self._infer_gtr: # the GTR configuration has been read in the __init__
            session_state['todo'].append('Reconstruct ancestral states, infer evolutionary model')
        else:
            session_state['todo'].append('Reconstruct ancestral states')


        if 'reroot' in config_dic and 'reroot' is not None:
            session_state['todo'].append('Find better root for the input tree')
            self._root = config_dic['reroot']


        if 'mu' in config_dic and config_dic['mu']!=0 and config_dic['mu'] is not None:
            try:
                self._mutation_rate=float(config_dic['mu'])
            except:
                self.logger("Cannot parse mutation rate from the configuration dictionary!", 1, warn=True)
                self._mutation_rate = None

            if self._mutation_rate is not None and self._root == 'best' :
                self.logger("Warning! the mutation rate cannot be used together "
                    "with the tree rooting to optimize molecular clock! Setting mutation rate to None", 1, warn=True)
                self._mutation_rate = None

        # in any case, we set the make time tree as todo
        session_state['todo'].append('Make time tree')

        if 'resolve_poly' in config_dic and config_dic['resolve_poly']:
            self._resolve_polytomies = True
            session_state['todo'].append("Resolve multiple mergers")
        else:
            self._resolve_polytomies = False


        if 'coalescent' in config_dic:
            if config_dic['coalescent'] == 0:
                self._Tc = None
            else:
                try:
                    self._Tc = float(config_dic['coalescent'])
                    session_state['todo'].append("Model coalescent process")
                except:
                    self.logger("Cannot parse coalescent timescale from the config dictionary!", 1, warn=True)
                    self._Tc = None


        if 'relax_clock' in config_dic:
            try:
                rc = config_dic['relax_clock']
                slack, coupling = float(rc['slack']), float(rc['coupling'])
                if slack == 0. and coupling == 0.:
                    self._relaxed_clock = None
                    session_state['todo'].append('Relax molecular clock')
                else:
                    self._relaxed_clock = (slack, coupling)
            except:
                self.logger("Cannot parse relaxed clock parameters from the config dict!", 1, warn=True)
                self._relaxed_clock = None

        self._session_state = session_state

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

        if self._build_tree:
            self._advance_session_progress()
            tree_builder(self._in_aln)

        self._advance_session_progress()
        self.optimize_seq_and_branch_len(infer_gtr=self._infer_gtr, prune_short=True)

        if self._root is not None:
            self._advance_session_progress()
            self.reroot(root=self._root)

        self._advance_session_progress()
        self.make_time_tree(slope=self._mutation_rate)

        if self._resolve_polytomies:
            self._advance_session_progress()
            # if polytomies are found, rerun the entire procedure
            if self.resolve_polytomies():
                self.prepare_tree()
                self.optimize_seq_and_branch_len(prune_short=False)
                if self._root == 'best':
                    self.reroot(root=self._root)
                self.make_time_tree()

        # set root.gamma bc root doesn't have a branch_length_interpolator but gamma is needed
        if not hasattr(self.tree.root, 'gamma'):
            self.tree.root.gamma = 1.0

        # add coalescent prior
        if self._Tc is not None:
            self._advance_session_progress()
            from merger_models import coalescent
            self.logger('TreeTime.run: adding coalescent prior',2)
            coalescent(self.tree, Tc=self._Tc)
            self.make_time_tree()

        if self._relaxed_clock  and len(self._relaxed_clock)==2:
            self._advance_session_progress()
            slack, coupling = self._relaxed_clock
            # iteratively reconstruct ancestral sequences and re-infer
            # time tree to ensure convergence.
            niter = 0
            while niter<max_iter:
                if self._relaxed_clock:
                    # estimate a relaxed molecular clock
                    self.relaxed_clock(slack=slack, coupling=coupling)

                ndiff = self.reconstruct_anc('ml')
                if ndiff==0:
                    break
                self.make_time_tree()
                niter+=1

    def _read_metadata_from_file(self, infile):
        """
        @brief      Reads a metadata from file or handle.

        @param      self    The object
        @param      infile  The input file or file handle

        """

        try:

            import pandas
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


if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt

    root_dir = '../../sessions/AFZGTFAIIOGI'
    myTree = TreeTimeWeb(root_dir, lambda x: None, verbose=4)
    myTree.run()

