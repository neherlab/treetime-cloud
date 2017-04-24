# TreeTime validation project
This is the sumplementary project to the Treetime phylogeny package. If you are not yet familiar with the TreeTime itself, please read first about the main project ([GitHub page](https://github.com/neherlab/treetime)). You should also have the TreeTime package installed to run the code from this project.

This project comprises the boilerplate code for Treetime validation, tests and benchmarking.

# Table of Contents

* [Prerequisites](#prerequisites)
* [Overview](#overview)
   * [Resources](#resources)
      * [External binaries](#external-binaries)
      * [Initial data](#initial-data)
* [Configuration and run](#configuration-and-run)
   * [Simulated data](#simulated-data)
   * [Influenza H3N2 - missing dates](#influenza-h3n2-reconstruction-with-missing-dates-information)
   * [Influenza H3N2 - subtrees](#influenza-h3n2-subtrees-of-a-single-big-tree)

# Prerequisites
To run the code, you need python-2.7 or later to be installed. You will also need `numpy, scipy, pandas, biopython` python libraries. To compare the Treetime against other phylogenetic packages ([LSD](http://www.atgc-montpellier.fr/LSD/), and [BEAST](http://beast.bio.ed.ac.uk/)), you need them to be installed in your system (refer the [External binaries](#external-binaries) section for more details). To generate dataset, we also use [FastTree](http://www.microbesonline.org/fasttree/) and [FFpopSim](http://webdav.tuebingen.mpg.de/ffpopsim/). The latter requires compilation, so if you decide to generate the whole datasets yourselves, you will need the compilation tools: `g++-4.8` or later, `gsl`, `boost`. See detailed instructions in the [External binaries](#external-binaries) section.

# Overview
Basically, the validation workflow consists of the three independent parts:

 - Validation of the TreeTime results on the simulated data.
 - Validation of the results on the Influenza data (sampling subtrees).
 - Validation of the results for the Influenza data (with only fraction of leaf dates known).

Each part can be run independently. The code for the parts are separated into the files with common suffixes: `XXX_submit.py` and `XXX_run.py`. The files of first type create the parameter spaces for the simulations and for each set of the parameters, call the `XXX_run.py` files, which perform simulation for a given set of parameters.

All common functions and variables are defined in files `utility_functions_XXX.py`.

To plot the resutls, use the `plot_xxx_res.py` files. These files import the default plotting parameters from the `plot_defaults.py` to make all figures have the same style and colors.

Since the project relies on many external binaries, for convenience, they all registered in the `external_binaries.py` file.

## Resources

### External binaries
The comparison to the other packages requires some external binaries to be downloaded, installed to the host system, and registered by the project. To avoid the issues with the paths of the program installation, there is an `external_binaries.py` file, which lists all third-party binaries used by the project.

The two phylogeny packages used in the validation are the [LSD](http://www.atgc-montpellier.fr/LSD/), and [BEAST](http://beast.bio.ed.ac.uk/). Download them and install, using the instructions provided by the vendors. After all done, the treetime-validation project should be able to find the executables of the projects.  In addition, for tree generation, we use the [FastTree](http://www.microbesonline.org/fasttree/) package. Download and install following the instructions.

Add the paths to the phylogeny packages to the  `external_binaries.py` file:

```python
FAST_TREE_BIN = "/usr/bin/fasttree"
LSD_BIN  = "/usr/bin/lsd"
BEAST_BIN = "/opt/BEASTv1.8.4/lib/beast.jar"
```

The simulation of the evolution process also includes usage of the [FFpopSim](http://webdav.tuebingen.mpg.de/ffpopsim/) forward-time simulation library.

NOTE: to enable intermediate sampling in the population in the evolution process, we use the FFpopSim extension.

Download the FFpopSim from [GiHub page](https://github.com/neherlab/ffpopsim), and checkout to the branch `historical_samples`. Compile the library using the instructions provided with the code. Then, compile the cpp programs developed for the treetime validation. The cpp files are located in the `resources/src` directory. The compilation is performed as shown below:

NOTE: To compile FFpopsim, you should have gnu-scientific library installed

```bash
$g++ -o ffpopsim --std=c++11 FFpopSim.cpp -I <path-to-ffpopsim-headers> [-L <path-to-ffpopsim-lib>] -lFFPopSim -lgsl -lgslcblas
```

There is another version used to generate data for skyline validation:

```bash
$g++ -o ffpopsim_skyline --std=c++11 FFpopSim_skyline.cpp -I <path-to-ffpopsim-headers> [-L <path-to-ffpopsim-lib>] -lFFPopSim -lgsl -lgslcblas
```

Then, register the compiled binaries to the `external_binaries.py` file:

```python

FFPOPSIM_BIN  = "<path-to-ffpopsim binary>"
FFPOPSIM_SKYLINE_BIN = "<path-to-ffpopsim_skyline binary>"
```

### Initial data
The data usde for the simulation are located in the `resources` folder. Basically, the validation scripts need the Influenza tree and alignment with substantial number of leaves (ideally, around 5000 leaves). We have chosen the influenza sequences for the period 2011-2013 years, segment HA. The sequences were downloaded from [Influenza Research Database](https://www.fludb.org), re-aligned. The tree for the alignment has been built using FastTree. In principle, any other tree can be used for the alignment. In case another tree is used, make sure you specify the name of the tree accordingly where needed.

To run Beast, we need som configuration template. As we work with influenza trees, we use the previously tested and published configuration from the [Bedford, Russell, et al](http://www.nature.com/nature/journal/v523/n7559/full/nature14460.html?WT.ec_id=NATURE-20150709&spMailingID=49054266&spUserID=MjA1NjIyNTk4MQS2&spJobID=720986702&spReportId=NzIwOTg2NzAyS0). The file is located in the `resources` folder. Note that we use it only as a configuration template. Than means, in each simulation we specify the actual sequences, initial tree and leaf dates, whereas all the configurations stay the same accross simulations.

# Configuration and run


## Simulated data
To validate the TreeTime results, we use the trees, simulated by the FFPopSim forward-time evolution simulation package. Before run the analysis in this section, please make sure you have downloadad and compiled FFPpopSim as described above. The validation is separated into two parts - the generation of the simulated dataset + its analysis, and plotiing the results.
The results of the TreeTime are compared against the Beast and LSD packages. If you use them in the analysis, please make sure you have donloaded the binaries, and configured the paths accordingly. The dataset generation and analysis is made by the two python scripts. First, configure the `generate_simulated_dataset_submit.py`.

#### Single-point simulations (Run script)
This script runs single-point simulations for a given set of parameters. It should not require any configuration except enabling/disabling the simulation steps.

```python

    # if true, everything will be regenerated from scratch. otherise, the
    # existing dataset will be used
    GENERATE_SIMULATED_DATA = True

    RUN_TREETIME = True
    RUN_LSD = True
    RUN_BEAST = True
```
#### Whole dataset generation (Submit script)
This script creates the range of the parameters used and then for each set of the input parameters calls `generate_simulated_dataset_run.py` script. First, define the output directories and filenames for the generated data:

```python
    # Directory to store results
    res_dir = "./simulated_data/dataset"
    # File prefix to store the formatted output for further processing/comparison
    outfile = "./simulated_data/2017-04-19"
```
Then, you should provide the range of the parameters to simulate:

```python
    # FFPopSim simulation parameters
    L = 1e4  # sequence length
    N = 100  # population size
    SAMPLE_VOL = 10  # number of sequences sampled to be included in the tree
    SAMPLE_NUM = 20  # number of samples (defines the total tree length)
    SAMPLE_FREQS = [10, 20, 50]  # number of generations between sampling.
                                 # Defines the total tree length.
                                 # In this case, tot evo time = T/N =
                                 # [2, 4, 10] coalescence times

    # mutation rates used in the simulations
    MUS = [7e-6, 1e-5, 2e-5, 5e-5, 7e-5, 1e-4, 2e-4, 5e-4, 7e-4, 1e-3, 2e-3]
    # number of simulations per set of parameters
    N_POINTS = 20

    # staring repetiion number. If you have some dataset, setting this number
    # to a bigger value will cause the sets concatenation rather then overwrite
    N_0 = 0
```

When all done, you should decide, whether the simulations will be run remotely on the cluster, or locally on your computer. In the first case, you should set the cluster flag to ON and edit the call command. Below, there is  a configuration for the Sun Grid Engine cluster.

```python

    CLUSTER = True

    if CLUSTER:
        call = ['qsub', '-cwd', '-b','y',
               '-l', 'h_rt=23:59:0', # BEAST might run long
                #'-o', './stdout.txt',
                  #'-e', './stderr.txt',
                '-l', 'h_vmem=50G', # BEAST requires A LOT
                 './generate_simulated_dataset_run.py']
    else:
        call = ['./generate_simulated_dataset_run.py']
```

NOTE: the time and the amount of memory is set for the Beast run, which uses Java virtual machine and requires a lot of memory just to be started. If you are not using Beast, you could decrease the runtime to 20 mins, and the amount of memory to 3G.


## Influenza H3N2 - reconstruction with missing dates information
The main goal of this part is to analyze the influence of the missing temopral information on the precision of the Tmrca reconstruction. The validation is sseparated into two parts: the dataset generation+analysis, and the results plotting.

### Dataset generation
To generate dataset and to run the simulations on the generated data, configure and run the two python scripts: `generate_flu_missingDates_submit.py` and `generate_flu_missingDates_run.py`. The first script generates the range of the parameters to be used in the simulations, and calls the second script for each combination of the input parameters. The seond script performs datapreparation for a single set of parameters, and then runs the treetime simulations for the specified parameters.

#### Single-point simulations (Run script)
The script performs rnadom sampling of the leaves on a given tree, erases the temporal information for these leaves, and then perform the treetime simulations. normally, it requeres no configuration.

#### Whole dataset generation (Submit script)
Before running the scipt, it should be configured. First of all, specify the output directories where the results will be placed.

```python
    #output directories
    out_dir = "./flu_H3N2/missing_dates/"
    subtree_dir = os.path.join(out_dir, "subtrees")
```

Then, define the name format of the output files:

```python
    # file formats
    resfile_fmt  = os.path.join(out_dir, "./H3N2_HA_2011_2013_{}seqs_res.csv")
    resfile_dates_fmt  = os.path.join(out_dir, "./H3N2_HA_2011_2013_{}seqs_dates_res.csv")
    treefile_fmt = os.path.join(subtree_dir, "./H3N2_HA_2011_2013_{}seqs.nwk")
    alnfile_fmt  = os.path.join(subtree_dir, "./H3N2_HA_2011_2013_{}seqs.fasta")
```

In the end, you should specify the range of the used parameters:

```python
    # for the simulations, use the tree with the specified number of
    # sequences. The trees should be produced beforehand and placed in the s
    # subtrees folder defined above.
    nseqs = [100]
    # Set of the fraction of dates known.
    dates_knonwn_fraction = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                             0.6, 0.7, 0.8, 0.9, 1.0]
    # Numbe of repetitions for each subtree and each known fraction. Multiple
    # repetitions allow to estimate the statistical error. The good estimate
    # for this parameter is between 10 and 50.
    Npoints = 20
```

In the end, you should decide, whether the simulations will be run remotely on the cluster, or locally on your computer. In the first case, you should set the cluster flag to ON and edit the call command. Below, there is  a configuration for the Sun Grid Engine cluster.

```python
    CLUSTER = True
    call = ['qsub', '-cwd', '-b','y',
                             '-l', 'h_rt=1:59:0',
                             #'-o', './stdout.txt',
                             #'-e', './stderr.txt',
                             '-l', 'h_vmem=3G',
                             './generate_flu_missingDates_dataset_run.py']
```


After all done, run the dataset generation running

```bash
$python generate_flu_missingDates_dataset_run.py
```
### Plotting the results
To plot the results, first you shold configure the `plot_flu_missing_dates.py` script. Normally, you should only specify the filenames, where the results of the above simulations are stored. It is only necessary, if the filenames or paths were modified. optionally, set flag to save figures.

```python
    save_fig = True

    work_dir = './flu_H3N2/missing_dates'
    fname_format = 'H3N2_HA_2011_2013_{}seqs_res.csv'
    fname_dates_format = 'H3N2_HA_2011_2013_{}seqs_res.csv_dates.csv'
```

Then, run the script:

```bash
$python ./plot_flu_missing_dates.py
```

## Influenza H3N2 - subtrees of a single big tree
This part describes the datasetgeneration and processing for the subtrees of a single Influenza tree. The main purpose of this validation run is to prove the stability of the TreeTime inferrence on the sample size. The results are also compared against the two best competitors - Beast and LSD. To show the stability of the inferred Tmrca on the sample size, we perform sampling of subtrees from the bigger Influenza tree. The latter is store in the resources folder in the repository. The sampling is done so that the root of every sampled tree is the same as of the initial tree, which implies that the expected inferred Tmrca date should be the same on every suubtree.

### Dataset generation
The dataset generation is performed by the two scripts: one is to run the simulations for a given set of parameters (`./generate_flu_subtrees_dataset_run.py`). Another  script (`./generate_flu_subtrees_dataset_submit.py`)  creates the range of the parameters, and for each combination of the parameters, calls the 'run' script. It also decides whether to run the simulations directly on the local computer, or it should be submitted to a remote cluster.

#### Single point simulation (Run script)
The configuration of the run script is basically to toggle switches to set which simulations should be performed:

```python
RUN_TREETIME = True
RUN_LSD = True
RUN_BEAST = True
```

Please make sure you have included the right binaries in the `external_bins.py` file. In addition, if you are using custom tree and alignment, specify the path, where the script can find both:

```python
aln_name = "./resources/flu_H3N2/H3N2_HA_2011_2013.fasta"
tree_name = "./resources/flu_H3N2/H3N2_HA_2011_2013.nwk"
```

NOTE: tree should be in newick format, the alignment should be in the fasta format.

#### Whole dataset generation (Submit script)
To run the validation, first edit the `generate_flu_subtrees_dataset_submit.py` file to adapt to your run conditions. The example of the configuration is shown below:


First, configure the output directory and filenames to store the results:

```python
# root dir to store all results
work_dir = "./flu_H3N2/subtree_samples"
#sub-directories and file names
out_dir = os.path.join(work_dir, "2017-04-20")
treetime_res_file = os.path.join(work_dir, "2017-04-20_treetime_res.csv")
lsd_res_file = os.path.join(work_dir, "2017-04-20_lsd_res.csv")
```

If the LSD simulation will be enabled, specify the parameters for the LSD run
(only in case they will differ from the defaults). The default parameters are set as shown below:

```python
# LSD run configuration
lsd_parameters = ['-c', '-r', 'a', '-v']
```

In the end, configure the simulation parameters. Specify the range of the subtree sample sizes:

```python
N_leaves_array = [20, 50, 100, 200, 500, 750, 1000, 1250, 1500, 1750, 2000]
```

and give the number of the repetitions for every subtree (each subtree size will be sampled tese number of times). The recommended `n_iter` size is between 10 and 50 to get enough statistics. This will give us good estimation for the statistical error of the simulation percision.

```python
n_iter = 20
```

Finally, decide whether you will run the simulations on a cluster in parallel, or on a local computer. In the former case, you should set `CLUSTER=True` and configure the cluster submit command. The example below shows the configuration for the Sun Grig Engine cluster.

```python

        ...

        if CLUSTER:
            call = ['qsub', '-cwd', '-b','y',
                    '-l', 'h_rt=23:59:0',
                    #'-o', './stdout.txt',
                    #'-e', './stderr.txt',
                    '-l', 'h_vmem=50G',
                    './generate_flu_subtrees_dataset_run.py']
        else:
            call = ['./generate_flu_subtrees_dataset_run.py']
        ...
```

NOTE: the time and the amount of memory is set for the Beast run, which uses Java virtual machine and requires a lot of memory just to be started. If you are not using Beast, you could decrease the runtime to 20 mins, and the amount of memory to 3G.

### Plotting the results
To plot the results, first edit the `./plot_flu_subtrees_res.py` file. The configuration includes specifying the filenames, where the results should be found:

```python
    #  directory to search for the result tables:
    res_dir = './flu_H3N2/subtree_samples/'
    treetime_res = os.path.join(res_dir, '2017-04-20_treetime_res.csv')
    lsd_res = os.path.join(res_dir, '2017-04-20_lsd_res.csv')
    beast_log_dir = os.path.join(res_dir, '2017-04-20/beast_out')
    beast_tree_dir = os.path.join(res_dir, '2017-04-20/subtrees')
```

And turning on and off the switches to set the data, whoch should be plotted:

```python
    PLOT_TREETIME = True
    PLOT_LSD = True
    PLOT_BEAST = True
```

After all done, just run the script:

```bash
$python plot_flu_subtrees_res.py
```




