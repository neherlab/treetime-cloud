# TreeTime validation project
This is the sumplementary project to the Treetime phylogeny package. If you are not yet familiar with the TreeTime itself, please read first about the main project ([GitHub page](https://github.com/neherlab/treetime)).

This project comprises the boilerplate code for Treetime validation, tests and benchmarking.

## Organization
Basically, the validation workflow is separated into two major parts: the dataset generation, and the processing the generated datasets followed by the results analysis and plotting. These two parts are intentionally implemented so that they can be run separately and independently.

## Configuration and run

### Influenza H3N2 - reconstruction with missing dates information
The main goal of this part is to analyze the influence of the missing temopral information on the precision of the Tmrca reconstruction. The validation is sseparated into two parts: the dataset generation+analysis, and the results plotting.

#### Dataset generation
To generate dataset and to run the simulations on the generated data, configure and run the two python scripts: `generate_flu_missingDates_submit.py` and `generate_flu_missingDates_run.py`. The first script generates the range of the parameters to be used in the simulations, and calls the second script for each combination of the input parameters. The seond script performs datapreparation for a single set of parameters, and then runs the treetime simulations for the specified parameters.

##### Single-point simulations (Run script)
The script performs rnadom sampling of the leaves on a given tree, erases the temporal information for these leaves, and then perform the treetime simulations. normally, it requeres no configuration.

##### Whole dataset generation (Submit script)
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
#### Plotting the results
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

### Influenza H3N2 - subtrees of a single big tree
This part describes the datasetgeneration and processing for the subtrees of a single Influenza tree. The main purpose of this validation run is to prove the stability of the TreeTime inferrence on the sample size. The results are also compared against the two best competitors - Beast and LSD. To show the stability of the inferred Tmrca on the sample size, we perform sampling of subtrees from the bigger Influenza tree. The latter is store in the resources folder in the repository. The sampling is done so that the root of every sampled tree is the same as of the initial tree, which implies that the expected inferred Tmrca date should be the same on every suubtree.

#### Dataset generation
The dataset generation is performed by the two scripts: one is to run the simulations for a given set of parameters (`./generate_flu_subtrees_dataset_run.py`). Another  script (`./generate_flu_subtrees_dataset_submit.py`)  creates the range of the parameters, and for each combination of the parameters, calls the 'run' script. It also decides whether to run the simulations directly on the local computer, or it should be submitted to a remote cluster.

##### Single point simulation (Run script)
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

##### Whole dataset generation (Submit script)
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

#### Plotting the results
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




