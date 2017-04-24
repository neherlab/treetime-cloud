# TreeTime validation project

This is the sumplementary project to the Treetime phylogeny package. If you are not yet familiar with the TreeTime itself, please read first about the main project ([GitHub page](https://github.com/neherlab/treetime)).

This project comprises the boilerplate code for Treetime validation, tests and benchmarking.

## Treetime validation

## Validation organization

Basically, the validation workflow is separated into two major parts: the dataset generation, and the processing the generated datasets followed by the results analysis and plotting. These two parts are intentionally implemented so that they can be run separately and independently.

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

In addition, if you are using custom tree and alignment, specify the path, where the script can find both:

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


### Analysis of the generated dataset



