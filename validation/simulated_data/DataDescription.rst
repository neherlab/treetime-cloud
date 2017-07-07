***********************************************
Dataset for TreeTime validation: simulated data
***********************************************

This folder contains the whole dataset used in the TreeTime validation procedure using the data from the simulations for the forward-time population evolution. The folder contains also the results to compare TreeTime with similar software: BEAST and LSD. In the following, there is description of the simulated raw data and the output files of the LSD and BEAST.


Simulated data. Folder: dataset/
================================

The simulations were performed using the FFpopSim package. The population was initialized with genome of length L=10000, N=100 individual genomes. The evolution was simulated in the scope of Wright-Fischer model, no selection was applied. The mutation rate Mu was set individually for each simulation. For simplicity, the time scale was selected so that 1 generation in Wright-Fischer model corresponds to one year. Every Ts generations, the Nv=10 genomes are sampled. From these samples, the tree is constructed, and the sampling dates are also recorded and saved to the names of the tree nodes and sampled genomes. Therefore, each simulation consists of the constant following parameters: genome length (L), population size (N), sample volume (Nv),  number of sampling events (Ns); and the variable parameters: mutation rate(Mu), sampling frequency in generations (Ts). In order to assess the statistical error rates, we performed several simulations for each set of the input parameters. The results of the simulations are stored in the **dataset**folder. The files are named as follows:

FFpopSim_L10000_N100_Ns20_Ts<sampling_freq>_Nv10_Mu<mutation_rate>_<point_number>.<filetype>

where the point_number separates simulations with similar parameters from each other, and the <filetype> defines the type of the data produced. There are following files produced in the simulations:

  * <filetype> = nwk : simulated newick tree with branch lengths in generations

  * <filetype> = nuc.fasta : alignment reconstructed from the genomes sampled. Worth noting that the genotypes used in simulations are binary (0,1), which translated to the nucleotides by substituting 0 with 'A' and 1 with 'C'.

  * <filetype> = opt.nwk : The simulated newick tree with the branch lengths optimized using maximum-likelihood algorithm.

  * <filetype> = ft.nwk : The newick tree, reconstructed from the sampled alignment. The reconstruction has been performed with FastTree program, which implements NJ followed with maximum-likelihood algorithm.

  * <filetype> = treetime.nwk : The tree resulting form the TreeTime reconstruction, with internal nodes set according to the sampling dates constraints. The original tree for the treetime.nwk files are the .opt.nwk trees (true trees, with optimized branch lengths). These results can be used to assess the error rates introduced by the tree reconstruction. For validation procedure, these results were not used.

  * <filetype> = treetime.ft.nwk : The tree resulting form the TreeTime reconstruction, with internal nodes set according to the sampling dates constraints. The original tree for the treetime.ft.nwk files are the trees reconstructed with FastTree from the simulated alignment. Since this is the usual situation (alignment is known, the tree is not), these results only were used in the validation.

BEAST data. Folder: _beast/
===========================

 For each simulated tree and alignment, the BEAST reconstruction was performed. The BEAST resources and results are all stored in the _beast folder. For each BEATS simulation, we used the FastTree reconstructed tree (./dataset/<simulation_name>.ft.nwk files), the alignments (dataset/<simulation_name>.nuc.fasta). The dates of the internal nodes were extracted from their names.

 To automatize the BEAST simulation run, we use on template (file resources/beast/template_bedford_et_al_2015.xml) filling it with the tree, alignment, sampling dates and output log file. The resulting config is saved under <simulation_name>.config.xml name. The BEAST output is stored as <simulation_name>.log.txt, and the tree, sampled along BEAST simulations, are stored in the <simulation_name>.trees.txt files.
 Therefore, for each datapoint in the dataset/ folder, there are three files produced in the _beast/ folder:

  * <simulation_name>.config.xml
  * <simulation_name>.log.txt
  * <simulation_name>.trees.txt

 These files are enough to reproduce all BEAST results.

LSD data. Folder: _lsd/
=======================

 LSD is another program used to reconstruct the time of the most-recent-common-ancestor using the sampling dates of the tree leaves.

 To run LSD simulations, the file with sampling dates in a special format is needed. For each datapoint in dataset/ folder, the corresponding file is produced and stored in the

  * _lsd/<simulation_name>.lsd_dates.txt file.


 All data to reproduce and analyze LSD simulations are stored in the _lsd/ folder. For each datapoint with <simulation_name> in the dataset/ folder, the following files are produced:

  * _lsd/<simulation_name> (without extension) : the general results file, summarizing the input parameters used, and showing the results of the LSD reconstruction.

  * _lsd/<simulation_name>.date.newick : output newick tree with branch lengths in units of years.

  * _lsd/<simulation_name>.newick : output newick tree after the LSD reconstruction. Internal nodes positions set according to the LSD optimization algorithm. Branch lengths are in the units of the input tree branch lengths.

  * _lsd/<simulation_name>.nexus : output nexus tree after the LSD reconstruction. Internal nodes positions set according to the LSD optimization algorithm.


 For LSD simulations, both true trees and FastTree-reconstructed trees were used hence producing two types of files. Those using FastTree trees as input have "_fasttree" suffix in their names, and only those were used in the TreeTime validation and comparison procedure.


Results tables
==============

 To facilitate the further analysis of the results, we provide the CSV tables with the pre-processed results for TreeTime, LSD, BEAST reconstructions. These tables can be used to directly reproduce the TreeTime validation results.


TreeTime CSV table. File: _treetime_fasttree_res.csv
----------------------------------------------------

 Contains the following information:

  * File: the name of the simulation datapoint (as shown in the dataset/ folder)

  * Tmrca_real: Tmrca from the FFpopSim simulations

  * Tmrca: Tmrca, as reconstructed by TreeTime

  * Mu: Mutation rate as reconstructed by TreeTime (real mutation rate is encoded in the File name)

  * R^2(initial_clock): The regression coefficient of the leaf sampling dates vs root-to-tip distances. Used to assess the quality of the initial clock used for TreeTime reconstruction.

  * R^2(internal_nodes): The regression coefficient of the internal nodes dates vs node-to-tip distances in the reconstructed tree. It is used to assess the quality of the internal nodes arrangement after the reconstruction.

LSD CSV table. File: _lsd_fasttree_res.csv
-------------------------------------------

 Contains the following information:

  * File: the name of the simulation datapoint (as shown in the dataset/ folder)

  * Tmrca_real: Tmrca from the FFpopSim simulations

  * Tmrca: Tmrca, as reconstructed by LSD

  * Mu: Mutation rate as reconstructed by LSD (real mutation rate is encoded in the File name)

  * Objective: value of the objective function from the LSD optimization algorithm. NOTE: the latest versions of the LSD do not output the objective function values in the results file. In this case, it is set to 0.

BEAST CSV table. File: _beast_res.csv
-------------------------------------

 Contains the following information:

  * Filename: the name of the simulation datapoint (as shown in the dataset/ folder)

  * PopSize: population size decoded from the Filename

  * Tmrca_real: Tmrca from the FFpopSim simulations

  * ClockRate_real: Mutation rate used in FFpopSim simulations. Decoded from the Filename

  * SamplesNum: Number of samples taken in the FFpopSim simulations. Decoded from the Filename

  * SampleFreq: Sampling frequency in generations. Decoded from the Filename

  * TotEvoTime(SampleNum*SampleFreq): Total evolution time in generations

  * Nmu: PopSize * Mutation rate

  * LH: Tree Likelihood

  * LH_std: Standard deviation of the Tree Likelihood in a single BEAST run after the algorithm converged.

  * Tmrca: Reconstructed Tmrca

  * Tmrca_std: Standard deviation of the Tmrca in a single BEAST run after the  algorithm converged.

  * Mu: reconstructed mutation rate

  * Mu_std: Standard deviation of the mutation rate in a single BEAST run after the  algorithm converged.


Plotting the results
====================

 To plot the results, make sure first that the treetime_validation python project is installed and the simulated_data (this archive) is unpacked to the root folder of the project. For detailed instructions, see the manual in the root folder.


Tmrca, Mu
---------

 To plot the results of the Tmrca and mu reconstruction, run the script plot_simulated_data_tmrca_mu.py from the treetime_validation project.

 The output plots show the accuracy of the Tmrca and mutation rate  reconstruction in dependence of the mutation rate. (or, more precisely, N*mu product)

 Besides the mutation rate, there is another free parameter, used in the FFpopSim simulation, which is sampling frequency. This parameter controls the total tree depth T. Since the accuracy of the Tmrca reconstruction normally is within one coalescence time, we relate T to the population size (N) to get the tree depth in units of the coalescent time. The plot script is configured so that it shows the accuracy of the reconstruction for a single value of T/N ratio, as the reconstruction accuracy is different for trees of different depths. In the plot script, you can set a particular value of T/N ratio. In the default dataset the following possible ratios are defined:

 T/N = 2,4,10 (tree depth is from 2 to 10 coalescent times)


Accuracy of the internal nodes positions
----------------------------------------

 We also provide the script to show the accuracy of the internal node positions reconstruction.

 plot_simulated_data_bl_corr.py

 This script will parse trees produced by FastTree, BEAST, TreeTime, find similar splits, and plot the corresponding branch length related to the real branch length as simulated by FFpopSim. The script has no configuration. It only needs to access the output trees of the named methods.

