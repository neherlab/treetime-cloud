***********************************************
Dataset for TreeTime validation: flu H3N2, HA
***********************************************

 This folder contains the whole dataset used in the TreeTime validation procedure with the Flu H3N2 (HA segment) dataset. The data (multiple sequence alignment) was downloaded from the flu database, and realigned. The data were restricted to those sequences, which have sampling dates specified between 2011 and 2013 with month precision. In total, there were 3585 sequences in the alignment. The sequences are in the resources/flu_H3N2 folder (attached as a separate archive). The tree for the alignment was reconstructed using FastTree program.

 The flu validation consists of the two parts: (i) estimating the stability of the TreeTime reconstruction on the tree size and (ii) estimating of the treetime sensitivity to the partial data missing.


**********************************************************************
Validation for the trees of different sizes. Folder: subtree_samples/
**********************************************************************

 In order to probe the TreeTime, LSD, BEAST stability on the sizes of the tree, the following procedure has been developed. From the original Flu tree, containing 3585 leaves, several trees of different number of leaves were sampled. The samples were made so that (i) the root of the original tree is the same as of the subtrees and (ii) the number of samples per year sustained uniform where possible. The sampled sequences were chosen randomly. To estimate the statistical error rates, subtree construction procedure has been repeated 10 times hence resulting in to independent reconstructions.

 Since all subtrees share the same root, the reconstruction by all methods should result in the same value of the root node date and hence the same mutation rate for all subtrees.

 The data consists of the following parts, located in the different folders:

 dataset/subtrees/ : subtree samples from the initial tree

 dataset/beast_out/ : BEAST configuration and output files to reproduce all presented BEAST results

 dataset/LSD_out/ : LSD input dates file and results.


Subtree samples. Folder: /dataset/subtrees
==========================================

 This folder contains subtrees of the initial tree. The name of the files in this folder are:

 H3N2_HA_2011_2013_<Nleaves>_<sample_num>.nwk

 The filename prefix indicates the name of the initial tree, the <Nleaves> is the number of leaves in the subtree, the <sample_num> is the number of the subtree with the particular number of leaves.
 In the following, the string  H3N2_HA_2011_2013_<Nleaves>_<sample_num> is referred to as subtree_name

BEAST results. Folder: /dataset/beast_out
=========================================

 For each subtree, we have performed BEAST reconstruction. The log file used as BEAST input was created from the template, located in resources/beast/template_bedford_et_al_2015.xml (attached as a separate archive). This template was filled with the particular subtree. The dates of the leaves were extracted from the leaf names, the relevant sequences were taken from the initial alignment. The Beast configuration file is located in:

 /dataset/beast_out/<subtree_name>.config.xml

 The BEAST output includes the log file and the trees sampled along the BEAST computations. These are located in the

 /dataset/beast_out/<subtree_name>.log.txt
 /dataset/beast_out/<subtree_name>.trees.txt

 respectively.

LSD data. Folder: /dataset/LSD_out
==================================

 This folder contains the files required to perform LSD computations on each subtree as well as the results of the LSD reconstruction.

 Input files for LSD are the subtree file (from the subtrees folder) and the dates file in a special format. The latter is accessible in the LSD_out folder under the name:

 /dataset/LSD_out/<subtree_name>.lsd_dates.txt

 The results of the LSD reconstruction are stored in the following files:

 /dataset/LSD_out/<subtree_name> (without extension) : the general results file, summarizing the input parameters used, and showing the results of the LSD reconstruction.

 /dataset/LSD_out/<subtree_name>.date.newick : output newick tree with branch lengths in units of years.

 /dataset/LSD_out/<subtree_name>.newick : output newick tree after the LSD reconstruction. Internal nodes positions set according to the LSD optimization algorithm. Branch lengths are in the units of the input tree branch lengths.

 /dataset/LSD_out/<subtree_name>.nexus : output nexus tree after the LSD reconstruction. Internal nodes positions set according to the LSD optimization algorithm.


Results tables
==============

 To facilitate the further analysis of the results, we provide the CSV tables with the pre-processed results for TreeTime, LSD, BEAST reconstructions. These tables can be used to directly reproduce the TreeTime validation results.


TreeTime CSV table. File: treetime_res.csv
------------------------------------------

 Contains the following information:

  * File: the name of the simulation datapoint (as shown in the dataset/ folder)

  * N_leaves: number of leaves in the subtree

  * Tmrca: Tmrca, as reconstructed by TreeTime

  * Mu: Mutation rate as reconstructed by TreeTime

  * R^2(initial_clock): The regression coefficient of the leaf sampling dates vs root-to-tip distances. Used to assess the quality of the initial clock used for TreeTime reconstruction.

  * R^2(internal_nodes): The regression coefficient of the internal nodes dates vs node-to-tip distances in the reconstructed tree. It is used to assess the quality of the internal nodes arrangement after the reconstruction.

  * Runtime: TreeTime total runtime, in seconds


LSD CSV table. File: _lsd_fasttree_res.csv
-------------------------------------------

 Contains the following information:

  * File: the name of the simulation datapoint (as shown in the dataset/ folder)

  * Nleaves: number of leaves in the subtree

  * Tmrca: Tmrca, as reconstructed by LSD

  * Mu: Mutation rate as reconstructed by LSD

  * Runtime: LSD total runtime, in seconds

  * Objective: value of the objective function from the LSD optimization algorithm. NOTE: the latest versions of the LSD do not output the objective function values in the results file. In this case, it is set to 0.


BEAST CSV table. File: beast_res.csv
------------------------------------

 Contains the following information:

  * Filename: the name of the simulation datapoint (as shown in the dataset/ folder)

  * Nleaves: number of leaves in the subtree

  * LH: Tree Likelihood

  * LH_std: Standard deviation of the Tree Likelihood in a single BEAST run after the algorithm converged.

  * Tmrca: Reconstructed Tmrca

  * Tmrca_std: Standard deviation of the Tmrca in a single BEAST run after the  algorithm converged.

  * Mu: reconstructed mutation rate

  * Mu_std: Standard deviation of the mutation rate in a single BEAST run after the  algorithm converged.


Plotting the results
====================

 To plot the results, make sure first that the treetime_validation python project is installed and the flu_H3N2 (this archive) is unpacked to the root folder of the project. For detailed instructions, see the manual in the upper directory.

 To plot the results, run the plot_flu_subtrees_res.py script, no further configuration required


**********************************
Validations with the missing data
**********************************

 To assess the TreeTime stability to the missing data, the following procedure was used: Smaller subtree, containing 100 nodes has been created from the initial flu tree. In this subtree, a fraction of node's dates were erased, so the TreeTime reconstruction run with incomplete data. Then, the two types of the validations were made: (i) the stability of the Tmrca reconstruction vs fraction of missing leaf dates and (ii) the precision of the unknown  leaf dates reconstruction

 The procedure is as follows: from a given tree, we choose randomly the nodes, which dates are known (the fraction of know dates is given). After that, the TreeTime and Beast run the reconstruction of Tmrca knowing the dates of only the fraction of nodes. Each reconstruction (incl. dates erasure) is repeated 10 times to assess the statistical error of the reconstruction.

 The data is located in the /missing_dates/ folder. The /missing_dates/subtrees/ folder contains the subtrees used for the validation.

 The /missing_dates/beast_out folder contains all data required to run BEAST reconstruction. The data files are named as follows:

 /missing_dates/beast_out/H3N2_HA_2011_2013_100seqs_Nk<known_dates_fraction>_<run_number>.config.xml

 /missing_dates/beast_out/H3N2_HA_2011_2013_100seqs_Nk<known_dates_fraction>_<run_number>.log.txt

 where the run_number is the number of the simulation.


Results tables
==============

The results tables are of two types: the reconstruction of Tmrca (with BEAST and TreeTime), and the reconstruction of the leaf dates.


Reconstruction of Tmrca, TreeTime. File: H3N2_HA_2011_2013_100seqs_res.csv
---------------------------------------------------------------------------


 This file contains the following information:

 * Filename: subtree filename

 * KnownDatesFraction: The fraction of dates known

 * Tmrca: the date of the most-recent common ancestor reconstructed by TreeTime

 * Mu: mutation rate reconstructed by TreeTime

 * R^2(initial_clock): The regression coefficient of the leaf sampling dates vs root-to-tip distances. Used to assess the quality of the initial clock used for TreeTime reconstruction.

 * R^2(internal_nodes): The regression coefficient of the internal nodes dates vs node-to-tip distances in the reconstructed tree. It is used to assess the quality of the internal nodes arrangement after the reconstruction.

 * Runtime(sec): TreeTime total runtime, in seconds


Reconstruction of Tmrca, BEAST. File: H3N2_HA_2011_2013_100seqs_beast_res.csv
------------------------------------------------------------------------------

 This file contains the following information:

 * Filename: subtree filename

 * KnownDatesFraction: The fraction of dates known

 * LH: Tree Likelihood

 * LH_std: Standard deviation of the Tree Likelihood in a single BEAST run after the algorithm converged.

 * Tmrca: Reconstructed Tmrca

 * Tmrca_std: Standard deviation of the Tmrca in a single BEAST run after the  algorithm converged.

 * Mu: reconstructed mutation rate

 * Mu_std: Standard deviation of the mutation rate in a single BEAST run after the  algorithm converged.


Reconstruction of the unknown leaf dates. File: H3N2_HA_2011_2013_100seqs_beast_res.csv
-----------------------------------------------------------------------------------------------------------------------

This file contains the following information:

 * LeafName: name of the leaf

 * Known_dates_fraction: fraction of the dates known

 * Tmrca: date of the most-recent common ancestor (reconstructed with TreeTime)

 * LeadDate_real: real date of the leaf

 * LeafDate_rec: reconstructed date of the leaf

 * DateError: error in the date reconstruction (real-rec)


Plotting the results:
----------------------

 To plot the results, make sure first that the treetime_validation python project is installed and the flu_H3N2 (this archive) is unpacked to the root folder of the project. For detailed instructions, see the manual in the upper directory.

 To plot the Tmrca reconstruction results, run the plot_flu_missing_dates_res.py script, no further configuration required

 To plot the unknown leaf dates reconstruction results, run the plot_flu_missing_dates_leafDatesReconstruction.py script, no further configuration required
