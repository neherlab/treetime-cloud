**********************************************************
Descirption of the dataset used in the TreeTime validation
**********************************************************

 This archive contains the data sets required to reproduce all TreeTime validation results.

 The data is supplementary to the article:

 To use the data located in this archiove, first download or clone the TreeTime-Validation project (python code) from

 https://pausag@gitlab.com/pausag/treetime_validation.git

 For Linux, open the treminal, navigate to the folder where to put the projectand type

 $git clone https://pausag@gitlab.com/pausag/treetime_validation.git.


 Then, unpack the archive parts in the root folder of the TreeTime-Validation project.

 For the detailed description of the data, refer to the DataDescrition.rst files  inside the data archives. The brief descriptionb is given below.

 The Archive consists of the following parts:

 * Simulated data.
 The validation of the treeTime package using the forward-time simulation for the evolution of a synthetic population. The evolving population was sampled regularly. From the sample data (sampling times and alignmets), the phylogenetic tree was reconstructed, which then used as the input for the TreeTime to reconstruct the time of the most-recent-common ancestor and similar parameters. The reconstructed data were compared to the simulation parameters. Reconstruction has been also made by other methods (BEAST, LSD) in order to compare the results to the TreeTime.

 * flu_H3N2
 The validation of the TreeTime package using the flu H3N2 sequences data (HA segment). The reconstruction for the mutation rates and time of the most-recent-common-ancestor were made by BEAST, LSD, TreeTime. The results were compared. The stability of the methods against the missing data and the sample size were also tested. The former were made by erasing the sampling dates from some of the leaves in the input tree. The latter - by sampling subtrees from the given phylogenetic tree such that all subtrees share the same root.

 * Skyline
 The estimation of the population size of the fluctuating population on the simulated data. (see also skyline tests on the ebola data in the validation project (python code). The evolution of the populsation was simulated using the FFpopSim package. The population size changed periodically. The population was sampled in a few points along the evolution timeline, and from these samples, the population size was reconstructed by the TreeTime.

 Besides that, there are resources and auxilliary files required to run the analysis from scratch.


