*******************************************************
Dataset for TreeTime validation: skyline reconstruction
*******************************************************

 This data set contains data required to reproduce the TreeTime results for the population size reconstruction. To perform the validation, a synthetic population has been evolved by the FFPopSim. The population size was modulated so that it changed sinusoidally with given period and amplitude. During the evolution, the population was sampled and then reconstructed by the TreeTime.

 The archive contains the simulation results and the reconstruction made by the TreeTime.


Simulation results (./simulated_data folder)
============================================

   The simulations have been made for different set of parameters. For each simulation point, there are three files produced:

   * <simulation_name>.nwk : phylogenetic tree, with branch lengths in units of population generations
   * <simulation_name>.opt.nwk : phylogenetic tree with branch lengths set to the maximum-likelihood values
   * <simulation_name>.nuc.fasta : multiple sequence alignment formed from sequences sampled during the evolution


    Each <simulation_name> encodes the critical simulation parameters to reproduce the simulations.

    Simulation_name =  FFpopSim_L<SeqLen>_N<MaxPopSize>_Ns<#Samplings>_Ts<SampleFreq>_Nv<SampleVol>_Mu<SubstitutionRate>_Amp<Amplitude>_Tfluct<FluctPer>_fluct, where the parameters are:

     * SeqLen: Sequence length
     * MaxPopSize: Maximal population size
     * #Samplings: Number of samples taken during the evolution
     * SampleFreq: Sampling frequency, in generations
     * SampleVol: Number of individual sampled (chosen randomly from population)
     * SubstitutionRate: substitution rate
     * Amplitude: PopSize fluctuation amplitude. MaxPopSize*Amplitude is the minimal population size (bottleneck)
     * FluctPer: Period of pop size fluctuations, in fraction of coalescent times (MaxPopSize)

TreeTime reconstruction results
===============================

 For each set of the simulation parameters, there is a TreeTime reconstruction file, which comprises the time-resolved population sizes (both real and reconstructed). Each file contains the header, which shows period and amplitude of the simulations. The table columns are following:
   * x: time, in generations
   * y: reconstructed pop size
   * trueY: true pop size

Plotting the results:
=====================

 Download the validation project, and run plot_skyline.py script. If needed,
 set the period and amplitude parameters in the script.







