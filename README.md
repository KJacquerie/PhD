
# Code for the article 
    'Switches from tonic to burst firing enable memory consolidation through late-phase synaptic plasticity    
by Kathleen Jacquerie, Danil Ty mankov, Pierre Sacre, Guillaume Drion.  

The scripts are mainly written by Kathleen Jacquerie. 
Caroline Minne, J iette Ponnet, and Nora Benghalem worked on this code for their master theses.


Each Figure from the article is associated with a corresponding folder. 
Each folder is divided into the same structure:  
/data - contains the parameter used to generate the sim ation 
/  - main scripts to run the sim ations.  
All computational experiments are launched thanks to a *Simu_... .jl* calling a function with the differential equations of the conductance-based models (neuronal model) and the synaptic weight change (plasticity r es) found in *model_... .jl*.  
/matlab - codes to analyze and plot the network sim ations.  

disclaimers: 
l_ij is indicated at g_AMPA_ij 

### Fig1
Simulate a network of neurons with / /Simu_scenario_GB2012.jl 
Inputs
gion.dat: intrinsic parameters of the neurons 
gsyn.dat: initial parameters of the late weights 
neurons_freq.dat: frequencies associated with each neuron at the different states.

Outputs
- Voltage recordings
- Firing pattern properties
- LFP recordings
- Synaptic weight recordings

MATLAB: Plot traces and quantify firing pattern properties with /matlab/Fig1_plotV.m 

### Fig2_MNIST_long
Simulate the pattern learning task in different configurations of switches and plasticity:
Figure 2B: Simu_MNIST_GB2012.jl
Figure 2C: Simu_MNIST_GB2012_learn.jl
Inputs
digit_presented.dat: contains the sequence of digits sequentially learned during the simulation. 
idx_presented.dat: contains the index of the sample learned at each state picked from the dataset. 
neurons_freq.dat: frequencies associated with each neuron at the different states.
neurons_freq_learn.dat: same but for the additional tonic firing states (figure 2C)

   Outputs     
w_state.dat: values of the early weights at the end of each state.  
g_state.dat: values of the late weights at the end of each state.  
 
MATLAB:  
RF_MNIST.m: code to obtain the weight matrix in the image 
RF_MNIST_corr.m: code to compute the correlation between the mean dataset (saved in mean_MNIST_20.dat for Figure 2D and mean_MNIST_learn.dat for Figure 2E) 
 
### Fig3 
Simulate the evolution of the early weight for several switches of firing activities.  
Simu_scenario_GB2012_RESET.jl associated with Figure 3AB 
Simu_TUNE.jl associated with Fig SI for one tonic firing state followed by three burst firing states with different hyperpolarizing currents.  
 
Inputs 
gion.dat: intrinsic parameters of the neurons 
gsyn.dat: initial parameters of the late weights 
neurons_freq.dat: frequencies associated with each neuron at the different states. 
 
  Outputs     
- Voltage recordings 
- Firing pattern properties 
- LFP recordings 
- Synaptic weight recordings 
 
MATLAB 
Fig3AB.m: code to obtain w(t), histograms, and scatter plot - different forms of normalization are encoded.  
Fig3C_tune.m: code associated with SIMU_TUNE.jl 
 
 
### Fig3_Demo 
Simu_Graupner2016.jl  code used to sim ate the evolution of the early weight for the model of Graupner et al. 2016.  
Simu_PairBased.jl code used to sim ate the evolution of the early weight in the Pair-Based model.  
This sim ation can be done with parameter sets obtained from the cortex (CTX) or the hippocampus (HPC) [region] depending on the paper. The weight dependency can be soft-bound (SB) or hard-bound (HB) [bound_type].  
 
Inputs 
gion.dat 
gsyn.dat 
neurons_freq.dat 
 
Outputs     
- Voltage recordings 
- Firing pattern properties 
- Synaptic weight recordings 
 
MATLAB 
Fig3_Demo_SB.m for Fig SI is associated with the analytical prediction of the converging point.  
set new_curves =1 to generate new traces. The code takes time since it loads all the data.  
This code traces the scatter plot. It also compares the analytical prediction and the sim ation values.  
 
###  Fig3_Model 

Simu_XX.jl runs the different sim ations for the model     XX     (can be Graupner, Deperrois,  '85).  
 
MATLAB 
Fig3_scatter_SB.m retrieves the scatter plot for soft-bound models 
Fig3_scatter_HB.m retrieves the scatter plot for hard-bound models 
 
 
### Fig4 

Simu_scenario_GB2012.jl: code for Fig4B 
Simu_scenario_GB2012_RESET.jl: code for Fig4C 
Simu_scenario_GB2012_noBURST.jl: code for Fig4D 
 
 Inputs     
gion.dat 
gsyn.dat 
neurons_freq.dat 
 
Outputs     
- Voltage recordings 
- Firing pattern properties 
- Synaptic weight recordings 
 
MATLAB 
Fig4_wl.m plot traces of the early weights, late weights, and total weights.  
 
### Fig4_MNIST 
go to folder Fig2_MNIST_long because it replicates these figures with different combinations of switches and plasticity.  
Figure 4E is associated with the code Simu_MNIST_GB2012_RESET.jl 
Figure 4F is accociated with the code Simu_MNIST_GB2012_noBURST.jl 
 
 
### Fig5 
Simu_SNR_GB2012_VAR.jl sim ates the network with different initial conditions for eta.  
Disclaimer: eta=1/tau or 1/tauG 
and l(0) is gCAMPA.  
 
Inputs     
gion.dat 
gsyn.dat 
neurons_freq.dat 
Outputs 
w_BACK: early-weight at the end of each state.  
g_BACK: late-weight at the end of each state. 
 
MATLAB 
Fig5_var_quantif.m plots the evolution of the SNR for the different parameter conditions and compute the category between the four microscopic trends.  
 
### Fig 6 
Comparison between 2012 model (folder: Graupner2012) and 2016 model (folder: Graupner2016) between non-overlapping patterns (folder: nOL) and overlapping patterns (folder: OL).  

Inputs     

Outputs 
w_state8: early-weight at the end of each state.  
g_state8: late-weight at the end of each state. 
 
 
MATLAB 
Fig6_RF.m: code to obtain Figure Eiii, Eiv, Fiii, Fiv for the weight matrices. Change the name of the experiment  [model] (Graupner2012, Graupner2016) or the [experiment_name]  
Fig6_quantif.m: code to obtain the correlation traces.  
 

}
