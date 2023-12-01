# Loads packages

using Plots
using DelimitedFiles
using Statistics
using Images, ImageView
using DifferentialEquations
using DataFrames
using Printf
using CSV
using LinearAlgebra
using Distributions
using DataStructures

## Include modelT

include("model_SNR_ju.jl")

# Network parameters
const ncellsI = 1
const nPre = 3  #10 #Presynaptic cells
const nPost = 2 #1 #Postsynaptic cells
const ncellsC = nPre+nPost
const ncells = ncellsI+ncellsC

## Simulation parameters 
const dt = 0.01

const N_states= 4 #10
const Duration_state = 15000 #7500
const Tdt_state = convert(Int64, Duration_state/dt)

const T = N_states*Duration_state
const Tdt = convert(Int64, T / dt) 
const t = range(dt, T, length = Tdt)

const state = zeros(Tdt,1)


#=Depending on the scenario we want to test, we define the nature of our states:
tonic (without structural plasticity) --> 0 
burst (with structural plasticity) --> 1
inactive (without strucutral plasticity) --> -1
tonic_l (with strucutral plasticity) --> 2
inactive_l (with structural plasticity) --> -2
=#



for idx_state=1:2:N_states
    state[(Tdt_state*(idx_state-1)+1):Tdt_state*idx_state].=0
    state[((Tdt_state*idx_state)+1):(Tdt_state*(idx_state+1))].=1
end 




## Neurons' model parameters

# Global parameters
const C = 1
const VNa = 50
const VK = -85
const VCa = 120
const Vl = -55
const VH = -20
const Kd = 170

# Cells parameters
const gl = 0.055
const gNa = 170.0
const gKd = 40
const k1 = 1.e-1
const k2 = 0.1e-1
const gH = 0.01
const gKCa = 4
const gCaT = 0.55

const gamma = 0. #1  #10% of variability in the network put 0.1 to get 10% etc 
const gl_cells = rand((gl-(gamma*gl):0.001:gl+(gamma*gl)),ncells)
const gNa_cells = rand((gNa-(gamma*gNa):0.001:gNa+(gamma*gNa)),ncells)
const gKd_cells = rand((gKd-(gamma*gKd):0.001:gKd+(gamma*gKd)),ncells)
const k1_cells = rand((k1-(gamma*k1):0.001:k1+(gamma*k1)),ncells)
const k2_cells = rand((k2-(gamma*k2):0.001:k2+(gamma*k2)),ncells)
const gH_cells = rand((gH-(gamma*gH):0.001:gH+(gamma*gH)),ncells)
const gKCa_cells = rand((gKCa-(gamma*gKCa):0.001:gKCa+(gamma*gKCa)),ncells)
const gCaT_cells = rand((gCaT-(gamma*gCaT):0.001:gCaT+(gamma*gCaT)),ncells)



## Current parameters
const IappI = 3.
const IappC = 0.

const spike_duration = 3
const IstepI = -1.2-IappI
const IstepC = 50.

IstepI_cell = IstepI .* ones(ncellsI)
IstepC_cell = IstepC .* ones(ncellsC)
Istep_cell = [IstepI_cell; IstepC_cell]

neurons_freq = zeros(N_states, ncells)
neurons_freq[:,1] .= 1 # inhibitory cell

for idx_state=1:1:N_states
    if (state[((idx_state-1)*Tdt_state+1)]==-1||state[((idx_state-1)*Tdt_state+1)]==-2)
        neurons_freq[idx_state,2:end]= rand([0.1 0.5 1], size(neurons_freq[idx_state,2:end]))

    else

        #For Gonzalez's inspired task 
        #=possible_freq =  [ 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2]
        neurons_freq[:,2:(end-1)] = rand(possible_freq, size(neurons_freq[:,2:(end-1)]))
        neurons_freq[:,2] .= 60 # associative input (1)

        #neurons_freq[:,2:6] .= 60 # associative input
     

        # neurons_freq[:,13] .= 60 # associative input (12)
        #neurons_freq[:,20] .= 60 # associative input (19)
        #neurons_freq[:,31] .= 60 # associative input (30)
        # neurons_freq[:,70] .= 60 # associative input (69)
        #neurons_freq[:,82] .= 60 # associative input (81)
        
        neurons_freq[:,end].= 40 # postsynaptic cell=#

        if(idx_state==1||idx_state==2)
            neurons_freq[idx_state,2]=70
            neurons_freq[idx_state,3]=5 
            neurons_freq[idx_state,4]= 20 
            neurons_freq[idx_state,5]= 30 
            neurons_freq[idx_state,6]=40 
        end

        if(idx_state==3||idx_state==4)
            neurons_freq[idx_state,2]=70 
            neurons_freq[idx_state,3]= 5 
            neurons_freq[idx_state,4]= 20
            neurons_freq[idx_state,5]= 30 
            neurons_freq[idx_state,6]= 40
        end

        #=if(idx_state==5||idx_state==6)
            neurons_freq[idx_state,2]=70
            neurons_freq[idx_state,3]= 5 #10
            neurons_freq[idx_state,4]= 20 #0.2
            neurons_freq[idx_state,5]= 30 #30
            neurons_freq[idx_state,6]= 40 #40
        end=#

    end
end

Iapp_cell = zeros(ncells, Tdt)

for idx_state=1:1:N_states
    if ((state[((idx_state-1)*Tdt_state+1)]==0)||(state[((idx_state-1)*Tdt_state+1)]==-1)||(state[((idx_state-1)*Tdt_state+1)]==2||(state[((idx_state-1)*Tdt_state+1)]==-2)))
        Iapp_cell[:,(((idx_state-1)*Tdt_state+1):(idx_state*Tdt_state))]=get_Iapp(Duration_state, dt, neurons_freq[idx_state,:], spike_duration)
    end
end
Iapp_cell[1:ncellsI,:] .= IappI

const StateTime = zeros(1,N_states)
for idx=1:1:N_states
    StateTime[idx]=(idx-1)*Duration_state
end

 
## synaptic plasticity
const expm = "Control"
const tau_Ca = 22.27212 #[ms]
const C_Pre  = 0.84410
const C_Post = 1.62138
const D_pre  = 9.53709 #[ms]

const tau_w    = 520.76129e3 #[s>ms]
const gamma_p  = 597.08922
const gamma_d  = 137.7586
const n_lin = 1.
const theta_p = 2.009289
const theta_d = 1.0

const Omega_p = gamma_p/(gamma_p+gamma_d)
const Omega_sleep = 0.75
const Omega_d = 0.
const Omega_0 = 0.

const tauw_p = tau_w /(gamma_p + gamma_d)
const tauw_d = tau_w / gamma_d
const tauw_0 = 0.
const zeta = tauw_d / tauw_p
const tau_x = 1  

const tau_l = 10. #3e1 #5e1  #5e5
const tauINCR = 100
l_init = 0.1*ones(nPre,nPost) # initial connectivity (initial late-weight)
w_init = 0.5*ones(nPre,nPost) # initial early-weight


##   CONNECTIVITY

# Connectivity of the inhibitory neuron
const gIGABAA = (2. * ones(ncellsC))./ ncellsI 
const gIGABAB = (1.5 * ones(ncellsC))./ ncellsI


# Pre to post

idxPrePost = [CartesianIndex(1,ncells-2)] # idxPrePost shows us the connections between neurons
for idx_q = 2:1:nPre
  push!(idxPrePost, CartesianIndex(idx_q, ncells-2)) 
end

for idx_q = 1:1:nPre
    push!(idxPrePost, CartesianIndex(idx_q, ncells-1)) 
end




@time (w, l_plot, l_incr, w_dot_plot, l_conv, VV) = simulateTOY_ncellsScenarioNMOD(
    ncells,
    ncellsI,
    ncellsC,
    Iapp_cell,
    Istep_cell,
    l_init,
    idxPrePost
)



#writedlm(@sprintf("/Users/magis/Documents/Justine/TFE/DataFig/FigNewRule/Tonic-Burst_l/w.dat"), w, header=false)
#writedlm(@sprintf("/Users/magis/Documents/Justine/TFE/DataFig/FigNewRule/Tonic-Burst_/l.dat"), l_plot, header=false)


p1 =plot(w[1:end,:]', legend= (frameon=false))
p2 = plot(l_plot[1:end,:]',legend= (frameon=false))
p3 = plot(w[1:end,:]'.*l_plot[1:end,:]',legend= (frameon=false))   
#p4 = plot(VV[1:end,:]', legend= (frameon=false))

plot(p1,p2,p3, layout=(3,1))
#plot(p1,p2,p3,p4, layout=(4,1))
