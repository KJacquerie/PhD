# Defines output directory
const directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig6"
const task_ID="noisy"

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

## Include model

#include("model_patterns_AMPA.jl")
include("model_GB2012.jl")
model = "Graupner2012" #or OLD

const SD = "no"
const gINIT = 0.02 #range(0.55, 0.75,5): va de 0.55 à 0.75 en 5 steps #range(0.05, 0.5,10)
const new_circuit=1

#Pattern characterization
const data = readdlm("dataset_2BarrePattern_OL.dat")
const NB_pixels = size(data, 1)
const NB_patterns = convert(Int64, size(data,2)) 
const NB_classe = 2
const NB_patterns_classe = convert(Int64, size(data,2)/NB_classe) #nb de patterns par classe 

experiment_name = "OL"



# Network parameters
const ncellsI = 1
const nPre= NB_pixels
const nPost= NB_classe
const ncellsC = nPre+nPost
const ncells = ncellsI+ncellsC

## Simulation parameters
const dt = 0.01

const N_cycles = 8
const Duration_cycle = 30000
const Tdt_cycle = convert(Int64, Duration_cycle/dt)

const N_states = N_cycles*2
const Duration_state = convert(Int64, Duration_cycle/2)
const Tdt_state = convert(Int64, Duration_state/dt)

const N_samples = 1
const Duration_sample = convert(Int64, Duration_state/N_samples)
const Tdt_sample = convert(Int64, Duration_sample / dt)


const T = N_cycles*Duration_cycle 
const Tdt = convert(Int64, T / dt)
const t = range(dt, T, length = Tdt)

const Duration_set = Duration_cycle 
const Tdt_set = convert(Int64, Duration_set/dt)

## States 
#=
    # METHODE 1 : definition of the states in function of the time 
    const state = zeros(Tdt,1)
    for idx_cycle=1:1:N_cycles
    # wake
    state[(idx_cycle-1)*Tdt_set+1: (idx_cycle-1)*Tdt_set+Tdt_state] .= 0
    # sleep
    state[(idx_cycle-1)*Tdt_set+Tdt_state+1:(idx_cycle-1)*Tdt_set+Tdt_cycle ] .= 1
    
end
=#
 
    # METHODE 2: Pour chancun des states, on dit ce qu'il est :
    #=
    tonic-->0
    burst--> 1
    rest --> -1
    tonic_g --> 2
    rest_g --> -2
    =#

const state = zeros(Tdt,1)

for i=1:1:N_states
    if(mod(i,2)==1)
        state[((Tdt_state*(i-1))+1):(Tdt_state*i)].=0
    else
        state[((Tdt_state*(i-1))+1):(Tdt_state*i)].=1
    end
end 


#=
state[1:Tdt_state,1].=0
state[(Tdt_state+1):(Tdt_state*2),1].=-1
state[((Tdt_state*2)+1):(Tdt_state*3),1].=0
state[((Tdt_state*3)+1):(Tdt_state*4),1].=-1
=#

#Definition of the state time

const StateTime = zeros(1,N_cycles*2)
let idx_count
    idx_count=1
    for idx_cycle=1:1:N_cycles
        StateTime[idx_count]   = (idx_cycle-1)*Duration_set+Duration_state
        StateTime[idx_count+1] = (idx_cycle-1)*Duration_set+Duration_cycle
        idx_count = idx_count+2
    end
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

const gamma = 0.0 #1  #10% of variability in the network put 0.1 to get 10% etc
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


    
freq_Pattern = 40 #Normal(55, 10)
freq_OutPattern = 0.01 #Normal(1, 0.01)
freq_post_ok= 55
freq_post_ko = 0.0001

freq_rest = Normal(1, 0.1) #old 0.1; 0.01


neurons_freq  = zeros(N_states, ncells)
neurons_freq[:,1].=1
pattern_presented = convert(Matrix{Int64}, zeros(N_states, 1))


#LEARNING SET
    #activity of input & output neurons 

for idx_state= 1:1:N_states #(*N_samples)

    if (state[((idx_state-1)*Tdt_state+1)]==-1)
        neurons_freq[idx_state,2:end] = abs.(rand(freq_rest, size(neurons_freq[idx_state,2:end])))
    else 

      #Pattern presented 

    #On présente 1 barre de chaque classe (aléatoire entre les 2) à chaque cycle
    #=
    if (mod(idx_state,4) == 1||mod(idx_state,4)==2) #4 car 4 patterns --> à adapter si plus
        pattern_presented[idx_state] = rand(1:2,1)[1] 
    else 
        pattern_presented[idx_state] = rand(3:4,1)[1]
    end 
    =#
            pattern_presented[1] = 1[1]
            pattern_presented[2] = 1[1]
            pattern_presented[3] = 3[1]
            pattern_presented[4] = 3[1]
            pattern_presented[5] = 2[1]
            pattern_presented[6] = 2[1]
            pattern_presented[7] = 4[1]
            pattern_presented[8] = 4[1]


            pattern_presented[9] = 1[1]
            pattern_presented[10] = 1[1]
            pattern_presented[11] = 3[1]
            pattern_presented[12] = 3[1]
            pattern_presented[13] = 2[1]
            pattern_presented[14] = 2[1]
            pattern_presented[15] = 4[1]
            pattern_presented[16] = 4[1]

    # Répéter le pattern choisi durant le 2e state du même cycle 
    if(mod(idx_state,2)==0)
        neurons_freq[idx_state,2:end] .= neurons_freq[(idx_state - 1),2:end]
    else    

    #pattern_presented[idx]=rand(1:NB_patterns,1)[1]




        idxPattern= getindex.(findall(x->x>0.5, data[:,pattern_presented[idx_state]]),1) .+1

        idxOutPattern= getindex.(findall(x->x<0.5, data[:,pattern_presented[idx_state]]),1) .+1
        #neurons_freq[idx,idxPattern] .= abs.(rand(freq_Pattern, size(neurons_freq[idx,idxPattern])))
        
        neurons_freq[idx_state,idxPattern] .= freq_Pattern
        #neurons_freq[idx,idxOutPattern] = abs.(rand(freq_OutPattern, size(neurons_freq[idx,idxOutPattern])))
        neurons_freq[idx_state, idxOutPattern] .= freq_OutPattern

        neurons_freq[idx_state, (ncellsI+nPre+1):end] .= freq_post_ko
        neurons_freq[idx_state, ncellsI+nPre+1+floor(Int64, (pattern_presented[idx_state]-1)/NB_classe)] = freq_post_ok 
        
        #si plusieurs neurones post: 
        #neurons_freq[idx, (ncellsI+nPre+1):end] .= 0.01
        #neurons_freq[idx, ncellsI+nPre+1+floor(Int64, (pattern_presented[idx]-1))] = 30
      #end 
    end 
    end 
    # freq neurone inhibiteur
    neurons_freq[:,1].=1  
    #writedlm(@sprintf("%s/data/task_%s/neurons_freq.dat",directory_name, task_ID), neurons_freq, header=false)


end

   

# creation of the applied current (square waves) 

Iapp_cell = zeros(ncells, Tdt)

for idx_state= 1:1:N_states
    #for idx=1:1:N_samples
        T1 = (idx_state-1)*Tdt_state +1
        T2 = (idx_state)*Tdt_state
        Iapp_cell[:, T1:T2] = get_Iapp(Duration_sample, dt, neurons_freq[idx_state,:], spike_duration)
   # end
end
Iapp_cell[1:ncellsI,:] .= IappI


const  BurstTime = zeros(1,N_cycles)
for idx=1:1:N_cycles
     BurstTime[idx] = Duration_state + (idx-1)*Duration_set
end
const BurstDuration = Duration_state #20000

#StateTime???

## synaptic plasticity
const expm = "Control"
#=
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
const Omega_sleep = 0.6
const Omega_d = 0.3
const Omega_0 = 0.

const tauw_p = tau_w /(gamma_p + gamma_d)
const tauw_d = tau_w / gamma_d
const tauw_0 = 0.
const zeta = tauw_d / tauw_p
const tau_x = 1
=#
const tau_Ca = 22.6936 #[ms]
const C_Pre  = 0.56#17539
const C_Post = 1.24#23964
const D_pre  = 4.60#98 #[ms]

const tau_w    = 346.3615e3 #[s>ms]
const gamma_p  = 725.085*1.1
const gamma_p_burst  = 725.085*0.95
const gamma_d  = 331.909
#const n_lin = 1.
const theta_p = 1.3
const theta_d = 1.
const wfix = 0.5

## LATE-PHASE PLASTICITY
const tauL  = 200 # 1e6


##   CONNECTIVITY

const gIGABAA = (2. * ones(ncellsC))./ ncellsI
const gIGABAB = (1.5 * ones(ncellsC))./ ncellsI


# Pre to post
idxPrePost = [CartesianIndex(1,ncellsI+nPre)]
for idx_pattern = 1:1:NB_classe
    if(idx_pattern==1)
        idx_q_start=2
    else 
        idx_q_start=1
    end
    for idx_q=idx_q_start:1:nPre
        push!(idxPrePost,CartesianIndex(idx_q,ncellsI + nPre+idx_pattern-1))
    end   
end 


gCAMPA = copy(gINIT)*ones(nPre,nPost) #copy(gINIT_mat[idx_gINIT-10])/4*ones(nPre,nPost)  
w_init = 0.5*ones(nPre, nPost) #zeros(ncellsC,ncellsC)
    
#ici divisé par 4 car 4 pixels 


    @time (w, gCAMPA_plot, w_state, g_state) = simulateTOY_ncellsScenarioNMOD(
        ncells,
        ncellsI,
        ncellsC,
        Iapp_cell,
        Istep_cell,
        gCAMPA,
        w_init, 
        idxPrePost )


writedlm(@sprintf("%s/data/%s/%s/w_state%s.dat", directory_name, model, experiment_name, N_cycles), w_state, header = false)
writedlm(@sprintf("%s/data/%s/%s/g_state%s.dat", directory_name, model, experiment_name, N_cycles), g_state, header = false)
