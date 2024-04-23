# Defines output directory
directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig5"

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
using DSP

# STRUCTURAL plasticity
#const tauG = 400 # (strong) or 800 (weak)
#gCAMPA = 0.001*ones(nPre,nPost)


## Include model

include("model_SNR_GB2012_VAR.jl")
include("PARAMS_cycle.jl")
const new_network=0
const fMax = 50
const bound_type = "SB"
const experiment_type = "low_dep"
experiment_name = "Graupner2012_VAR"



# Network parameters
const ncellsI = 1
const nPre = 100
const nPost = 1
const ncellsC = nPre+nPost
const ncells = ncellsI+ncellsC


## Simulation parameters
const dt = 0.01

const N_cycles=5
const Duration_cycle = 30000
const Tdt_cycle = convert(Int64, Duration_cycle/dt)

const N_patterns = 1
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

const state = zeros(Tdt,1)
for idx_cycle=1:1:N_cycles
    # wake
    state[(idx_cycle-1)*Tdt_set+1: (idx_cycle-1)*Tdt_set+Tdt_state] .= 0
    # sleep
    #if(SD=="yes")
    #    state[(idx_cycle-1)*Tdt_set+Tdt_state+1:(idx_cycle-1)*Tdt_set+Tdt_cycle ] .= -1
    #else
        state[(idx_cycle-1)*Tdt_set+Tdt_state+1:(idx_cycle-1)*Tdt_set+Tdt_cycle ] .= 1
    #end
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



if(new_network==1 && fMax==50)
    const gamma = 0.10 #1  #10% of variability in the network put 0.1 to get 10% etc
    const gl_cells = rand(Uniform(gl*(1-gamma),gl*(1+gamma)),ncells)
    const gNa_cells = rand(Uniform(gNa*(1-gamma),gNa*(1+gamma)),ncells)
    const gKd_cells = rand(Uniform(gKd*(1-gamma),gKd*(1+gamma)),ncells)
    const k1_cells = rand(Uniform(k1*(1-gamma),k1*(1+gamma)),ncells) #k1*ones(ncells) 
    const k2_cells = rand(Uniform(k2*(1-gamma),k2*(1+gamma)),ncells) #k2*ones(ncells) 
    const gH_cells = rand(Uniform(gH*(1-gamma),gH+(gamma*gH)),ncells)
    const gKCa_cells = rand(Uniform(gKCa*(1-gamma),gKCa*(1+gamma)),ncells)
    const gCaT_cells = rand(Uniform(gCaT*(1-gamma),gCaT*(1+gamma)),ncells)
    const g_cond = [gNa_cells'; gKd_cells'; gCaT_cells'; gH_cells'; gKCa_cells'; gl_cells'; k1_cells'; k2_cells']'
    writedlm(@sprintf("%s/data/%s/gion.dat",directory_name, experiment_name), g_cond, header=false)
else 
    gion = readdlm(@sprintf("%s/data/%s/gion.dat",directory_name, experiment_name))
    const gl_cells = gion[:,6]
    const gNa_cells = gion[:,1]
    const gKd_cells = gion[:,2]
    const k1_cells = gion[:,end-1]
    const k2_cells = gion[:,end]
    const gH_cells = gion[:,4]
    const gKCa_cells = gion[:,5]
    const gCaT_cells = gion[:,3]
end



const IappI = 3.
const IappC = 0.

const spike_duration = 3 
const IstepI = -1.2-IappI
const IstepC = 50.0 

IstepI_cell = IstepI .* ones(ncellsI)
IstepC_cell = IstepC .* ones(ncellsC)
Istep_cell = [IstepI_cell; IstepC_cell]

const N_states = N_cycles*2
if(new_network==1)
    gamma=0.1
    neurons_freq  = zeros(N_samples*N_states, ncells)   
    neurons_freq[:,1] .= 1 # inhibitory cell

    neurons_freq[:,2:6] = round.(rand(Uniform(73,76),N_samples*N_states,5))
    neurons_freq[:,7:end-1] = round.(rand(Uniform(0.1,5),N_samples*N_states,95))
    neurons_freq[:,end] .= 25.00
    writedlm(@sprintf("%s/data/%s/neurons_freq.dat", directory_name, experiment_name), neurons_freq, header=false)
else
    neurons_freq= readdlm(@sprintf("%s/data/%s/neurons_freq.dat",directory_name, experiment_name))
end
Iapp_cell = zeros(ncells, Tdt)

for idx_cycle= 1:1:N_cycles
    for idx=1:1:N_samples
        T1 = (idx_cycle-1)*Tdt_set +(idx-1)*Tdt_sample+1
        T2 = (idx_cycle-1)*Tdt_set + idx*Tdt_sample
        Iapp_cell[:, T1:T2] = get_Iapp(Duration_sample, dt, neurons_freq[(idx_cycle-1)*N_samples+idx,:], spike_duration)
    end
end
Iapp_cell[1:ncellsI,:] .= IappI


const  BurstTime = zeros(1,N_cycles)
for idx=1:1:N_cycles
     BurstTime[idx] = Duration_state + (idx-1)*Duration_set
end
const BurstDuration = Duration_state #20000

const StateTime = zeros(1,N_cycles*2)
let idx_count
    idx_count=1
    for idx_cycle=1:1:N_cycles
        StateTime[idx_count]   = (idx_cycle-1)*Duration_set+Duration_state
        StateTime[idx_count+1] = (idx_cycle-1)*Duration_set+Duration_cycle
        idx_count = idx_count+2
    end
end

## synaptic plasticity
const expm = "Control"

# SJO param
const tau_Ca = 22.6936 #[ms]
const C_Pre  = 0.56#17539
const C_Post = 1.24#23964
const D_pre  = 4.60#98 #[ms]

const tau_w    = 346.3615e3 #[s>ms]
const gamma_p  = 725.085*1.1#*0.7
const gamma_p_sleep  = 725.085*0.95
const gamma_d  = 331.909#*0.7
const theta_p = 1.3
const theta_d = 1.
const wfix = 0.5


##   CONNECTIVITY

const gIGABAA_unit = 2.0
const gIGABAB_unit = 1.5
if(new_network==1 && fMax==50)
    const gIGABAA = rand(Uniform(gIGABAA_unit*(1-gamma),gIGABAA_unit*(1+gamma)),ncellsC)./ ncellsI
    const gIGABAB = rand(Uniform(gIGABAB_unit*(1-gamma),gIGABAB_unit*(1+gamma)),ncellsC)./ ncellsI
    const g_syn = [gIGABAA'; gIGABAB']'
    writedlm(@sprintf("%s/data/%s/gsyn.dat",directory_name, experiment_name), g_syn, header=false)
else
    gsyn = readdlm(@sprintf("%s/data/%s/gsyn.dat",directory_name, experiment_name))
    const gIGABAA = gsyn[:,1]
    const gIGABAB = gsyn[:,2]
end


w_init = 0.5*ones(nPre,nPost)



#idx_wl = convert(Matrix{Int64}, readdlm("idx_wl.dat"))

gCAMPA_init_mat = [0.0001 0.0005 0.001 0.002 0.005 0.008 0.01]
tauG_mat = range(100, 1000, length = 10)
for idx_tauG=1:1:length(tauG_mat)
    #for idx_gCAMPA_init=1:1:length(gCAMPA_init_mat)
    idx_gCAMPA_init=2
        println("tauG=", tauG_mat[idx_tauG])
        println("gCAMPA0=", gCAMPA_init_mat[idx_gCAMPA_init])
        @time () = simulateTOY_ncellsScenarioNMOD(
            ncells,
            ncellsI,
            ncellsC,
            Iapp_cell,
            Istep_cell,
            gCAMPA_init_mat[idx_gCAMPA_init], 
            tauG_mat[idx_tauG], 
            idx_gCAMPA_init
        )
    #end
end

