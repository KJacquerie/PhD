# Defines output directory
const directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig2_MNIST_prediction"
#const task_ID="noisy"

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
include("model_MNIST_GB2012_fair.jl")
const experiment_name = "Graupner2012"
const gval_adjusted = 1.#[1 4.3 4.3 5]

const data = readdlm(@sprintf("%s/julia/testset_MNIST.dat", directory_name)) 

new_circuit=0

#idx_presented = collect(1:1:500)
idx_presented = [collect(1:1:20);collect(51:1:70);collect(101:1:120);collect(151:1:170)]#; collect(201:1:240); collect(251:1:270); collect(271:1:290); collect(301:1:340); collect(351:1:390); collect(401:1:440); collect(451:1:490)]

gCAMPA_state = readdlm(@sprintf("%s/data/%s/g_state.dat",directory_name, experiment_name))
w_init_state = readdlm(@sprintf("%s/data/%s/w_state.dat",directory_name, experiment_name))

gCAMPA_GB2012 = readdlm(@sprintf("%s/data/Graupner2012/g_state.dat",directory_name))
w_init_GB2012 = readdlm(@sprintf("%s/data/Graupner2012/w_state.dat",directory_name))


STATE1 =  97  #size(w_init_state)[2]-22 # - 11 # 1
STATE2 =  100#size(w_init_state)[2]    # -    # 6


const SD = "no"

#Pattern characterization



const NB_digits = 10
const NB_pixels = size(data, 1)

const NB_samples_tot = convert(Int64, size(data,2))
const NB_samples_digit = convert(Int64, size(data,2)/NB_digits)#nb de patterns par classe 


# Network parameters
const ncellsI = 1
const nPre= NB_pixels
const nPost= NB_digits
const ncellsC = nPre+nPost
const ncells = ncellsI+ncellsC

## Simulation parameters
const dt = 0.01


const N_cycles = 1 #31
const Duration_cycle = 2000
const Tdt_cycle = convert(Int64, Duration_cycle/dt)

const N_states = 1#N_cycles*2
const Duration_state = copy(Duration_cycle)#convert(Int64, Duration_cycle/2)
const Tdt_state = copy(Tdt_cycle)#convert(Int64, Duration_state/dt)

const N_samples = 1
const Duration_sample = convert(Int64, Duration_state/N_samples)
const Tdt_sample = convert(Int64, Duration_sample / dt)


const T = N_cycles*Duration_cycle 
const Tdt = convert(Int64, T / dt)
const t = range(dt, T, length = Tdt)

const Duration_set = Duration_cycle 
const Tdt_set = convert(Int64, Duration_set/dt)

const Tdt_transient= convert(Int64, 500 / dt)

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
const IappI = 0. #3
const IappC_pre = 0.
const IappC_post = 6.5

const spike_duration = 3
const IstepI = -1.2-IappI
const IstepC = 50.

IstepI_cell = IstepI .* ones(ncellsI)
IstepC_cell = IstepC .* ones(ncellsC)
Istep_cell = [IstepI_cell; IstepC_cell]


    
freq_Pattern = 55 #35 #Normal(55, 10)
freq_OutPattern = 0.0001 #Normal(1, 0.01)
freq_post_ok= 0.0001 #### PREDICTION CHANGED HERE
freq_post_ko = 0.0001

freq_rest = Normal(1, 0.1) 
freq_noise = Normal(15,5)



#pattern_presented = convert(Matrix{Int64}, zeros(N_states, 1))
#idx_presented =  readdlm(@sprintf("%s/data/idx_presented.dat",directory_name))
digit_presented = convert(Matrix{Int64}, zeros(N_samples*N_states, 1))


neurons_freq  = zeros(length(idx_presented), ncells)
neurons_freq[:,1].=1

#TESTING SET

#idx_kat=2

if(new_circuit==1)
    for idx_dgt= 1:1:length(idx_presented)
      
        idxPattern= getindex.(findall(x->x>0.5, data[:,convert(Int64,idx_presented[idx_dgt])]),1) .+1
        nb_pixel_ON = size(idxPattern,1)
        idxOutPattern= getindex.(findall(x->x<0.5, data[:,convert(Int64,idx_presented[idx_dgt])]),1) .+1

        neurons_freq[idx_dgt,idxPattern] .= freq_Pattern
        neurons_freq[idx_dgt, idxOutPattern] .= freq_OutPattern

        neurons_freq[idx_dgt, (ncellsI+nPre+1):end] .= freq_post_ko
        neurons_freq[idx_dgt, ncellsI+nPre+1+floor(Int64, (idx_presented[idx_dgt]-1)/NB_samples_digit)] = freq_post_ok 

    end

    neurons_freq[:,1].=1  
    
    writedlm(@sprintf("%s/data/neurons_freq.dat",directory_name), neurons_freq, header=false)
    #writedlm(@sprintf("%s/data/digit_presented.dat",directory_name), digit_presented, header=false)
else
    neurons_freq = readdlm(@sprintf("%s/data/neurons_freq.dat",directory_name))
    #new_circuit=1digit_presented = readdlm(@sprintf("%s/data/digit_presented.dat",directory_name))
end

#idx_state=10
#convert(Int64,idx_presented[idx_dgt])


ImageView.imshow(reshape(neurons_freq[2,2:end-10], 22, 22)')
#ImageView.imshow(reshape(data[:,91], 22, 22)')


## synaptic plasticity
const expm = "Control"
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
const tauL  = 400 # 1e6


##   CONNECTIVITY

const gIGABAA = (2. * ones(ncellsC))./ ncellsI
const gIGABAB = (1.5 * ones(ncellsC))./ ncellsI


 
for idx_state= STATE1:1:STATE2
    println(idx_state)
    w_init = reshape(w_init_state[:,idx_state], nPre, nPost)
    gCAMPA = reshape(gCAMPA_state[:,idx_state], nPre, nPost)
    
    wg_max = 1#maximum(w_init.*gCAMPA) #1
    wg_min = 0#minimum(w_init.*gCAMPA) #0

    wg_scale = 1#maximum(w_init_GB2012[:,idx_state].*gCAMPA_GB2012[:,idx_state]) #1

    Spkt_mat = zeros(NB_digits, length(idx_presented) )
    for idx_digit_tested = 1:1:length(idx_presented)
        println(convert(Int64,floor((idx_presented[idx_digit_tested]-1)/NB_samples_digit)))
        @time (Spkt) = simulateTOY_ncellsScenarioNMOD(
            ncells,
            ncellsI,
            ncellsC,
            Istep_cell,
            gCAMPA,
            w_init, 
            neurons_freq[idx_presented[idx_digit_tested],:], # TRAIN #neurons_freq[idx_digit_tested,:], # TEST
            wg_min, 
            wg_max, 
            wg_scale
        )
        Spkt_mat[:,idx_digit_tested]=Spkt
    end
    writedlm(@sprintf("%s/data/%s/spkt_mat_%d.dat",directory_name, experiment_name, idx_state), Spkt_mat, header=false)
      
end

#println(convert(Int64,floor((idx_presented[2]-1)/NB_samples_digit)))
