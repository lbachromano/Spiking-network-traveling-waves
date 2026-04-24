# =============================================================================
# Spiking neural network simulation using Brian2
#
# Models a 2D cortical sheet with excitatory (pyramidal) and inhibitory
# (interneuron) populations.
#
# This code has external input sigmoid + noise.
# For a more general code with constant external input, see the file
# Simulations_const_input.py
#
# Usage: python myscript.py <index_print> <seed>
#   index_print : simulation index (used for output file naming)
#   seed        : random seed for reproducibility
#
# Output: .mat files saved in the 'results/' directory, including
#         LFP estimates and spike times/IDs for E and I populations.
# =============================================================================

from brian2 import *
prefs.codegen.target = 'numpy'  # use the Python fallback (can also use cython)
import numpy as np
import scipy.io as sio
import os
import sys
import random

if len(sys.argv) < 5:
    print("-" * 30)
    print("ERROR: Missing argument.")
    print("Usage: python Simulations_sigmoidal_external_rate.py <index_print> <seed> <nu_ext> <g>")
    print("Example: python Simulations_sigmoidal_external_rate.py 1 2")
    print("-" * 30)
    sys.exit(1)
index_print = int(sys.argv[1])
print("dx_e =", index_print) # simulation number
seed_value = int(sys.argv[2])
nu_ext_input=float(sys.argv[3])
g_input=float(sys.argv[4])
random.seed(seed_value)
np.random.seed(seed_value) # Ensure NumPy is seeded too
seed(seed_value)

# Output directory
results_dir = "results_f1_f2/examples"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
index_list=1 #network index. for plots 1 and 2 fixed to 1



######################################
# Parameters: Grid and neuron number #
######################################
Size_x=2050 #size of the array
Size_y=Size_x
Volume=Size_x*Size_y*umeter**2
grid_dist_r=400*umeter
grid_nr=int((Size_x-400 + grid_dist_r/umeter)//(grid_dist_r/umeter))
rows_r, cols_r = grid_nr,grid_nr
Total_nr_neurons=28000*(Volume/mm**2)
N_E=int(Total_nr_neurons*0.8)   # Pyramidal cells
N_I=int(Total_nr_neurons*0.2)   # Interneurons
N_R=grid_nr**2                  # Recording sites
# Setting up recording sites
eqs_r='''x : meter
         y : meter
         index : integer (constant)
'''
Rec = NeuronGroup(rows_r * cols_r, eqs_r)
Rec.x = '(i // rows_r) * grid_dist_r - (rows_r-1)/2.0  * grid_dist_r'
Rec.y = '(i % cols_r) * grid_dist_r - (cols_r-1)/2.0  * grid_dist_r'
# print the recording sites
sio.savemat(os.path.join(results_dir, 'Rec_x.mat'), {'data': Rec.x[:]})
sio.savemat(os.path.join(results_dir, 'Rec_y.mat'), {'data': Rec.y[:]})
##################
# External input #
##################
tau_x=50*ms
duration_ext=1800*ms
sigma_n = 3.3 * Hz/sqrt(ms)  # NOISE LEVEL VERY IMPORTANT
tau_n=150*ms
nu_baseline = nu_ext_input * Hz      # baseline rate
t0 = 1200 * ms               # sigmoid midpoint time
A_sigm = 1200 * Hz    # sigmoid step height from baseline to ceiling
k = 0.02 / ms                # sigmoid slope (steepness)
noise_sigma=0               # NO private noise!!!
##############
# Duration   #
##############
duration=2000*ms
duration_ext=duration
##################
# Single neuron  #
##################
V_th=18*mV
v_reset=11*mV
tau_mE=(20)*ms
tau_mI=(10)*ms
tau_ref_E=2*ms
tau_ref_I=1*ms
tau_rAI=0.7*ms
tau_rAE=0.7*ms
tau_dAE =3.5*ms
tau_dAI =3.5*ms
tau_rG=1*ms
tau_dG=18*ms
#####################
#  Synapses         #
#####################
vel_lat=0.3*mm/ms
E_rho = 0.6      #anisotropy of connectivity: E to E connections
I_rho = 0.0      #anisotropy of connectivity: all other connections
lam_dist=[380,250,250,250]  #spread of connectivity (u m)
EEprob_z=0.15   #parameter influencing the number of connections per neuron
J_ext_I=0.85 #strength of external connectoins to I
J_ext_E=0.455 #strength of external connectoins to E
g=g_input
J_EE=0.1575
J_EI=0.2625
J_IE=1.07*g*J_EE
J_II=g*J_EI

##################################
#  Defining the dynamics         #
##################################

start_scope()

# First, let's pre-generate the external rate time series
dt_val = defaultclock.dt
steps = int(duration_ext/dt_val)
time_array = np.arange(steps) * dt_val
x_global  = np.zeros(steps) * Hz
nu_global = np.zeros(steps) * Hz
en_global = np.zeros(steps) * Hz
en_global[0] = np.random.randn() * sigma_n * np.sqrt(tau_n)
for idx_i in range(steps):
    t = time_array[idx_i]
    sig_val  = A_sigm  / (1 + np.exp(-k * (t - t0)))
    rate_q   = nu_baseline + sig_val       # Quantity (Hz)
    if idx_i > 0:
        en_global[idx_i] = en_global[idx_i-1] + (-en_global[idx_i-1]/tau_n) * dt_val + sigma_n * np.sqrt(2 * dt_val) * np.random.randn()
    x_global [idx_i] = rate_q
    nu_global[idx_i] = max(0.0, rate_q + en_global[idx_i])
print('Global noise std:', np.std(en_global))
print('Mean absolute diff between noiseless and noisy:', np.mean(np.abs(nu_global - x_global)))


# Save external input rate to file
sio.savemat(os.path.join(results_dir, 'external_input_rate.mat'), {
    'time_ms':     (time_array / ms),
    'nu_global':   (nu_global / Hz),   # noisy rate
})
print(f"External input rate saved to '{results_dir}/external_input_rate.mat'")

# Create neuron groups: E
eqs_e='''x : meter
        y : meter
        dist_syn_e : meter
        dist_syn_i : meter
        rho : 1
        dv_E/dt = (-v_E + I_AE - I_GE)/tau_mE : volt (unless refractory)
        dI_AE/dt = (-I_AE + x_AE)/tau_dAE : volt
        dI_GE/dt = (-I_GE + x_GE)/tau_dG : volt
        dx_AE/dt = (-x_AE)/tau_rAE : volt
        dx_GE/dt = (-x_GE)/tau_rG : volt   
        # Create a noise variable with a differential equation  
        dprivate_noise/dt = -private_noise/ms + sqrt(2/ms)*xi_noise : 1
        input_rate = rate_array(t) + noise_sigma*private_noise*Hz : Hz  ##!!! noise_sigma is set to zero        
        noise_sigma : 1 (constant)
        index : integer (constant)
'''
E = NeuronGroup(N_E, eqs_e, threshold='v_E>V_th', reset='v_E=v_reset',
                refractory=tau_ref_E, method='euler')
E.I_AE=0 * volt
E.I_GE=0 * volt
E.x_AE = 0 * volt
E.x_GE = 0 * volt
E.rho=E_rho
E.x = np.random.uniform(low=-Size_x/2, high=Size_x/2, size=(1,N_E))*umeter
E.y = np.random.uniform(low=-Size_y/2, high=Size_y/2, size=(1,N_E))*umeter
E.dist_syn_e = lam_dist[0] * umeter
E.dist_syn_i = lam_dist[1] * umeter
E.noise_sigma = noise_sigma  ##!!! noise_sigma is set to zero, but this allows to introduce extra variability
E.private_noise = np.random.randn(N_E)

# Create neuron groups: I
eqs_i='''x : meter
         y : meter
         cc : 1
         dist_syn_e : meter
         dist_syn_i : meter
         rho : 1
         dv_I/dt = (-v_I + I_AI - I_GI)/tau_mI : volt (unless refractory)
         dI_AI/dt = (-I_AI + x_AI)/tau_dAI : volt
         dI_GI/dt = (-I_GI + x_GI)/tau_dG : volt
         dx_AI/dt = (-x_AI)/tau_rAI : volt
         dx_GI/dt = (-x_GI)/tau_rG : volt
         # Create a noise variable with a differential equation  
         dprivate_noise/dt = -private_noise/ms + sqrt(2/ms)*xi_noise : 1
         input_rate = rate_array(t) + noise_sigma*private_noise*Hz : Hz    ##!!! noise_sigma is set to zero                 
         noise_sigma : 1 (constant)
         index : integer (constant)
'''
I = NeuronGroup(N_I, eqs_i, threshold='v_I>V_th', reset='v_I=v_reset',
                refractory=tau_ref_I, method='euler')
I.I_AI=0 * volt
I.I_GI=0 * volt
I.x_AI = 0 * volt
I.x_GI = 0 * volt
I.x = np.random.uniform(low=-Size_x/2, high=Size_x/2, size=(1,N_I))*umeter
I.y = np.random.uniform(low=-Size_y/2, high=Size_y/2, size=(1,N_I))*umeter
I.cc= 2
I.dist_syn_e = lam_dist[2] * umeter
I.dist_syn_i = lam_dist[3] * umeter
I.rho=I_rho
I.noise_sigma = noise_sigma  # Hz, adjust this value to control variability between neurons
I.private_noise = np.random.randn(N_I)  # Initialize with random values

# Create a TimedArray for the global rate using the pre-generated series
rate_array = TimedArray(nu_global, dt=defaultclock.dt)
# Create Poisson spike generators for each neuron
PE = PoissonGroup(N_E, rates='input_rate')
PI = PoissonGroup(N_I, rates='input_rate')
# Link the PoissonGroup rates to the neuron group variables
PE.variables.add_reference('input_rate', E, 'input_rate')
PI.variables.add_reference('input_rate', I, 'input_rate')
# number of recorded neurons
N_REC_E=320
N_REC_I=80

##################################
#  Defining the synapses         #
##################################

S_PE = Synapses(PE, E, 'w : volt', on_pre='x_AE_post += w')
S_PE.connect(j='i')  # One-to-one connections
S_PE.w = 'J_ext_E*mV *tau_mE/tau_rAE'
S_PE.delay = 0*ms

S_PI = Synapses(PI, I, 'w : volt', on_pre='x_AI_post += w')
S_PI.connect(j='i')  # One-to-one connections
S_PI.w = 'J_ext_I*mV *tau_mI/tau_rAI'
S_PI.delay = 0*ms


EE = Synapses(E, E, 'w : volt', on_pre='x_AE_post += w')
EE.connect('i != j', p='0.91*0.084*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_pre))**2 + ((y_pre-y_post)/dist_syn_e_pre)**2 -2* rho_pre*((x_pre-x_post)/dist_syn_e_pre)*((y_pre-y_post)/dist_syn_e_pre))/sqrt(1-rho_pre*rho_pre))')
EE.w ='J_EE*mV *tau_mE/tau_rAE'
EE.delay = '''0.1*ms + sqrt((x_pre - x_post)**2 + (y_pre - y_post)**2) / vel_lat'''

# E to I. rho post, i.e. 0
# but dist i pre, i.e. 250
EI = Synapses(E, I, 'w : volt', on_pre='x_AI_post += w')
EI.connect(p='0.1206*exp(-sqrt(((x_pre-x_post)/(dist_syn_i_pre))**2 + ((y_pre-y_post)/dist_syn_i_pre)**2 -2* rho_post*((x_pre-x_post)/dist_syn_i_pre)*((y_pre-y_post)/dist_syn_i_pre))/sqrt(1-rho_post*rho_post))')
EI.w ='J_EI*mV  *tau_mI/tau_rAI'
EI.delay = '''0.1*ms + sqrt((x_pre - x_post)**2 + (y_pre - y_post)**2) / vel_lat'''

# I to E
IE = Synapses(I, E, 'w : volt', on_pre='x_GE_post += w')
IE.connect(p='0.1199*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_pre))**2 + ((y_pre-y_post)/dist_syn_e_pre)**2 -2* rho_pre*((x_pre-x_post)/dist_syn_e_pre)*((y_pre-y_post)/dist_syn_e_pre))/sqrt(1-rho_pre*rho_pre))')
IE.w ='J_IE*mV *tau_mE/tau_rG'
IE.delay = '''0.1*ms + sqrt((x_pre - x_post)**2 + (y_pre - y_post)**2) / vel_lat'''

# I to I
II = Synapses(I, I, 'w : volt', on_pre='x_GI_post += w')
II.connect('i != j', p='0.1199*exp(-sqrt(((x_pre-x_post)/(dist_syn_i_pre))**2 + ((y_pre-y_post)/dist_syn_i_pre)**2 -2* rho_pre*((x_pre-x_post)/dist_syn_i_pre)*((y_pre-y_post)/dist_syn_i_pre))/sqrt(1-rho_pre*rho_pre))')
II.w ='J_II*mV *tau_mI/tau_rG'
II.delay = '''0.1*ms + sqrt((x_pre - x_post)**2 + (y_pre - y_post)**2) / vel_lat'''

print(np.mean(EE.N_incoming[:]))
print(np.mean(EI.N_incoming[:]))
print(np.mean(IE.N_incoming[:]))
print(np.mean(II.N_incoming[:]))

avg_delay_EE = np.mean(EE.delay[:]) / ms
avg_delay_EI = np.mean(EI.delay[:]) / ms
avg_delay_IE = np.mean(IE.delay[:]) / ms
avg_delay_II = np.mean(II.delay[:]) / ms

##################################
#  Run simulations         #
##################################

# Set up monitors
M_AMPA_P = StateMonitor(E, 'I_AE', record=True, dt=1*ms)
M_GABA_P = StateMonitor(E, 'I_GE', record=True, dt=1*ms)

 
# Run the simulation with a single run command for the entire network
run(duration)


# COMPUTE LFPs
# Attention! this works because I am recording every ms
delay_AMPA = int(round((6*ms) / M_AMPA_P.clock.dt))
N = M_AMPA_P.I_AE.shape[1]
T_common = N - delay_AMPA
g_linear = M_AMPA_P.I_AE[:, :T_common] - 1.65 * M_GABA_P.I_GE[:, delay_AMPA:]
g_absum = np.abs(M_AMPA_P.I_AE[:, :T_common]) + np.abs(M_GABA_P.I_GE[:, :T_common])
def spatial_mask(k,ex,ey,Rx,Ry,LFP_g):
    R_x=150*umeter
    R_r=R_x*1.85
    f_x=(R_x/umeter)**2/sqrt(R_x/umeter)
    temp_dist=np.sqrt((ex[:]-Rx[k])**2 + (ey[:]-Ry[k])**2)
    f_r=f_x/((temp_dist/umeter)**2)
    f_r[temp_dist<=R_x]=1/sqrt(temp_dist[temp_dist<R_x]/umeter)
    norm_f=np.sum(f_r)
    LFP_t=np.transpose(LFP_g)
    LFP_k=np.dot(LFP_t[:,temp_dist<R_r], f_r[temp_dist<R_r])/norm_f  # Size = t * 1
    return LFP_k
def compute_LFP_for_g(LFP_g_matrix):
    LFP_out = np.zeros((T_common, N_R))
    for k in range(N_R):
        LFP_out[:, k] = spatial_mask(k, E.x, E.y, Rec.x, Rec.y, LFP_g_matrix)
    return LFP_out

# --- COMPUTE LFPs FOR EACH VARIANT ---
LFP_linear = compute_LFP_for_g(g_linear)
LFP_absum  = compute_LFP_for_g(g_absum)
def save_lfp(mat_name, LFP_matrix):
    adict = {}
    adict['data'] = LFP_matrix[:, :]  # sampled
    full_path = os.path.join(results_dir, mat_name)
    sio.savemat(full_path, adict)

base = f"i{index_print}_J{index_list}"
save_lfp(f"LFP_linearDelay_{base}.mat", LFP_linear)
save_lfp(f"LFP_absSum_{base}.mat", LFP_absum)
 
