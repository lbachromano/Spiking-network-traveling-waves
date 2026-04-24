#

# =============================================================================
# Spiking neural network simulation using Brian2
#
# Models a 2D cortical sheet with excitatory (pyramidal) and inhibitory
# (interneuron) populations.
#
# This code has constant external input.
#
# Distance between recordings sites here is grid_dist_r = 141*umeter
#
# Output: .mat files saved in the 'results/' directory, including
#         LFP estimates and spike times/IDs for E and I populations.
# =============================================================================

from brian2 import *
prefs.codegen.target = 'numpy'  # use the Python fallback (can also use cython)
import numpy as np
import scipy.io as sio
import os
import random

start_scope()


if len(sys.argv) < 7:
    print("-" * 30)
    print("ERROR: Missing argument.")
    print("Usage: python Simulations_const_ext_rate.py <index_print> <seed> <n_trials> <homogeneous> <nu_ext> <rho_EE>")
    print("Example: python Simulations_const_ext_rate.py 1 42")
    print("-" * 30)
    sys.exit(1)
index_print = int(sys.argv[1])
seed_value = int(sys.argv[2])   # random seed
nr_iterations = int(sys.argv[3])   # defines the length of the simulation
homogeneous = int(sys.argv[4])   # flag for the connectivity structure
nu_ext = float(sys.argv[5])
EE_rho = float(sys.argv[6])
random.seed(seed_value)
np.random.seed(seed_value) # Ensure NumPy is seeded too
seed(seed_value)
print('External field value is', nu_ext)

#parameters:  anisotropy
I_rho = 0.0 #both I to E and I to I
EI_rho = 0.0
C0 = 0.1203 # setting the number of incoming connection
C0_EE= 0.0609 # setting the number of incoming connection

if homogeneous == 1:
    results_dir = 'results_f3'
else:
    results_dir = 'results_f4_f5'
    
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
index_list=1 #network index. for plots 1 and 2 fixed to 1


# Grid and neuron parameters (unchanged)
Size_x=2050
Size_y=Size_x
Volume=Size_x*Size_y*umeter**2

vel_lat=0.3*mm/ms
size_x = Size_x * umeter
size_y = Size_y * umeter

grid_dist_r = 141*umeter
grid_nr=int((Size_x-400 + grid_dist_r/umeter)//(grid_dist_r/umeter))
rows_r, cols_r = grid_nr,grid_nr

Total_nr_neurons=28000*(Volume/mm**2)
N_E=int(Total_nr_neurons*0.8)   # Pyramidal cells
N_I=int(Total_nr_neurons*0.2)   # Interneurons
N_R=grid_nr**2                  # Recording sites

# Setting up recording sites (unchanged)
eqs_r='''x : meter
         y : meter
         index : integer (constant)
'''
Rec = NeuronGroup(rows_r * cols_r, eqs_r)
Rec.x = '(i // rows_r) * grid_dist_r - (rows_r-1)/2.0  * grid_dist_r'
Rec.y = '(i % cols_r) * grid_dist_r - (cols_r-1)/2.0  * grid_dist_r'

sio.savemat(os.path.join(results_dir, 'Rec_x.mat'), {'data': Rec.x[:]})
sio.savemat(os.path.join(results_dir, 'Rec_y.mat'), {'data': Rec.y[:]})




# draw N_REC_E unique random neuron-indices from [0, N_E)


##################
# External input #
##################
 
nu_baseline =nu_ext*2195 * Hz # nu_{theta}=2195
warmup=200*ms
duration=1000*ms
duration_ext=warmup+nr_iterations*duration

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

#################
#  Synapses
#################

lam_dist=[380,250,250,250] #conectivity spatial spread (u m)
EEprob_z=0.15

J_ext_I=0.85
J_ext_E=0.455

g=4.1
J_EE=0.1575
J_EI=0.2625
J_IE=1.07*g*J_EE
J_II=g*J_EI


# Create the excitatory neurons with private external input
eqs_e='''x : meter
        y : meter
        dist_syn_e : meter
        dist_syn_i : meter
        rhoe : 1
        rhoi : 1
        C_EE : 1
        C_EtoI : 1
        dv_E/dt = (-v_E + I_AE - I_GE)/tau_mE : volt (unless refractory)
        dI_AE/dt = (-I_AE + x_AE)/tau_dAE : volt
        dI_GE/dt = (-I_GE + x_GE)/tau_dG : volt
        dx_AE/dt = (-x_AE)/tau_rAE : volt
        dx_GE/dt = (-x_GE)/tau_rG : volt
        index : integer (constant)

'''
E = NeuronGroup(N_E, eqs_e, threshold='v_E>V_th', reset='v_E=v_reset',
                refractory=tau_ref_E, method='euler')
E.I_AE=0 * volt
E.I_GE=0 * volt
E.x_AE = 0 * volt
E.x_GE = 0 * volt
E.x = np.random.uniform(low=-Size_x/2, high=Size_x/2, size=(1,N_E))*umeter
E.y = np.random.uniform(low=-Size_y/2, high=Size_y/2, size=(1,N_E))*umeter
E.dist_syn_e = lam_dist[0] * umeter
E.dist_syn_i = lam_dist[1] * umeter
E.rhoe=EE_rho
E.rhoi=EI_rho
E.C_EE = C0_EE / np.sqrt(1 - EE_rho**2) #correction due to anisotropy
E.C_EtoI=C0
 

# Create the inhibitory neurons with private external input
eqs_i='''x : meter
         y : meter
         cc : 1
         dist_syn_e : meter
         dist_syn_i : meter
         rho : 1
         C_II : 1
         C_ItoE : 1
         dv_I/dt = (-v_I + I_AI - I_GI)/tau_mI : volt (unless refractory)
         dI_AI/dt = (-I_AI + x_AI)/tau_dAI : volt
         dI_GI/dt = (-I_GI + x_GI)/tau_dG : volt
         dx_AI/dt = (-x_AI)/tau_rAI : volt
         dx_GI/dt = (-x_GI)/tau_rG : volt
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
I.C_ItoE=C0
I.C_II=C0
  


# Create Poisson spike generators for each neuron
PE = PoissonGroup(N_E, rates='nu_baseline')
PI = PoissonGroup(N_I, rates='nu_baseline')

# Connect the external inputs to the neurons
S_PE = Synapses(PE, E, 'w : volt', on_pre='x_AE_post += w')
S_PE.connect(j='i')  # One-to-one connections
S_PE.w = 'J_ext_E*mV *tau_mE/tau_rAE'
S_PE.delay = 0*ms

S_PI = Synapses(PI, I, 'w : volt', on_pre='x_AI_post += w')
S_PI.connect(j='i')  # One-to-one connections
S_PI.w = 'J_ext_I*mV *tau_mI/tau_rAI'
S_PI.delay = 0*ms

if homogeneous == 1:
    N_REC_E=800
    N_REC_I=200
    selected_E = np.random.choice(N_E, N_REC_E, replace=False).tolist()
    selected_I = np.random.choice(N_I, N_REC_I, replace=False).tolist()

    EE = Synapses(E, E, 'w : volt', on_pre='x_AE_post += w')
    EE.connect('i != j', p='0.0086')
    EE.w ='J_EE*mV *tau_mE/tau_rAE'
    EE.delay = '2.57*ms'
    # E to I
    EI = Synapses(E, I, 'w : volt', on_pre='x_AI_post += w')
    EI.connect(p='0.0086')
    EI.w ='J_EI*mV  *tau_mI/tau_rAI'
    EI.delay = '1.725*ms'
    # I to E
    IE = Synapses(I, E, 'w : volt', on_pre='x_GE_post += w')
    IE.connect(p='0.0086')
    IE.w ='J_IE*mV *tau_mE/tau_rG'
    IE.delay = '1.725*ms'
    # I to I
    II = Synapses(I, I, 'w : volt', on_pre='x_GI_post += w')
    II.connect('i != j', p='0.0086')
    II.w ='J_II*mV *tau_mI/tau_rG'
    II.delay = '1.725*ms'

else:
    
    N_REC_E=80000
    N_REC_I=20000
    selected_E = np.random.choice(N_E, N_REC_E, replace=False).tolist()
    selected_I = np.random.choice(N_I, N_REC_I, replace=False).tolist()
    
    EE = Synapses(E, E, 'w : volt', on_pre='x_AE_post += w')
    EE.connect('i != j',p='C_EE_pre*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_pre))**2 + ((y_pre-y_post)/(dist_syn_e_pre))**2 - 2*rhoe_pre*((x_pre-x_post)/dist_syn_e_pre)*((y_pre-y_post)/dist_syn_e_pre))/sqrt(1-rhoe_pre*rhoe_pre))')
    EE.w ='J_EE*mV *tau_mE/tau_rAE'
    EE.delay = '''0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2) / vel_lat'''

    # E to I
    EI = Synapses(E, I, 'w : volt', on_pre='x_AI_post += w')
    EI.connect(p='C_EtoI_pre*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_post))**2 + ((y_pre-y_post)/(dist_syn_e_post))**2-2*rhoi_pre*((x_pre-x_post)/dist_syn_e_post)*((y_pre-y_post)/dist_syn_e_post))/sqrt(1-rhoi_pre*rhoi_pre))')
    EI.w ='J_EI*mV  *tau_mI/tau_rAI'
    EI.delay = '''0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2) / vel_lat'''

    # I to E
    IE = Synapses(I, E, 'w : volt', on_pre='x_GE_post += w')
    IE.connect(p='C_ItoE_pre*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_pre))**2 + ((y_pre-y_post)/(dist_syn_e_pre))**2 - 2*rho_pre*((x_pre-x_post)/dist_syn_e_pre)*((y_pre-y_post)/dist_syn_e_pre))/sqrt(1-rho_pre*rho_pre))')
    IE.w ='J_IE*mV *tau_mE/tau_rG'
    IE.delay = '''0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2) / vel_lat'''

    # I to I
    II = Synapses(I, I, 'w : volt', on_pre='x_GI_post += w')
    II.connect('i != j', p='C_II_pre*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_pre))**2 + ((y_pre-y_post)/(dist_syn_e_pre))**2 - 2*rho_pre*((x_pre-x_post)/dist_syn_e_pre)*((y_pre-y_post)/dist_syn_e_pre))/sqrt(1-rho_pre*rho_pre))')
    II.w ='J_II*mV *tau_mI/tau_rG'
    II.delay = '''0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2) / vel_lat'''
    
    x_sel = (E.x[selected_E] / umeter).reshape(-1)   # shape (N_REC_E,)
    y_sel = (E.y[selected_E] / umeter).reshape(-1)
    full_path = os.path.join(results_dir,'selected_E_locations.mat')
    sio.savemat(full_path,
        {
            'selected_ids': selected_E,   # neuron indices in E
            'x_um':          x_sel,       # x‐positions [μm]
            'y_um':          y_sel        # y‐positions [μm]
        }
    )

    Ix_sel = (I.x[selected_I] / umeter).reshape(-1)   # shape (N_REC_E,)
    Iy_sel = (I.y[selected_I] / umeter).reshape(-1)
    full_path = os.path.join(results_dir,'selected_I_locations.mat')
    sio.savemat(full_path,
        {
            'selected_ids': selected_I,   # neuron indices in E
            'x_um':          Ix_sel,       # x‐positions [μm]
            'y_um':          Iy_sel        # y‐positions [μm]
        }
    )
 

print(np.mean(EE.N_incoming[:]))
print(np.mean(EI.N_incoming[:]))
print(np.mean(IE.N_incoming[:]))
print(np.mean(II.N_incoming[:]))

iteration=1

run(warmup)

for iteration in range(1, nr_iterations + 1):
     
    # Set up monitors
    M_AMPA_P = StateMonitor(E,'I_AE', record=True, dt=1*ms)
    M_GABA_P = StateMonitor(E,'I_GE', record=True, dt=1*ms)
    
    M_E = SpikeMonitor(E, record=selected_E)
    M_I = SpikeMonitor(I, record=selected_I)
    
    # Run the simulation with a single run command for the entire network
    run(duration)
    N = M_AMPA_P.I_AE.shape[1]
   
    def spatial_mask(k,ex,ey,Rx,Ry,LFP_g):
        R_x=130*umeter
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
        LFP_out = np.zeros((N, N_R))
        for k in range(N_R):
            LFP_out[:, k] = spatial_mask(k, E.x, E.y, Rec.x, Rec.y, LFP_g_matrix)
        return LFP_out

    g_absum = np.abs(M_AMPA_P.I_AE[:, :N]) + np.abs(M_GABA_P.I_GE[:, :N])
    LFP_absum  = compute_LFP_for_g(g_absum)

    # PRINT
    dict = {}
    dict['data'] = LFP_absum[:, :]  # I save LFPs every ms
    name='LFP_i{0}_J{1}.mat'.format(iteration,index_print)
    full_path = os.path.join(results_dir, name)
    sio.savemat(full_path, dict)
    
    name='E_spike_times_i{0}_J{1}.mat'.format(iteration-1,index_print)
    full_path = os.path.join(results_dir, name)
    sio.savemat(full_path, {'SpikeTimes': M_E.t[:]})
    name='E_spike_id_i{0}_J{1}.mat'.format(iteration-1,index_print)
    full_path = os.path.join(results_dir, name)
    sio.savemat(full_path,{'SpikeId': M_E.i[:]})
    name='I_spike_times_i{0}_J{1}.mat'.format(iteration-1,index_list)
    full_path = os.path.join(results_dir, name)
    sio.savemat(full_path, {'SpikeTimes': M_I.t[:]})
    name='I_spike_id_i{0}_J{1}.mat'.format(iteration-1,index_print)
    full_path = os.path.join(results_dir, name)
    sio.savemat(full_path,{'SpikeId': M_I.i[:]})
        
    del M_AMPA_P, M_GABA_P, M_E, M_I
   


