from brian2 import *
prefs.codegen.target = 'numpy'  # use the Python fallback (can also use cython)
import numpy as np
import scipy.io as sio
import os
import random
start_scope()

if len(sys.argv) < 4:
    print("-" * 30)
    print("ERROR: Missing argument.")
    print("Usage: python Simulations_const_ext_rate.py <index_print> <nu_{ext}/nu_{theta}> <g>")
    print("Example: python Simulations_dynamics_isotropic.py 1 0.92 3.3")
    print("-" * 30)
    sys.exit(1)
index_print = int(sys.argv[1])
nu_input= float(sys.argv[2])
g_here=float(sys.argv[3])

seed_value = index_print
random.seed(seed_value)
np.random.seed(seed_value) # Ensure NumPy is seeded too
seed(seed_value)


index_list=1
list_index=1
I_rho = 0.0
E_rho = 0.0

results_dir = 'results_isotropic_f5_g{0}_nu{1}'.format(g_here,nu_input)
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

    
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
N_Pr=1500                       # Nr of neurons whose spike trains I will record

# Setting up recording sites (unchanged)
eqs_r='''x : meter
         y : meter
         index : integer (constant)
'''
Rec = NeuronGroup(rows_r * cols_r, eqs_r)
Rec.x = '(i // rows_r) * grid_dist_r - (rows_r-1)/2.0  * grid_dist_r'
Rec.y = '(i % cols_r) * grid_dist_r - (cols_r-1)/2.0  * grid_dist_r'

sio.savemat('Rec_x.mat', {'data': Rec.x[:]})
sio.savemat('Rec_y.mat', {'data': Rec.y[:]})

##################
# External input #
##################
 
nu_baseline =  nu_input*2195 * Hz      # baseline rate
nr_iterations=5
warmup=200*ms
duration=500*ms
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

#############################################################################
#  Synapses
#############################################################################

lam_dist=[380,250,250,250]
EEprob_z=0.15

J_ext_I=0.85
J_ext_E=0.455

g=g_here
J_EE=0.1575
J_EI=0.2625
J_IE=1.07*g*J_EE
J_II=g*J_EI

##################



# Create the excitatory neurons with private external input
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
 

# Create the inhibitory neurons with private external input
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

EE = Synapses(E, E, 'w : volt', on_pre='x_AE_post += w')
EE.connect('i != j',
           p='0.06314*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_pre))**2 + ((y_pre-y_post)/(dist_syn_e_pre))**2 - 2*rho_pre*((x_pre-x_post)/dist_syn_e_pre)*((y_pre-y_post)/dist_syn_e_pre))/sqrt(1-rho_pre*rho_pre))')
EE.w = 'J_EE*mV * tau_mE/tau_rAE'
EE.delay = '0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2)/vel_lat'

EI = Synapses(E, I, 'w : volt', on_pre='x_AI_post += w')
EI.connect(p='0.120777*exp(-sqrt(((x_pre-x_post)/(dist_syn_i_pre))**2 + ((y_pre-y_post)/(dist_syn_i_pre))**2 - 2*rho_post*((x_pre-x_post)/dist_syn_i_pre)*((y_pre-y_post)/dist_syn_i_pre))/sqrt(1-rho_post*rho_post))')
EI.w = 'J_EI*mV * tau_mI/tau_rAI'
EI.delay = '0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2)/vel_lat'

IE = Synapses(I, E, 'w : volt', on_pre='x_GE_post += w')
IE.connect(p='0.12*exp(-sqrt(((x_pre-x_post)/(dist_syn_e_pre))**2 + ((y_pre-y_post)/(dist_syn_e_pre))**2 - 2*rho_pre*((x_pre-x_post)/dist_syn_e_pre)*((y_pre-y_post)/dist_syn_e_pre))/sqrt(1-rho_pre*rho_pre))')
IE.w = 'J_IE*mV * tau_mE/tau_rG'
IE.delay = '0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2)/vel_lat'

II = Synapses(I, I, 'w : volt', on_pre='x_GI_post += w')
II.connect('i != j',
           p='0.12*exp(-sqrt(((x_pre-x_post)/(dist_syn_i_pre))**2 + ((y_pre-y_post)/(dist_syn_i_pre))**2 - 2*rho_pre*((x_pre-x_post)/dist_syn_i_pre)*((y_pre-y_post)/dist_syn_i_pre))/sqrt(1-rho_pre*rho_pre))')
II.w = 'J_II*mV * tau_mI/tau_rG'
II.delay = '0.1*ms + sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2)/vel_lat'


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

    # Run the simulation with a single run command for the entire network
    run(duration)
    N = M_AMPA_P.I_AE.shape[1]
   
    def spatial_mask(k,ex,ey,Rx,Ry,LFP_g):
        R_x=100*umeter
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
    adict = {}
    adict['data'] = LFP_absum[:, :]  # I save LFPs every ms
    name='LFP_i{0}_J{1}.mat'.format(iteration,index_print)
    full_path = os.path.join(results_dir, name)
    sio.savemat(full_path, adict)
 
            
    del M_AMPA_P, M_GABA_P
 
