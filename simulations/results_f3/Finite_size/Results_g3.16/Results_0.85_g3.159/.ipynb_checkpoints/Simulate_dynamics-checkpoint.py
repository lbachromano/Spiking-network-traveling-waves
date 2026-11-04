from brian2 import *
prefs.codegen.target = 'numpy'  # use the Python fallback (can also use cython)
import numpy as np
import scipy.io as sio
import os



# check there is an argument after the script name
if len(sys.argv) < 3:
    print("Usage: python myscript.py <dx_e>")
    sys.exit(1)

index_list= int(sys.argv[1])   # convert from string to float
print("dx_e =", index_list)

random_seed= int(sys.argv[2])   # convert from string to float
print("random_seed =", random_seed)

Size_list = [500, 750, 1000, 1250, 1500,1750, 2000, 2250, 2500]
Size_x=Size_list[index_list-1] #2050;

JXE = 2195*0.85 #0.85

###############
#  Size grid  #
###############

Size_y = Size_x
Volume = Size_x * Size_y * umeter ** 2

Size_ref = 2050.0
scale_factor = (Size_ref / Size_x)**2

###################
#  'Utah arrays'  #
###################
# At least 150*umeter from the array to the border of the grid
grid_dist_r = 400 * umeter
grid_nr = int((Size_x - 400 + grid_dist_r / umeter) // (grid_dist_r / umeter))  # // floor division
rows_r, cols_r = grid_nr, grid_nr

################
#  Nr neurons  #
################
Total_nr_neurons = 28000 * (Volume / mm ** 2)
N_E = int(Total_nr_neurons * 0.8)  # Pyramidal cells
N_I = int(Total_nr_neurons * 0.2)  # Interneurons
N_R = grid_nr ** 2                 # Recording sites
N_Pr = 1500                        # Nr of neurons whose spike trains I will record

# Position of neurons to be recorded
eqs_r = '''
    x : meter
    y : meter
    index : integer (constant)
'''
Rec = NeuronGroup(rows_r * cols_r, eqs_r)
Rec.x = '(i // rows_r) * grid_dist_r - (rows_r-1)/2.0  * grid_dist_r'
Rec.y = '(i % cols_r) * grid_dist_r - (cols_r-1)/2.0  * grid_dist_r'

##################
# External input #
##################
nu_signal_exc = JXE * Hz
nu_signal_inh = JXE * Hz

sigma_n = 0.12 / ms
tau_n = 150 * ms

##############
# Duration   #
##############
number_iterations = 4
iteration = 2
duration1 = 100 * ms
duration = 500 * ms
duration_ext = duration1 + duration*(number_iterations)

##################
#  Single neuron #
##################
rho_d = 0.6
V_th = 18 * mV
v_reset = 11 * mV

tau_mE = (20) * ms
tau_mI = (10) * ms
tau_ref_E = 2 * ms
tau_ref_I = 1 * ms

tau_rAE = 0.7 * ms
tau_rAI = 0.7 * ms
tau_dAE = 3.5 * ms
tau_dAI = 3.5 * ms

tau_rG = 1 * ms
tau_dG = 18 * ms

#############################################################################
#  Synapses
#############################################################################

lam_dist = [380, 250, 250, 250]
EEprob_z = 0.15

J_ext_I = 0.85
J_ext_E = 0.455

g = 3.159
J_EE = 0.1575
J_EI = 0.2625
J_IE = g*1.07 * J_EE
J_II = g * J_EI

##################

start_scope()

P_E = PoissonGroup(N_E, rates=nu_signal_exc)
P_I = PoissonGroup(N_I, rates=nu_signal_exc)

# Excitatory neurons
eqs_e = '''
    x : meter
    y : meter
    dist_syn_e : meter
    dist_syn_i : meter
    rho : 1
    dv_E/dt = (-v_E + I_AE - I_GE)/tau_mE : volt (unless refractory)
    dI_AE/dt   = (-I_AE + x_AE)/tau_dAE : volt
    dI_GE/dt   = (-I_GE + x_GE)/tau_dG : volt
    dx_AE/dt   = (-x_AE)/tau_rAE : volt
    dx_GE/dt   = (-x_GE)/tau_rG : volt
    index : integer (constant)
'''
E = NeuronGroup(N_E, eqs_e, threshold='v_E>V_th', reset='v_E=v_reset', refractory=tau_ref_E, method='rk4')
E.I_AE = 0 * volt
E.I_GE = 0 * volt
E.x_AE = 0 * volt
E.x_GE = 0 * volt

E.x = np.random.uniform(low=-Size_x / 2, high=Size_x / 2, size=(N_E)) * umeter
E.y = np.random.uniform(low=-Size_y / 2, high=Size_y / 2, size=(N_E)) * umeter
E.dist_syn_e = lam_dist[0] * umeter
E.dist_syn_i = lam_dist[1] * umeter
E.rho = rho_d

# Inhibitory neurons
eqs_i = '''
    x : meter
    y : meter
    cc : 1
    dist_syn_e : meter
    dist_syn_i : meter
    rho : 1
    dv_I/dt = (-v_I +  I_AI - I_GI)/tau_mI : volt (unless refractory)
    dI_AI/dt   = (-I_AI + x_AI)/tau_dAI : volt
    dI_GI/dt   = (-I_GI + x_GI)/tau_dG : volt
    dx_AI/dt   = (-x_AI)/tau_rAI : volt
    dx_GI/dt   = (-x_GI )/tau_rG : volt
    index : integer (constant)
'''
I = NeuronGroup(N_I, eqs_i, threshold='v_I>V_th', reset='v_I=v_reset', refractory=tau_ref_I, method='rk4')
I.I_AI = 0 * volt
I.I_GI = 0 * volt
I.x_AI = 0 * volt
I.x_GI = 0 * volt
I.x = np.random.uniform(low=-Size_x / 2, high=Size_x / 2, size=(N_I)) * umeter
I.y = np.random.uniform(low=-Size_y / 2, high=Size_y / 2, size=(N_I)) * umeter
I.cc = 2
I.dist_syn_e = lam_dist[2] * umeter
I.dist_syn_i = lam_dist[3] * umeter
I.rho = rho_d

# Increase the number of neurons to record for firing rate calculations
N_REC_E = 9000  # Number of excitatory neurons to record for rates
N_REC_I = 9000  # Number of inhibitory neurons to record for rates

 
# Randomly select neurons to record
  
 
PE = Synapses(P_E, E, 'w : volt', on_pre='x_AE_post += w')
PE.connect(j='i')
PE.w = 'J_ext_E*mV *tau_mE/tau_rAE '
PE.delay = 0*ms

PI = Synapses(P_I, I, 'w : volt', on_pre='x_AI_post += w')
PI.connect(j='i')
PI.w = 'J_ext_I*mV *tau_mI/tau_rAI'
PI.delay = 0*ms



# Reference connection probabilities at Size_x = 2050 um
p_E = 0.0084838996
p_I = 0.0084551686

# Scale probability to preserve mean in-degree
p_E *= scale_factor
p_I *= scale_factor

# E to E
EE = Synapses(E, E, 'w : volt', on_pre='x_AE_post += w')
EE.connect('i != j', p=p_E)
EE.w = 'J_EE*mV *tau_mE/tau_rAE'
EE.delay = '2.57*ms'

# E to I
EI = Synapses(E, I, 'w : volt', on_pre='x_AI_post += w')
EI.connect(p=p_E)
EI.w = 'J_EI*mV *tau_mI/tau_rAI'
EI.delay = '1.725*ms'

# I to E
IE = Synapses(I, E, 'w : volt', on_pre='x_GE_post += w')
IE.connect(p=p_I)
IE.w = 'J_IE*mV *tau_mE/tau_rG'
IE.delay = '1.725*ms'

# I to I
II = Synapses(I, I, 'w : volt', on_pre='x_GI_post += w')
II.connect('i != j', p=p_I)
II.w = 'J_II*mV *tau_mI/tau_rG'
II.delay = '1.725*ms'

print(np.mean(EE.N_incoming[:]))
print(np.mean(EI.N_incoming[:]))
print(np.mean(IE.N_incoming[:]))
print(np.mean(II.N_incoming[:]))

# Initial transient run


    
run(duration1)

#  RUN
for iteration in range(1, number_iterations+1):
    # DELETE AND RECREATE STATE MONITORS
    M_AMPA_P = StateMonitor(E, 'I_AE', record=True, dt=1*ms)
    M_GABA_P = StateMonitor(E, 'I_GE', record=True, dt=1*ms)

    # Population rate monitors
    E_pop = PopulationRateMonitor(E)
    I_pop = PopulationRateMonitor(I)

    # Run the simulation
    run(duration)

    N = M_AMPA_P.I_AE.shape[1]
   
    def spatial_mask(k,ex,ey,Rx,Ry,LFP_g):
        R_x=150*umeter
        R_r=R_x*1.85
        f_x=(R_x/umeter)**2/sqrt(R_x/umeter)
        temp_dist=sqrt((ex[:]-Rx[k])**2 + (ey[:]-Ry[k])**2)
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
    name = 'LFP_i{0}_J{1}_{2}.mat'.format(iteration - 1, index_list-1,random_seed)
    sio.savemat(name, adict)
 
 
    name = 'RATE_E_i{0}_J{1}_{2}.mat'.format(iteration-1,index_list-1,random_seed)
    sio.savemat(name, {'data':E_pop.smooth_rate(window='flat', width=5*ms)/Hz})
    name = 'RATE_I_i{0}_J{1}_{2}.mat'.format(iteration-1, index_list-1, random_seed)
    sio.savemat(name, {'data':I_pop.smooth_rate(window='flat', width=5*ms)/Hz})
 