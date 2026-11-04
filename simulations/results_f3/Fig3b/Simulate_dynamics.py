from brian2 import *
prefs.codegen.target = 'numpy'  # use the Python fallback (can also use cython)
import numpy as np
import scipy.io as sio
import os
import random


# check there is an argument after the script name
if len(sys.argv) < 3:
    print("Usage: python myscript.py <dx_e>")
    sys.exit(1)

dx_e = int(sys.argv[1])   # convert from string to float
print("dx_e =", dx_e)

random_seed = int(sys.argv[2])

print("dx_e =", dx_e)
print("seed =", random_seed)

# Set random seeds
seed(random_seed)           # Brian2
np.random.seed(random_seed) # NumPy
random.seed(random_seed)    # Python random


index_list=dx_e

List_syn_dist = [0.8,0.833,0.866,0.9,0.933,0.966,1,1.33]
 
JXE = 2195*List_syn_dist[index_list-1]

###############
#  Size grid  #
###############
Size_x=2050;
Size_y = Size_x
Volume = Size_x * Size_y * umeter ** 2

scale_factor = (2050 / 2600)**2

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

print("defaultclock.dt =", defaultclock.dt)
print("defaultclock.dt in ms =", defaultclock.dt / ms)

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
rec_indices_E = np.random.choice(N_E, N_REC_E, replace=False)
rec_indices_I = np.random.choice(N_I, N_REC_I, replace=False)
 
PE = Synapses(P_E, E, 'w : volt', on_pre='x_AE_post += w')
PE.connect(j='i')
PE.w = 'J_ext_E*mV *tau_mE/tau_rAE '
PE.delay = 0*ms

PI = Synapses(P_I, I, 'w : volt', on_pre='x_AI_post += w')
PI.connect(j='i')
PI.w = 'J_ext_I*mV *tau_mI/tau_rAI'
PI.delay = 0*ms



# E to E
EE = Synapses(E, E, 'w : volt', on_pre='x_AE_post += w')
EE.connect('i != j', p='0.9873*0.8273*2*0.0043*0.7511/0.6217')
EE.w ='J_EE*mV *tau_mE/tau_rAE'
EE.delay = '2.57*ms'

# E to I
EI = Synapses(E, I, 'w : volt', on_pre='x_AI_post += w')
EI.connect(p='0.9873*0.8273*2*0.0043*0.7511/0.6217')
#EI.connect('i != j', p='0.5 *exp(-(((x_pre-x_post)/dist_syn_i_pre)**2 + ((y_pre-y_post)/dist_syn_i_pre)**2 - 2*rho_pre*(x_pre-x_post)*(y_pre-y_post)/(dist_syn_i_pre*dist_syn_i_pre))/(2*(1 - rho_pre*rho_pre)))')
#EI.w = 2*np.random.lognormal( -0.84,0.98   ,size(EI.w)) *mV  *tau_mI/tau_rAI     #logn opt 5
EI.w ='J_EI*mV  *tau_mI/tau_rAI'
EI.delay = '1.725*ms'
# I to E
IE = Synapses(I, E, 'w : volt', on_pre='x_GE_post += w')
IE.connect(p='0.98*0.8240*2*0.013*0.7511/(0.6217*3)')
#IE.connect('i != j', p='0.5 *exp(-(((x_pre-x_post)/dist_syn_e_pre)**2 + ((y_pre-y_post)/dist_syn_e_pre)**2 - 2*rho_pre*(x_pre-x_post)*(y_pre-y_post)/(dist_syn_e_pre*dist_syn_e_pre))/(2*(1 - rho_pre*rho_pre)))')
#IE.w =2*np.random.lognormal(-0.40,1.37,size(IE.w))*mV *tau_mE/tau_rG   #logn opt 6
IE.w ='J_IE*mV *tau_mE/tau_rG'
IE.delay = '1.725*ms'
# I to I
II = Synapses(I, I, 'w : volt', on_pre='x_GI_post += w')
II.connect('i != j', p='0.98*0.8240*2*0.013*0.7511/(0.6217*3)')
#II.w = 2* np.random.lognormal(-0.08,1.47,size(II.w))*mV *tau_mI/tau_rG   #logn opt 6
II.w ='J_II*mV *tau_mI/tau_rG'
II.delay = '1.725*ms'

print(np.mean(EE.N_incoming[:]))
print(np.mean(EI.N_incoming[:]))
print(np.mean(IE.N_incoming[:]))
print(np.mean(II.N_incoming[:]))

# Initial transient run


    
run(duration1)

#  RUN
for iteration in range(1, number_iterations+1):

    # Population rate monitors
    E_pop = PopulationRateMonitor(E)
    I_pop = PopulationRateMonitor(I)
    
    # Run the simulation
    run(duration)

   
    # PRINT
    adict = {}

    name = 'RATE_E_i{0}_J{1}_trial{2}.mat'.format(iteration-1,index_list-1,random_seed)
    sio.savemat(name, {'data':E_pop.smooth_rate(window='flat', width=5*ms)/Hz})
    name = 'RATE_I_i{0}_J{1}_trial{2}.mat'.format(iteration-1,index_list-1,random_seed)
    sio.savemat(name, {'data':I_pop.smooth_rate(window='flat', width=5*ms)/Hz})

