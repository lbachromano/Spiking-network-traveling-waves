
clear param
alpha=2.065;

%%
param.L_ee=1; %Homogeneous network; these parameters won't be used in this folder
param.L_ei=1;
param.L_ie=1;
param.L_ii=1;
%%
param.tau_me=20e-3;  
param.tau_mi=10e-3;

param.tau_L_eA=2.57e-3;
param.tau_L_eG=1.725e-3;
param.tau_L_iA=1.725e-3;
param.tau_L_iG=1.725e-3;

param.taurpe =2.e-3;
param.taurpi =1.e-3;

param.tau_d_Ae=3.5e-3;
param.tau_d_Ai=3.5e-3;
param.tau_r_Ae=0.7e-3;
param.tau_r_Ai=0.7e-3;

param.tau_d_Gi=18e-3;
param.tau_r_Gi=1e-3;
param.tau_d_Ge=param.tau_d_Gi;
param.tau_r_Ge=param.tau_r_Gi;

param.tilde_AMPA_e=param.tau_d_Ae+param.tau_r_Ae;
param.tilde_GABA_e=param.tau_d_Ge+param.tau_r_Ge;
param.tilde_AMPA_i=param.tau_d_Ai+param.tau_r_Ai;
param.tilde_GABA_i=param.tau_d_Gi+param.tau_r_Gi;



%% N

param.C_EE= 800;
param.C_IE= 800;
param.C_EI= 200;
param.C_II= 200;

%% mV

param.theta =18;
param.hvr= 11;


param.gi= 3.625;
param.ge= 1.07.*param.gi;
param.we=1.5*0.1050;
param.wi=1.5*0.1750;

param.Ji_ext = 0.85;
param.Je_ext = 0.455;
 
param.nutheta=param.C_EE*param.theta/(param.C_EE*param.tau_me*param.we);
param.nu_ext=param.nutheta*0.1;


%% nu_theta from analytic formula
% excitatory population (a = e)
tau_m_a    = param.tau_me;
tau_syn_a  = param.tilde_AMPA_e;   % or whichever synaptic time const. you want
J_a_ext    = param.Je_ext;
theta      = param.theta;

num = (alpha/2) * J_a_ext * sqrt(tau_syn_a) + ...
      sqrt( ((alpha/2) * J_a_ext * sqrt(tau_syn_a))^2 + ...
             4 * tau_m_a * J_a_ext * theta );

den = 2 * tau_m_a * J_a_ext;
nt_alpha_E = (num/den)^2;

tau_m_a    = param.tau_mi;
tau_syn_a  = param.tilde_AMPA_i;   % or whichever synaptic time const. you want
J_a_ext    = param.Ji_ext;
theta      = param.theta;

num = (alpha/2) * J_a_ext * sqrt(tau_syn_a) + ...
      sqrt( ((alpha/2) * J_a_ext * sqrt(tau_syn_a))^2 + ...
             4 * tau_m_a * J_a_ext * theta );

den = 2 * tau_m_a * J_a_ext;
nt_alpha_I = (num/den)^2;



param.nt_alpha = 0.8*nt_alpha_E+0.2*nt_alpha_I;
param.nutheta=param.nt_alpha;
