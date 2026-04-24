function results = Mean_field_rates_IC(param,initial_nu)

    lb=[0.000001,0.000001];
    ub=[500, 500];
    errfcn = @(x) Compute_nu(x,param); 
    options = optimoptions(@fmincon,'OptimalityTolerance',1.0000e-08);
    [x,fval] = fmincon(errfcn,initial_nu,[],[],[],[],lb,ub);
    
    %%
   
    nu_e=x(1);
    nu_i=x(2);

    mu_i_A=param.tau_mi *(param.Ji_ext*param.nu_ext+param.wi*param.C_IE*nu_e);
    mu_i_G=param.tau_mi *(param.wi*param.C_II* param.gi*nu_i);
    mu_e_A=param.tau_me *(param.Je_ext*param.nu_ext+param.we*param.C_EE*nu_e);
    mu_e_G=param.tau_me *(param.we*param.C_EI* param.ge*nu_i);
 
    sig_i_A=param.tau_mi *(param.Ji_ext^2*param.nu_ext+param.wi^2*param.C_IE*nu_e);
    sig_i_G=param.tau_mi *param.wi^2*param.gi^2*param.C_II*nu_i;  
    sig_e_A=param.tau_me *(param.Je_ext^2*param.nu_ext+param.we^2*param.C_EE*nu_e);
    sig_e_G=param.tau_me *param.we^2*param.ge^2*param.C_EI*nu_i;  
 
    t_syn_i=(sig_i_A+sig_i_G)/((sig_i_A/param.tilde_AMPA_i)+(sig_i_G/param.tilde_GABA_i));
    t_syn_e=(sig_e_A+sig_e_G)/((sig_e_A/param.tilde_AMPA_e)+(sig_e_G/param.tilde_GABA_e));
    
    results=zeros(6,1);
    results(1)=mu_e_A-mu_e_G; %mu_E
    results(2)=mu_i_A-mu_i_G; %mu_I
    results(3)=sig_e_A+sig_e_G; %sigma_E squared
    results(4)=sig_i_A+sig_i_G; %sigma_I squared
    results(5)=nu_e;
    results(6)=nu_i; 
    results(7)=t_syn_e; 
    results(8)=t_syn_i;


end