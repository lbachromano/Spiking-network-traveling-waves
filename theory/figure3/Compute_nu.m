function f= Compute_nu(all_nu,param)
 
    alpha=sqrt(2)*abs(zeta(1/2));
    Vr=param.hvr;
    theta=param.theta;
    
    xnuE=all_nu(1);
    xnuI=all_nu(2);
   
    mu_i_A=param.tau_mi *(param.Ji_ext*param.nu_ext+param.wi*param.C_IE*xnuE);
    mu_i_G=param.tau_mi *(param.wi*param.C_II* param.gi*xnuI);
    mu_e_A=param.tau_me *(param.Je_ext*param.nu_ext+param.we*param.C_EE*xnuE);
    mu_e_G=param.tau_me *(param.we*param.C_EI* param.ge*xnuI);    
    mu_i=mu_i_A-mu_i_G;
    mu_e=mu_e_A-mu_e_G;
    
    sig_i_A=param.tau_mi *(param.Ji_ext^2*param.nu_ext+param.wi^2*param.C_IE*xnuE);
    sig_i_G=param.tau_mi *param.wi^2*param.gi^2*param.C_II*xnuI;  
    sig_e_A=param.tau_me *(param.Je_ext^2*param.nu_ext+param.we^2*param.C_EE*xnuE);
    sig_e_G=param.tau_me *param.we^2*param.ge^2*param.C_EI*xnuI;  
    sigma_i=sig_i_A+sig_i_G;
    sigma_e=sig_e_A+sig_e_G;
    
    
    t_syn_i=(sig_i_A+sig_i_G)/((sig_i_A/param.tilde_AMPA_i)+(sig_i_G/param.tilde_GABA_i));
    t_syn_e=(sig_e_A+sig_e_G)/((sig_e_A/param.tilde_AMPA_e)+(sig_e_G/param.tilde_GABA_e));

    E_lower=(Vr-mu_e)/sqrt(sigma_e)+alpha*sqrt(t_syn_e/param.tau_me)/2;
    E_upper=(theta-mu_e)/sqrt(sigma_e)+alpha*sqrt(t_syn_e/param.tau_me)/2;
    
    I_lower=(Vr-mu_i)/sqrt(sigma_i)+alpha*sqrt(t_syn_i/param.tau_mi)/2;
    I_upper=(theta-mu_i)/sqrt(sigma_i)+alpha*sqrt(t_syn_i/param.tau_mi)/2;

    fun1 = @(y) erfcx(-y); 
    nuE_new = param.taurpe+param.tau_me.*sqrt(pi).*integral(fun1,E_lower,E_upper);
    nuE_new=1./nuE_new;
    nuI_new = param.taurpi+param.tau_mi.*sqrt(pi).*integral(fun1,I_lower,I_upper);
    nuI_new=1./nuI_new;
 
    f= (nuI_new - all_nu(2)).^2+(nuE_new - all_nu(1)).^2;

end

 