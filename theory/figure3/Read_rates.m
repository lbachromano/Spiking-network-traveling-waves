clear 
Set_parameters;
%% Read the results from simulations (can be reproduce using code 
% simulations/run_const_external_rate.sh with homogeneous=1 and varying
% external inputs as in the variable input_drive below

foldername='Rates_from_simulations';   
input_drive=[1571,1714,1857,1999,2142,2285,2428,2571,2714,2857,2999,3142]./param.nutheta;
for index_J=1:length(input_drive)
    ttl=sprintf('%s/FiringRates_J%d.mat',foldername,index_J);
    rate(index_J)=load(ttl).mean_rate_E;
    st_rate(index_J)=load(ttl).std_rate_E;    
    rate_i(index_J)=load(ttl).mean_rate_I;
    st_rate_i(index_J)=load(ttl).std_rate_I;
end

%% Computing firing rates from mean field

list_varying_param=2195.*linspace(0.7, 1.4,15);
N_pt=length(list_varying_param);
initial_nu=[0.5,5];
for j=1:N_pt 
    
    param.nu_ext=list_varying_param(j);
    soluz_MF=Mean_field_rates_IC(param,initial_nu);
    initial_nu=[soluz_MF(5),soluz_MF(6)];
    nu_I(j,1) =  soluz_MF(6);
    nu_E(j,1) =  soluz_MF(5);
 
end

 
%% Plot

 colors=[     0    0.3470    0.5410;
    0.7500    0.2250    0.0980];

figure
errorbar(input_drive,rate,st_rate,'.','color',colors(1,:),'markersize',30)
hold on
errorbar(input_drive,rate_i,st_rate_i,'.','color',colors(2,:),'markersize',30)
hold on
plot(list_varying_param./param.nutheta,nu_E,'color',colors(1,:),'linewidth',3)
hold on
plot(list_varying_param./param.nutheta,nu_I,'color',colors(2,:),'linewidth',3)
xlabel('\nu_{ext} / \nu_{\theta}')
ylabel('Rate [Hz]')
legend('E','I')
set(gca,'fontsize',20)
xlim([0.7 1.4])
 