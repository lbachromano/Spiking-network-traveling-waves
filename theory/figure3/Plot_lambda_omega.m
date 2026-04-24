%plots the real part of the eigenvalue at g=3.625 for varying external field, 
% and the imaginary part compared to frequency of oscillations from simulations 


clear
Set_parameters;

figure;
lf=load(sprintf('Oscillation_LFP_f.mat'));
all_y=lf.mean_freq
all_std=lf.sem_freq
all_x=lf.input_drive;
errorbar(all_x,all_y,all_std,'o')


%% Load results from linear stability analysis
myname=sprintf('Linear_stability_results.mat');
lc=load(myname);
all_param=[];
all_minima=[];
for i=1:length(lc.local_minima_all)
    min_h=lc.local_minima_all{i};
    nr_minima=size(min_h,1);
    params_h=lc.list_varying_param(i).*ones(nr_minima,1);
    all_param=[all_param; params_h];
    all_minima=[all_minima; min_h];
end
 

%% Plot


plot(all_param , abs(all_minima(:,2)), '.k', 'linewidth',2);
hold on
errorbar(all_x,all_y,all_std,'.','markersize',40)
ylabel('\omega [Hz]')
xlabel('\nu_{ext}/\nu_{\theta}')
set(gca,'fontsize',18)
xlim([0.7 1.4])
