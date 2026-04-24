% =========================================================================
% Visualization of LFP activity around movement onset
%
% Plots z-scored LFP traces and autocorrelation functions (ACF) for a
% selected trial, separately for pre- and post-movement periods.
%
% Data available on Zenodo: 10.5281/zenodo.19420479
% The script will prompt you to select the directory where data is stored.
%
% =========================================================================

%% load lfp
addpath(fullfile(pwd, '..',  'functions'));

dataDir = uigetdir(pwd, 'Select the folder containing your .mat files');
if ischar(dataDir)
    chanlfp = load(fullfile(dataDir, 'rs1050225_MI_clean_LFP.mat'));
    load(fullfile(dataDir, 'Beh_matrix.mat'));
else
    error('Folder not selected. Script stopped.');
end

% extract channels
fields_lfp   = fieldnames(chanlfp);
units_lfp    = fields_lfp(startsWith(fields_lfp,'lfp'));
num_channels = length(units_lfp);
size_lfp=length(chanlfp.(units_lfp{1}));

%% Parameters
fs=1000;
order=2;
freqrange = [6 80];  % broad filter for visualization
n_chan_visualize=20; % number of LFPs to plot  
selected_trial=49; %select a trial to visualize
 
%% FILTERED
 selected_trial=[57,60];
idx=0;
for trial_nr = selected_trial 
    idx=idx+1;
    initial_time=ceil(beh(trial_nr,1)*fs);
    final_time=ceil(beh(trial_nr,6)*fs)-1;      
    start_move(idx)=ceil(beh(trial_nr,5)*fs)-initial_time;
    end_time(idx)=ceil(beh(trial_nr,6)*fs)-initial_time;
    del_t=final_time-initial_time;
    filtered_lfp=zeros(num_channels,1+del_t);
    for u = 1:num_channels
        lfp_temp = chanlfp.(units_lfp{u});  
        filtered_lfp(u,:)=(buttfilt(lfp_temp(initial_time:final_time),freqrange,fs,'bandpass',order));
    end   
    lfp_zscored=(filtered_lfp-mean(filtered_lfp,2))./(std(filtered_lfp,[],2));
    
    Preparation_LFP=filtered_lfp(:,start_move(idx)-1000:start_move(idx)); %go cue is at 1000 ms
    Execution_LFP=filtered_lfp(:,start_move(idx):end); 
       
    clear times_to_start;
    times_to_start=[1:1:end_time(idx)];
    times_to_start=times_to_start-start_move(idx);
    
    figure
    subplot(3,1,1)
    plot(times_to_start,lfp_zscored(1:n_chan_visualize,:)')  
    xlim([times_to_start(1) times_to_start(end)])
    xlabel('Time (ms)')
    ylabel('LFP amplitude')
    set(gca,'fontsize',12)
    subplot(3,1,2)
    for ll=1:n_chan_visualize
    [xc tt]=xcov(Preparation_LFP(ll,:),Preparation_LFP(ll,:),'normalized');
    plot(tt,xc) %LFPs are sampled every ms
    
    hold on
    end
    xlim([0 150])
    ylim([-1 1])
    xlabel('lag [ms]')
    ylabel('LFPs ACF')
    set(gca,'fontsize',12)

    subplot(3,1,3)
    for ll=1:n_chan_visualize
    [xc tt]=xcov(Execution_LFP(ll,:),Execution_LFP(ll,:),'normalized');
    plot(tt,xc) %LFPs are sampled every ms
    hold on
    end
    xlim([0 150])
    ylim([-1 1])
    xlabel('lag [ms]')
    ylabel('LFPs ACF')
    set(gca,'fontsize',12)

end

 