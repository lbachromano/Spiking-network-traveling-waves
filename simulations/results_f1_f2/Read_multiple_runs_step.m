% Loads and preprocesses simulated Local Field Potential (LFP) data by applying 
% a bandpass Butterworth filter and z-score normalization across multiple channels. 
% It then visualizes the results by plotting the filtered time-series signals 
% and comparing the auto-covariance functions of the signal during defined 
% "preparation" and "execution" phases.

clear
addpath(fullfile(pwd, '..', '..', 'functions'));

 
%parameters
folder_name='examples_paper';   %used as example here, replace
%folder_name='examples';
skip_times=60; %skip the first 100 ms of the simulations
Ex_T=1200; %after having skipped 100
T_end=1600;

freqrange = [6 250]; %gross filtering just not to see the step function
                     % The data are filtered in the same way
%parameters for filter
fs=1000;
order=2; 

for trial_number=1:10 %one selected example

  %used as example here, replace with the desired data list

    name = fullfile(folder_name, sprintf('LFP_absSum_i%d_J1.mat', trial_number));
    LFP_from_brian=load(name);
    lfp_zscored=(LFP_from_brian.data)';
    lfp_zscored(:,1:skip_times)=[];
    total_T=size(lfp_zscored,2);
    n_channel=size(lfp_zscored,1);
    
    filtered_lfp=NaN(size(lfp_zscored));
    for u=1:n_channel
        filtered_lfp(u,:)=(buttfilt(lfp_zscored(u,:),freqrange,fs,'bandpass',order));
    end
    lfp_zscored=(filtered_lfp-mean(filtered_lfp,2))./(std(filtered_lfp,[],2)); 

    prep_LFP=filtered_lfp(:,1:Ex_T);
    exec_LFP=filtered_lfp(:,Ex_T:T_end);


    ttt=1:total_T;
    figure('Position', [100, 100, 400, 400])  % [left, bottom, width, height] in pixels
    subplot(3,1,1)
    plot(ttt,lfp_zscored(1:n_channel,:)')
    xlim([0, T_end]) %show the first 1800 ms for comparison with data
    ylabel('LFP amplitude')
    xlabel('Time (ms)')
    
    subplot(3,1,2)
    for ll=1:n_channel
        [xc tt]=xcov(prep_LFP(ll,:),prep_LFP(ll,:),'normalized');
        plot(tt,xc)  
        hold on
    end
    ylabel('LFP ACF')
    xlabel('Time lag (ms)')
    xlim([0 150])
    
    subplot(3,1,3)
    for ll=1:n_channel
        [xc tt]=xcov(exec_LFP(ll,:),exec_LFP(ll,:),'normalized');
        plot(tt,xc)  
        hold on
    end
    ylabel('LFP ACF')
    xlabel('Time lag (ms)')
    xlim([0 150])
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 6 6]);  % 6x6 inches square

end

% saveas(gcf, sprintf('LFP_plot_trial%d.pdf', trial_number));