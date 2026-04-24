% Load and plot external input rate

nu_th=2195;
data = load('examples_original_code/external_input_rate.mat');
fs=0.1;
skip_times=60/fs; %ms, in line with the other processing scripts
time_to_plot=data.time_ms(1:end-skip_times); %in line with the other processing scripts
rate_to_plot = data.nu_global(skip_times+1:end);

figure;
plot(time_to_plot, rate_to_plot./nu_th, 'k-', 'LineWidth', 2);
xlim([0 1800])
xlabel('Time (ms)');
ylabel('Rate (Hz)');
set(gca,'fontsize',16)
ylim([0.7 1.7])
