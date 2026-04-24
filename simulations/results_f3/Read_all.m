% Read_all.m
% Loads simulation output files for a single trial and generates three figures:
%   1. Z-scored LFP traces over time
%   2. Autocorrelation function (ACF) of the LFP, with estimated oscillation frequency
%   3. Raster plot of spiking activity for E (800) and I (200) neurons
% External drive strength is read from the simulation script.


clear

fid = fopen('../run_const_external_rate.sh', 'r');
if fid == -1
    error('Could not open file: %s', '../run_const_external_rate.sh');
end
content = fscanf(fid, '%c');
fclose(fid);

token = regexp(content, 'external_input\s*=\s*([\d.]+)', 'tokens');
if isempty(token)
    error('Could not find external_input in the shell script.');
end
input_drive = str2double(token{1}{1});

T_max = 1000;
skip_times = 0;
N_E = 800;
N_I = 200;
trial_to_plot = 1; %plot one trial
sim_to_plot = 1; %plot one simulation

default_pos = get(0, 'DefaultFigurePosition');
new_width = default_pos(3) * 0.3;
new_height = default_pos(4) * 0.3;

name = sprintf('LFP_i%d_J%d.mat', trial_to_plot, sim_to_plot);
lv = load(name);
LFP_all = (lv.data)';
lfp_zscored = (LFP_all - mean(LFP_all, 2)) ./ std(LFP_all, [], 2);
n_channel = sqrt(size(lfp_zscored, 1));

lfp_zscored(:, 1:skip_times) = [];

% Panel 1: Plot the LFP
figure('Position', [default_pos(1), default_pos(2), new_width, new_height]);
plot(lfp_zscored')
xlim([0 300])
xlabel('Time [ms]')
axis square
ylabel('LFP ampl.')
label_str = sprintf('$\\nu_{\\mathrm{ext}} / \\nu_{\\theta} = %.2f$', input_drive);
set(gca, 'fontsize', 11)

% Panel 2: Plot ACF
figure('Position', [default_pos(1), default_pos(2), new_width, new_height]);
frequencies = [];

for ll = 1:n_channel
    [xc, tt] = xcorr(lfp_zscored(ll, :), 'normalized');
    plot(tt, xc)
    hold on

    pos_lags = tt(tt > 0);
    pos_acf = xc(tt > 0);
    [pks, locs] = findpeaks(pos_acf, 'MinPeakHeight', 0.1);

    if ~isempty(locs)
        period_ms = pos_lags(locs(1));
        freq = 1000 / period_ms;
        frequencies = [frequencies freq];
    end
end

mean_freq = mean(frequencies);
std_freq = std(frequencies);
text_str = sprintf('f = %.1f ± %.1f Hz', mean_freq, std_freq);

ylabel('LFP ACF')
xlim([0 150])
xlabel('lag [ms]')
set(gca, 'fontsize', 11)
axis square

% Panel 3: Raster plot
name = sprintf('E_spike_times_i%d_J%d.mat', trial_to_plot-1, sim_to_plot);
V = load(name);
Time_spikes_e = V.SpikeTimes(1:end);

name = sprintf('E_spike_id_i%d_J%d.mat', trial_to_plot-1, sim_to_plot);
V = load(name);
Id_spikes_e = V.SpikeId(1:end);

name = sprintf('I_spike_times_i%d_J%d.mat', trial_to_plot-1, sim_to_plot);
V = load(name);
Time_spikes_i = V.SpikeTimes(1:end);

name = sprintf('I_spike_id_i%d_J%d.mat', trial_to_plot-1, sim_to_plot);
V = load(name);
Id_spikes_i = V.SpikeId(1:end);

Spike_t = cell(N_E + N_I, 1);

for kk = 1:N_E
    spike_idx = find(Id_spikes_e == kk - 1);
    Spike_t{kk} = Time_spikes_e(spike_idx);
end

for ll = 1:N_I
    spike_idx = find(Id_spikes_i == ll - 1);
    Spike_t{N_E + ll} = Time_spikes_i(spike_idx);
end

raster = zeros(N_E + N_I, T_max);
bins = 0 : .001 : T_max / 1000;
for u = 1:(N_E + N_I)
    raster(u, :) = histcounts(Spike_t{u}, bins);
end

raster = flipud(raster);
[neuron_ids, spike_times] = find(raster);
neuron_ids = flip(neuron_ids);

figure('Position', [default_pos(1), default_pos(2), new_width, new_height]);
scatter(spike_times, neuron_ids, 3, 'k', 'filled')
ylabel('Neuron ID')
yticks([])
xlabel('Time [ms]')
set(gca, 'fontsize', 11)