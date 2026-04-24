% snapshot_spikes.m
% -----------------
% Shows spikes on the 2D array 
% in the trailing 1 ms window for each frame.
% Spikes from the E population in black
% Spikes from the I population in red

folder='Examples'; 
trial=1;
name = fullfile(folder, sprintf('LFP_i%d_J%d.mat',trial,1));
ld=load(name);

figure;plot(ld.data);
xlabel('Time (ms)')
ylabel('LFP')
title('Check LFP')
set(gca,'fontsize',18)

%% Parameters
t0      = 450;% start time [ms]
win_ms  = 3;           %  ms per frame
step_ms = 3;           %  ms step per frame
t_span  = 40;          % total duration [ms]
dt_ms   = 0.1;         % Brian defaultclock (0.1 ms)
time_list = t0:step_ms:(t0 + t_span);
T = numel(time_list);
% Apply margin 
mx = 400; % µm 
my =400; % µm 

%% File names
e_loc_file   = 'selected_E_locations.mat';
i_loc_file   = 'selected_I_locations.mat';
e_times_file = 'E_spike_times_i0_J1.mat';
e_id_file    = 'E_spike_id_i0_J1.mat';
i_times_file = 'I_spike_times_i0_J1.mat';
i_id_file    = 'I_spike_id_i0_J1.mat';

%% Load data
E_loc = load(e_loc_file);   % fields: selected_ids, x_um, y_um
I_loc = load(i_loc_file);
E_sp  = load(e_times_file); % field: SpikeTimes (seconds)
E_id  = load(e_id_file);    % field: SpikeId
I_sp  = load(i_times_file);
I_id  = load(i_id_file);

% Convert spike times to milliseconds
E_times = double(E_sp.SpikeTimes) * 1e3;  % ms
I_times = double(I_sp.SpikeTimes) * 1e3;  % ms

% Optional: snap to 0.1 ms grid to avoid fp jitter
E_times = round(E_times/dt_ms)*dt_ms;
I_times = round(I_times/dt_ms)*dt_ms;

% Precompute for mapping IDs->positions
sel_ids_E_d = double(E_loc.selected_ids);
sel_ids_I_d = double(I_loc.selected_ids);

% After loading E_loc/I_loc:
Xall = [E_loc.x_um(:); I_loc.x_um(:)];
Yall = [E_loc.y_um(:); I_loc.y_um(:)];

xmin = min(Xall); xmax = max(Xall);
ymin = min(Yall); ymax = max(Yall);

xin = [xmin + mx, xmax - mx];
yin = [ymin + my, ymax - my];

 
%% Plot
figure('Color','w');
for k = 1:T
    t = time_list(k);

    % spikes in half-open window  
    e_mask = (E_times > t - win_ms) & (E_times <= t);
    i_mask = (I_times > t - win_ms) & (I_times <= t);

    e_active_ids_d = double(E_id.SpikeId(e_mask));
    i_active_ids_d = double(I_id.SpikeId(i_mask));

    [liaE, locE] = ismember(e_active_ids_d, sel_ids_E_d);
    xE = E_loc.x_um(locE(liaE));  yE = E_loc.y_um(locE(liaE));

    [liaI, locI] = ismember(i_active_ids_d, sel_ids_I_d);
    xI = I_loc.x_um(locI(liaI));  yI = I_loc.y_um(locI(liaI));
   
    ax = subplot(ceil(T/10), 10, k); hold on;
    scatter(xE, yE, 5, 'k', 'filled');
    scatter(xI, yI, 5, 'r', 'filled');
    xlim(xin); ylim(yin);
    axis square
    box(ax,'on');                 % <- adds the top & right borders
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 1;

    xticklabels([]); yticklabels([]);
    title(sprintf('t = (%.0f, %.0f] ms', t - win_ms - t0, t - t0));
    set(ax,'LooseInset',[0 0 0 0]);  % (optional; beware it can clip borders when exporting)

end

