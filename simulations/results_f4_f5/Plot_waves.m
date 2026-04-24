  % =========================================================================
% TRAVELLING WAVE DETECTION IN LFP ARRAY FROM SIMULATIONS 
% =========================================================================
%
% DESCRIPTION:
%   This script detects and characterises planar travelling waves in
%   local field potential (LFP) from simulations (constant external input).
%   Detection is based on two criteria applied at each time frame:
%     (1) A high Rayleigh Z statistic, indicating that spatial phase
%         gradients are coherently oriented (consistent with a plane wave).
%     (2) A low plane-fit residual, indicating that the instantaneous
%         phase map is well described by a planar wavefront.
%   Thresholds for both criteria are derived empirically from a
%   conservative spatial surrogate null distribution (block-shuffled
%   phase maps), avoiding parametric assumptions.
%
% INPUTS:
%   Outputs of the simulation Simulations_const_ext_rate.py in the parent
%   directory (e.g. LFP_i1_J1.mat)
%
% OUTPUTS:
%   wavelengths_nu<X>.mat       - Per-wave spatial wavelength estimates (cm)
%   duration_nu<X>.mat          - Per-wave duration (samples)
%   where <X> is the value of nu_baseline read from the Python file Simulations_const_ext_rate.py.
%
% MAIN PROCESSING STEPS:
%   1. Load and z-score LFP data; map channels onto 2-D array geometry.
%   2. Identify the dominant LFP frequency using Welch PSD of the
%      spatial-mean signal; bandpass-filter around this frequency (±3 Hz).
%   3. Compute instantaneous phase via Hilbert transform; spatially unwrap.
%   4. At each time frame, compute:
%        - Rayleigh Z statistic from unit-normalised phase gradient vectors.
%        - RMS residual of the best-fit plane to the phase map.
%        - Spatial phase gradient magnitude (used for wave speed estimation).
%   5. Generate surrogate phase maps by spatial block-shuffling and compute
%      the same statistics to build a null distribution.
%   6. Detect wave epochs: contiguous frames where Z > threshold_Z (99.9999th
%      percentile of surrogate Z) AND residual < threshold_e (20th percentile
%      of surrogate residual), lasting at least 10 ms.
%   7. For each detected wave, estimate:
%        - Mean propagation speed (cm/s): |d phi/dt| / |grad phi|.
%        - Propagation direction (radians): from plane-fit wave-vector,
%          oriented by the sign of the mean temporal phase derivative.
%        - Spatial wavelength (cm): 2*pi / |k|, where k is the plane-fit
%          wave-vector magnitude.
%
% PARAMETERS (set at the top of the script):
%   max_nr_waves_to_plot - maximum number of waves to plot
%   threshol_input  - Percentile for Rayleigh Z threshold (default: 99.9999)
%   th_e            - Percentile for residual threshold (default: 20)
%   width_filter    - Half-bandwidth of bandpass filter in Hz (default: 3)
%   spatial_resolution - Electrode spacing in cm (default: 0.0141 cm)
%   K_surr_per_frame   - Surrogate samples generated per time frame (default: 5)
%   use_modes          - Surrogate type: 'spatialblock', 'phase', or 'rowcolshift'
%   block_sz           - Block size for spatial block shuffling (default: 1)

% =========================================================================


rng(1);

%PARAMETERS
max_nr_waves_to_plot=4;
folder='Examples'; %read the 50 trials in the folder Examples; change as desired
number_of_trials_per_simulations=1;
number_independent_simulations=5; 

width_filter=3;
threshold_input=99.9999 % percentile; High Z threshold from shuffled null
th_e=20; %percentile;  low residual threshold from shuffled
which_wave=1;
K_surr_per_frame=5;
spatial_resolution = 0.0141;   % cm (i.e. 0.4 mm electrode spacing)
dx = spatial_resolution;      % now dx = 0.04 cm
%this method will give me a threshold around 0.5 for the PGD
use_modes = {'spatialblock'};  % or {'phase'} or 'rowcolshift'
block_sz=1; %can be 1,2,5,10
 
addpath(fullfile(pwd, '..', '..', 'functions'));

Concatenate_trials(folder,number_of_trials_per_simulations,number_independent_simulations); 
chanlfp=load('Noisy_LFP_anis.mat');
lfp_temp_l = chanlfp.prep_LFP';
number_trials=1;
L_array_x=sqrt(size(lfp_temp_l,2));
L_array_y=sqrt(size(lfp_temp_l,2));


channel_map=zeros(L_array_y,L_array_x);
for i = 1:L_array_x
    channel_map(:,i) = ((i-1)*L_array_y + 1) : (i*L_array_y);  % ascending down the column
end
channel_map
%%

lfp_temp=lfp_temp_l;
size(lfp_temp)

lfp_temp=(lfp_temp-mean(lfp_temp,1))./std(lfp_temp,[],1);
duration=size(lfp_temp,1);
num_channels=size(lfp_temp,2);


%% FILTER AND APPLY HILBERT TRANSFORM

fs=1000;
order=4;
%% --- Choose LFP band center from LFP only ---
domFreq_LFP =Compute_dominant_frequency(lfp_temp');
freqrange = [max(0.5, domFreq_LFP - width_filter), min(fs/2 - 1e-6, domFreq_LFP + width_filter)];
fprintf('LFP-only band centered at %.2f Hz.\n', domFreq_LFP);

%%
LFP_phases=zeros(L_array_y,L_array_x,duration);
Array_x=zeros(L_array_y,L_array_x);
Array_shuffled=zeros(L_array_y,L_array_x);
LFP_data=zeros(L_array_y,L_array_x,duration); 

%%

lfp_temp(isnan(lfp_temp)==1)=0;
% Initialize array to store frequencies
 
%%
 %full-array pseudoinverse for plane fit 
dx_cm = dx;  % already in cm
[Nr, Nc] = deal(L_array_y, L_array_x);
[Xc, Yr] = meshgrid((0:Nc-1)*dx_cm, (0:Nr-1)*dx_cm);
D_full = [Xc(:), Yr(:), ones(Nr*Nc,1)];       % (Nr*Nc) design matrix
pinvD_full = pinv(D_full);  

Simulation_surr.Z_pool = [];
Simulation_surr.resid_pool = [];
%%
for trial=1:number_trials

     initial_time=1;
     final_time=duration;
    
    analytic_signal = zeros(L_array_y, L_array_x, duration);
    for u = 1:num_channels
        LFP_filt = buttfilt(lfp_temp(:,u), freqrange, fs, 'bandpass', order);    
        analytic = hilbert(LFP_filt);   % Hilbert    
        [r,c] = find(channel_map==u);
        analytic_signal(r,c,:) = analytic;
    end
    LFP_phases = angle(analytic_signal);
    LFP_data   = real(analytic_signal);

   % LFP_data: [Ny x Nx x T]
    Z = reshape(permute(analytic_signal,[3 1 2]), duration, []); % [T x Nch]
    m = mean(Z, 2);
    R0_t = abs(m).^2 ./ mean(abs(Z).^2, 2);
    R0   = mean(R0_t,'omitnan');

fprintf('R0 (k=0 power fraction) = %.3f\n', R0);
    
 
    e_resid  = nan(duration,1);
    speed    = nan(duration,1);
    Z_vals  = nan(duration,1);      % Rayleigh statistic per frame
    R_vals  = nan(duration,1);      % (optional) mean resultant length per frame
    % Local (per-trial) surrogate pools
    Z_S_pool_cell = cell(duration,1);
    e_resid_S_pool_cell = cell(duration,1);
    min_grad = 1e-3;                % rad/cm guard

    prev_Array_x = [];

     for tt=1: duration-1
        % Unwrap spatial phase (Code A helper)
        ph_t  = LFP_phases(:,:,tt);
        lfp_t = LFP_data(:,:,tt);
        Array_x = unwrap_phase(ph_t);
        LFP_data(:,:,tt)   = lfp_t;
        LFP_phases(:,:,tt) = Array_x - Array_x(ceil(Nr/2),ceil(Nc/2));


        % ---------- Code B detection pieces start ----------
        % (1) Spatial gradients WITH spacing -> rad/cm
        % Rayleigh Z from unit gradient directions
        [Z_val, r_val, Nvalid, mean_grad_mag] = rayleigh_Z_from_phase(Array_x, dx, min_grad);
        Z_vals(tt)   = Z_val;     % store Rayleigh statistic
        R_vals(tt)   = r_val;     % optional: store mean resultant length
        PGD_D(tt)    = mean_grad_mag;   % keep for speed calc only

        % (2) Plane fit residual (RMS of phase error, in radians)
        phi_vec      = Array_x(:);
        abc_full     = pinvD_full * phi_vec;                 % [a; b; c]
        phi_plane    = reshape(D_full * abc_full, Nr, Nc);
        e_resid(tt)  = sqrt(mean((Array_x(:) - phi_plane(:)).^2));

 
% (3) Conservative pooled surrogates (phase + spatial-block + row/col-shift)
% For thresholds we only need the pooled distributions; we don't need
% per-time storage. Collect to local buffers, then append after the trial.


for km = 1:numel(use_modes)
    mode_k = use_modes{km};
    for kk = 1:K_surr_per_frame
        switch mode_k
            case 'phase' % your original power-spectrum matched surrogate
                Array_surr = phase_scramble_like(Array_x);

            case 'spatialblock' % conservative: shuffle 2x2 (or 3x3) tiles
                Array_surr = spatial_block_shuffle(Array_x, block_sz);

            case 'rowcolshift' % conservative: random cyclic shifts of rows & cols
                Array_surr = rowcol_circshift(Array_x);

            otherwise
                Array_surr = phase_scramble_like(Array_x);
        end

        [Z_s_loc, r_s_loc, ~, ~] = rayleigh_Z_from_phase(Array_surr, dx, min_grad);
        Z_S_pool_cell{tt} = Z_s_loc;
        phi_vec_s   = Array_surr(:);
        abc_full_s  = pinvD_full * phi_vec_s;
        phi_plane_s = reshape(D_full * abc_full_s, Nr, Nc);
        e_resid_loc = sqrt(mean((Array_surr(:) - phi_plane_s(:)).^2));
        e_resid_S_pool_cell{tt} = e_resid_loc;
    end
end


        % Temporal phase derivative for speed (radians/ms -> rad/s after *fs)
        if tt > 1
            % robust per-sample phase diff in (-pi,pi]
            dphi   = angle(exp(1i*Array_x) .* conj(exp(1i*prev_Array_x)));
            dphidt = dphi * fs;       % rad/s
        else
            dphidt = zeros(size(Array_x));
        end
        prev_Array_x = Array_x;

        % Speed estimate; units cm/s because denominator PGD_D is rad/cm and numerator rad/s
        if PGD_D(tt) > min_grad
            speed(tt) = abs(nanmean(dphidt(:))) / PGD_D(tt);   % (rad/s) / (rad/cm) = cm/s
        else
            speed(tt) = NaN;
        end
    end


    % Store trial results (augment Code A with detection summaries)
    Simulation_data(trial).speed   = speed;
    Simulation_data(trial).LFP     = LFP_data;
    Simulation_data(trial).Ph      = LFP_phases;
    Simulation_data(trial).Z       = Z_vals;      % Rayleigh
    Simulation_data(trial).e_resid = e_resid;

    Z_S_pool = vertcat(Z_S_pool_cell{:});
    e_resid_S_pool = vertcat(e_resid_S_pool_cell{:});
    Simulation_surr.Z_pool     = [Simulation_surr.Z_pool;     Z_S_pool];
    Simulation_surr.resid_pool = [Simulation_surr.resid_pool; e_resid_S_pool];

     
end

%% --- Thresholds from pooled conservative surrogates ---
all_Z     = [];  all_resid = [];
for tt = 1:number_trials
    all_Z     = [all_Z;     Simulation_data(tt).Z];
    all_resid = [all_resid; Simulation_data(tt).e_resid];
end

all_Z_shuffled   = Simulation_surr.Z_pool;
all_e_resid_sh   = Simulation_surr.resid_pool;
 
    % High Z threshold from shuffled null, low residual threshold from shuffled null
    threshold_Z = prctile(all_Z_shuffled, threshold_input);   %  
    threshold_e = prctile(all_e_resid_sh,  th_e);            %  
 

%% Quick QA PLOTS
 
figure; hold on
histogram(all_Z,'Normalization','probability')
histogram(all_Z_shuffled,'Normalization','probability')
xline(threshold_Z,'k','LineWidth',2)
title('Rayleigh Z vs surrogate'); legend('Z','Z surrogate','threshold'); grid on


figure; hold on
histogram(all_resid,'Normalization','probability')
histogram(all_e_resid_sh,'Normalization','probability')
xline(threshold_e,'k','LineWidth',2)
title('Plane-fit residual vs surrogate'); legend('resid','resid surrogate','threshold'); grid on


%% --- Detect waves: PGD>thr AND residual<thr (Code B criterion) ---
starts = {}; stops = {}; duration_waves = {};
min_len_ms = 10;                  % your choice
min_len_samples = round(fs * min_len_ms / 1000);

for trial = 1:number_trials
    Z_vals     = Simulation_data(trial).Z(:);
    resid_vals = Simulation_data(trial).e_resid(:);
    wave_detected = (Z_vals > threshold_Z) & (resid_vals < threshold_e);


    starts{trial} = strfind([0 wave_detected'], [0 1]);
    stops{trial}  = strfind([wave_detected' 0], [1 0]);

    % keep only waves lasting at least 10 ms
    lens = stops{trial} - starts{trial} + 1;
    keep = lens >= min_len_samples;
    starts{trial} = starts{trial}(keep);
    stops{trial}  = stops{trial}(keep);
    duration_waves{trial} = stops{trial} - starts{trial} + 1;
end
dh = duration_waves{1};%access for trial 1sdh = duration_waves{1};
%% --- PLOTS PGD & residual with thresholds and wave starts (like A+B) ---

trial = 1;
fg1=figure; pos = get(fg1,'position'); set(fg1,'position',[pos(1) pos(2) 1200 420]);
 
subplot(2,1,1)
plot(Simulation_data(trial).Z(:),'Color',[.8 .6 .6],'LineWidth',1); hold on
plot([1 duration],[threshold_Z threshold_Z],'k')
ylabel('Rayleigh Z'); xlim([0 duration])
subplot(2,1,2)
plot(Simulation_data(trial).e_resid(:),'b'); hold on
plot([1 duration],[threshold_e threshold_e],'k')
ylabel('Residual'); xlabel('Time [ms]'); xlim([0 duration])

 
figure; hold on
for ll=1:numel(starts{trial})
    plot([starts{trial}(ll) starts{trial}(ll)],[-1.5 1.5],'k','LineWidth',2)
end
for jj=1:5
    for ii=1:5
        plot(squeeze(Simulation_data(trial).LFP(ii,jj,:)),'Color',[.6 .6 .8],'LineWidth',1);
        hold on
    end
end
ylabel('LFP'); xlabel('Time [ms]'); xlim([0 duration])

 
%%  Wave velocities / directions / wavelengths  
wave_selected_speed=[];
wave_velocities   = [];
wave_directions   = [];
wave_wavelengths  = [];
 

for trial = 1:number_trials
    for wave = 1:length(starts{trial})

        wave_start = starts{trial}(wave);
        wave_stop  = stops{trial}(wave);

        % average speed for this wave (as in A)  
        wave_speeds      = Simulation_data(trial).speed(wave_start:wave_stop);
        avg_wave_speed   = mean(wave_speeds,'omitnan');
        wave_velocities  = [wave_velocities; avg_wave_speed];
        wave_selected_speed = [wave_selected_speed; wave_speeds(:)];
 

        %  direction (use gradient with dx so units match) 
        phase_event   = Simulation_data(trial).Ph(:,:,wave_start:wave_stop);
        n_time_event  = size(phase_event,3);
        [ny,nx,~] = size(phase_event);


        % Preallocate arrays to hold the spatial gradients at each time point
        FX_event = zeros(L_array_x, L_array_x, n_time_event);
        FY_event = zeros(L_array_x, L_array_x, n_time_event);
        
        % Loop over the time frames of the wave event and compute spatial gradients.

        % --- Plane-fit direction and wavelength ---
        kx_list = nan(n_time_event,1);
        ky_list = nan(n_time_event,1);
        
        for tt = 1:n_time_event
            phi_t = phase_event(:,:,tt);
            abc = pinvD_full * phi_t(:);  % [a; b; c] = [dphi/dx; dphi/dy; offset]
            kx_list(tt) = abc(1);
            ky_list(tt) = abc(2);
        end
        
        % Average across frames
        kx_mean = mean(kx_list, 'omitnan');
        ky_mean = mean(ky_list, 'omitnan');
        
        % Wave-vector direction
        theta_oriented = atan2(ky_list, kx_list);
        wave_directions = [wave_directions; theta_oriented];
        
        % Wavelength from plane-fit slope magnitude
        k_mag = hypot(kx_mean, ky_mean);
        if k_mag >= 1e-3
            lambda_plane = 2*pi / k_mag;  % in cm
        else
            lambda_plane = NaN;
        end
        wave_wavelengths = [wave_wavelengths; lambda_plane];
    end
end


% read entire file into one big string
txt = fileread(fullfile('..', 'run_const_external_rate.sh'));

% look for the pattern "nu_baseline = <number>"
tok = regexp(txt, 'external_input\s*=\s*([\d\.]+)', 'tokens');

% convert first token to number
if ~isempty(tok)
    nu_baseline = str2double(tok{1}{1});
    fprintf('external_input = %g\n', nu_baseline)
else
    error('Couldn''t find external_input in the file.');
end
% read entire file into one big string
 
title_file=sprintf('wavelengths_nu%.3f.mat', nu_baseline);
save(title_file, 'wave_wavelengths');

title_file=sprintf('duration_nu%.3f.mat', nu_baseline);
save(title_file,'dh');
 
figure
polarhistogram(wave_directions, 16)
title('Sim')
set(gca,'fontsize',15)


 
%% Plot waves
all_waves = [];  % each row: [trial, start, stop]
for trial = 1:number_trials
    for w = 1:length(starts{trial})
        all_waves = [all_waves; trial, starts{trial}(w), stops{trial}(w)];
    end
end
Ntot = size(all_waves,1);
Nplot = min(max_nr_waves_to_plot, Ntot);
idx = randperm(Ntot, Nplot);
selected_waves = all_waves(idx,:);
for kk = 1:size(selected_waves,1)     
            trial = selected_waves(kk,1);
            wave_start = selected_waves(kk,2);
            wave_stop  = selected_waves(kk,3);

        LFP_waves_sign = Simulation_data(trial).LFP(:,:,wave_start:wave_stop);

        % Consistent color scaling across frames
        minv = min(LFP_waves_sign,[],'all');
        maxv = max(LFP_waves_sign,[],'all');

        Twave = size(LFP_waves_sign,3);
        frames_per_fig = 50; 
        ncols = 10; 
        nrows = ceil(frames_per_fig / ncols);
        for start_idx = 1:frames_per_fig:Twave
            end_idx = min(start_idx + frames_per_fig - 1, Twave);
            figure('Name', sprintf('Wave %d', start_idx));
            for jj = start_idx:end_idx
                subplot_idx = jj - start_idx + 1;
                subplot(nrows, ncols, subplot_idx);
                imagesc(LFP_waves_sign(:,:,jj), 'Interpolation', 'bilinear');
                set(gca, 'YDir', 'normal');
                title(sprintf('t = %d ms', wave_start + jj - 1), 'FontSize', 8);
                xticklabels([]); yticklabels([]);
                caxis([minv, maxv]);
            end
            colormap(brewermap(256,'*RdBu'));
            han = axes('Position',[0 0 1 1],'Visible','off');
            cb = colorbar('Position',[0.92 0.11 0.02 0.77]); % 
        end
    end

    %% --- Helper: phase-scramble surrogate (from Code B) ---
    function Xs = phase_scramble_like(X)
        %PHASE_SCRAMBLE_LIKE  2-D surrogate with the same power spectrum as X.
        % Keeps |FFT2(X)|, replaces phases with those of a real white-noise field
        % (guarantees Hermitian symmetry), then IFFT -> real map.
        [nr, nc] = size(X);
        F  = fft2(X);
        A  = abs(F);
        R     = randn(nr, nc);
        Theta = angle(fft2(R));
        % Keep DC phase consistent
        Theta(1,1) = angle(F(1,1));
        Fs = A .* exp(1i*Theta);
        Xs = real(ifft2(Fs));
        % Match mean exactly
        Xs = Xs - mean(Xs(:)) + mean(X(:));
    end

    function Xs = spatial_block_shuffle(X, bsz)
    % Shuffle non-overlapping bsz×bsz tiles across the grid.
    % Preserves local within-block gradients; breaks global planes.
    [nr, nc] = size(X);
    assert(mod(nr,bsz)==0 && mod(nc,bsz)==0, 'Grid must be divisible by block size');
    Br = nr/bsz; Bc = nc/bsz;
    
    % Extract tiles
    tiles = cell(Br,Bc);
    for rr = 1:Br
        for cc = 1:Bc
            r0 = (rr-1)*bsz + 1;
            c0 = (cc-1)*bsz + 1;
            tiles{rr,cc} = X(r0:r0+bsz-1, c0:c0+bsz-1);
        end
    end
    
    % Permute tile positions
    perm_idx = randperm(Br*Bc);
    Xs = zeros(size(X));
    k = 1;
    for rr = 1:Br
        for cc = 1:Bc
            [pr, pc] = ind2sub([Br Bc], perm_idx(k));
            r0 = (rr-1)*bsz + 1; c0 = (cc-1)*bsz + 1;
            Xs(r0:r0+bsz-1, c0:c0+bsz-1) = tiles{pr,pc};
            k = k + 1;
        end
    end
    end
    
    function Xs = rowcol_circshift(X)
    % Random cyclic shifts of every row and every column (independently).
    % Preserves row/column marginals & local autocorr; breaks global alignment.
    [nr, nc] = size(X);
    Xs = X;
    
    % Random shift per row
    shifts_r = randi([0 nc-1], nr, 1);
    for r = 1:nr
        Xs(r,:) = circshift(Xs(r,:), shifts_r(r), 2);
    end
    
    % Random shift per column
    shifts_c = randi([0 nr-1], nc, 1);
    for c = 1:nc
        Xs(:,c) = circshift(Xs(:,c), shifts_c(c), 1);
    end
    end

 

function [Z, r, Nvalid, mean_grad_mag] = rayleigh_Z_from_phase(Array_x, dx, min_grad)
    % Gradients (rad/cm)
    [FX, FY] = gradient(Array_x, dx);
    mags = hypot(FX(:), FY(:));              % |grad| (rad/cm)
    ok = isfinite(mags) & (mags >= min_grad);
    Nvalid = sum(ok);
    if Nvalid == 0
        Z = 0; r = 0; mean_grad_mag = NaN; 
        return
    end
    % Unit gradient direction vectors u_i
    ux = FX(:)./mags; uy = FY(:)./mags;
    U  = [ux(ok), uy(ok)];                   % [Nvalid x 2]
    mU = mean(U,1);                          % mean direction vector
    r  = sqrt(sum(mU.^2));                   % mean resultant length in 2D
    Z  = Nvalid * r.^2;                      % Rayleigh statistic
    % Keep mean |grad| for your speed calc (cm/s = (rad/s)/(rad/cm))
    mean_grad_mag = mean(mags(ok));
end
