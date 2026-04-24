
% =========================================================================
%
% Detects travelling waves in macaque cortical LFP data using spatial phase
% gradient analysis. Waves are identified by high Rayleigh Z (coherent
% gradient direction) and low plane-fit residual (planar wavefront).
% Surrogate thresholds are estimated from spatially shuffled null data.
% 
% =========================================================================


function wave_detection_Macaque_2025_test_R(...
     chanlfp,beh,number_trials,duration,width_filter,threshold_input,th_e,max_nr_waves_to_plot,...
    use_modes, K_surr_per_frame, block_sz)

 
if nargin < 9 || isempty(use_modes),       use_modes = {'spatialblock'}; end
if nargin < 10 || isempty(K_surr_per_frame),K_surr_per_frame = 2; end
if nargin < 11 || isempty(block_sz),        block_sz = 2; end


rng(1); % reproducible
 
fields_lfp   = fieldnames(chanlfp);
units_lfp    = fields_lfp(startsWith(fields_lfp,'lfp'));
num_channels = length(units_lfp);
size_lfp     = length(chanlfp.(units_lfp{1}));
location_channels = chanlfp.MIchan2rc;
location_channels(num_channels+1 :end,:)=[];
tomask    = [2 25 66 67 74]; % Channels excluded prior to analysis (identified as noisy or artifactual)
 

%% --- Parameters 

dx   = 0.04;     % cm (electrode spacing)
dt   = 0.001;    % s
fs   = 1000;     % Hz
n_ch = 10;
order     = 4;          % IIR filter order
min_grad        = 1e-3;       % guard for tiny gradients in speed calc (


%% --- Build position matrix   ---
Matrice_posizione = NaN(10,10);
for u = 1:length(location_channels(:,1))
    Matrice_posizione(location_channels(u,1),location_channels(u,2)) = u;
end

%% --- Gather LFP per channel   ---
all_temp = zeros(size_lfp, num_channels);
for u = 1:num_channels
    all_temp(:,u) = chanlfp.(units_lfp{u});
end
% Remove masked channels from both data and mapping
all_temp(:,tomask) = []; 
num_channels       = size(all_temp,2);
location_channels(tomask,:) = [];
 
%% --- Plane-fit precomputation (for residuals) ---
dx_cm = dx;  % already in cm
[Nr, Nc] = deal(n_ch, n_ch);
[Xc, Yr] = meshgrid((0:Nc-1)*dx_cm, (0:Nr-1)*dx_cm);
D_full = [Xc(:), Yr(:), ones(Nr*Nc,1)];       % (Nr*Nc) x 3 design matrix
pinvD_full = pinv(D_full);                    % precompute pseudoinverse

% Frequency band
ave_LFP=mean(all_temp,2);

%% Let us concatenate 10 trials to estimate the dominant frequency
All_trials = [];
number_trials_f=min(number_trials,10); %10 trial are enough for reliable estimation
for jj = 1:number_trials_f
    initial_time = max(1, floor(beh(jj,1)*fs) - 200);
    final_time   = min(initial_time + duration, size(all_temp,1));
    
    seg = all_temp(initial_time:final_time, :);
    All_trials = [All_trials; seg];
end
domFreq_LFP = Compute_dominant_frequency(All_trials');

freqrange = [max(0.5, domFreq_LFP - width_filter), min(fs/2 - 1e-6, domFreq_LFP + width_filter)];
fprintf('LFP-only band centered at %.2f Hz.\n', domFreq_LFP);
clear All_trials
%% --- Main loop over trials   ---

LFP_phases_nan = NaN(n_ch,n_ch,duration);
LFP_data_nan   = NaN(n_ch,n_ch,duration);
LFP_phases     = NaN(n_ch,n_ch,duration);
LFP_data       = NaN(n_ch,n_ch,duration);


for trial = 1:number_trials

    initial_time = max(1, floor(beh(trial,1)*fs) - 200);
    final_time   = min(initial_time + duration, size(all_temp,1));
    seg_len      = final_time - initial_time + 1;
    if seg_len < 2
        warning('Trial %d skipped: too short segment.', trial);
    continue
    end
 

    % --- Filter & Hilbert per channel, place into 10x10 grid with NaNs where missing  ---
    for u = 1:num_channels
        lfp_temp            = all_temp(:,u);
        seg                 = lfp_temp(initial_time:final_time);
        LFP_data_temp       = buttfilt(seg, freqrange, fs, 'bandpass', order);
        wrapped_phase       = angle(hilbert(LFP_data_temp'));   % (1 x T)
        LFP_phases_nan(location_channels(u,1),location_channels(u,2),1:seg_len-1) = wrapped_phase(1:seg_len-1);
        LFP_data_nan(location_channels(u,1),location_channels(u,2),1:seg_len-1)   = LFP_data_temp(1:seg_len-1);
    end

    T = seg_len - 1;                 % instead of duration-1
    PGD     = zeros(T,1); PGD_D = zeros(T,1);
    PGD_S   = zeros(T,1); PGD_DS = zeros(T,1);
    e_resid = zeros(T,1); 
    speed   = NaN(T,1);
    % Local (per-trial) surrogate pools
    PGD_S_pool    = [];
    e_resid_S_pool = [];
    Z_vals  = zeros(T,1);      % Rayleigh statistic per frame
    R_vals  = zeros(T,1);      % (optional) mean resultant length per frame
    % Local (per-trial) surrogate pools
    Z_S_pool       = [];
    e_resid_S_pool = [];


    clear prev_Array_x
    for tt = 1:T
        % Inpaint to full grid 
        ph_t  = inpaint_nans(squeeze(LFP_phases_nan(:,:,tt)), 3);
        lfp_t = inpaint_nans(squeeze(LFP_data_nan  (:,:,tt)), 3);

        % Unwrap spatial phase  
        Array_x = unwrap_phase(ph_t);     % radians
        LFP_data(:,:,tt)   = lfp_t;       % keep full LFP grid
        LFP_phases(:,:,tt) = Array_x - Array_x(3,3);  % only for later plotting consistency

        % ----------   detection pieces start ----------
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
for km = 1:numel(use_modes)
    mode_k = use_modes{km};
    for kk = 1:K_surr_per_frame
        switch mode_k
            case 'phase' 
                Array_surr = phase_scramble_like(Array_x);

            case 'spatialblock' % conservative: shuffle 2x2 (or 3x3) tiles
                Array_surr = spatial_block_shuffle(Array_x, block_sz);

            case 'rowcolshift' % conservative: random cyclic shifts of rows & cols
                Array_surr = rowcol_circshift(Array_x);

            otherwise
                Array_surr = phase_scramble_like(Array_x);
        end


        [Z_s_loc, r_s_loc, ~, ~] = rayleigh_Z_from_phase(Array_surr, dx, min_grad);
        Z_S_pool   = [Z_S_pool; Z_s_loc];

        phi_vec_s   = Array_surr(:);
        abc_full_s  = pinvD_full * phi_vec_s;
        phi_plane_s = reshape(D_full * abc_full_s, Nr, Nc);
        e_resid_loc = sqrt(mean((Array_surr(:) - phi_plane_s(:)).^2));
        e_resid_S_pool = [e_resid_S_pool; e_resid_loc];
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

    % Store trial results  
    Macaque_data(trial).speed   = speed;
    Macaque_data(trial).LFP     = LFP_data;
    Macaque_data(trial).Ph      = LFP_phases;

    Macaque_data(trial).Z       = Z_vals;      % Rayleigh
    Macaque_data(trial).e_resid = e_resid;
 


    if ~exist('Macaque_surr','var')
        Macaque_surr.Z_pool     = [];
        Macaque_surr.resid_pool = [];
    end
    Macaque_surr.Z_pool     = [Macaque_surr.Z_pool;     Z_S_pool];
    Macaque_surr.resid_pool = [Macaque_surr.resid_pool; e_resid_S_pool];

    LFP_phases_nan(:) = NaN;
    LFP_data_nan(:)   = NaN;
    LFP_phases(:)     = NaN;
    LFP_data(:)       = NaN;

end  %% end loop over trials


%% --- Thresholds from pooled conservative surrogates ---

all_Z_shuffled   = Macaque_surr.Z_pool;
all_e_resid_sh   = Macaque_surr.resid_pool;


    % High Z threshold from shuffled null, low residual threshold from shuffled null
threshold_Z = prctile(all_Z_shuffled, threshold_input);   % e.g., 95–99.9
threshold_e = prctile(all_e_resid_sh,  th_e);            % e.g., 5–20
 

%% --- Detect waves: PGD>thr AND residual<thr   ---
starts = {}; stops = {}; duration_waves = {};
min_len_ms = 10;                  % your choice
min_len_samples = round(fs * min_len_ms / 1000);

for trial = 1:number_trials
    Z_vals     = Macaque_data(trial).Z(:);
    resid_vals = Macaque_data(trial).e_resid(:);
    wave_detected = (Z_vals > threshold_Z) & (resid_vals < threshold_e);


    starts{trial} = strfind([0 wave_detected'], [0 1]);
    stops{trial}  = strfind([wave_detected' 0], [1 0]);

    % keep only waves lasting at least 10 ms
    lens = stops{trial} - starts{trial} + 1;
    keep = lens >= min_len_samples;
    starts{trial} = starts{trial}(keep);
    stops{trial}  = stops{trial}(keep);
end

%% --- PLOTS PGD & residual with thresholds and wave starts  ---

trial = 1;
fg1=figure; pos = get(fg1,'position'); set(fg1,'position',[pos(1) pos(2) 1200 420]);
subplot(2,1,1)
plot(Macaque_data(trial).Z(:),'Color',[.8 .6 .6],'LineWidth',1); hold on
plot([1 duration],[threshold_Z threshold_Z],'k')
ylabel('Rayleigh Z'); xlim([0 duration])

subplot(2,1,2)
plot(Macaque_data(trial).e_resid(:),'b'); hold on
plot([1 duration],[threshold_e threshold_e],'k')
ylabel('Residual'); xlabel('Time [ms]'); xlim([0 duration])

figure; hold on
for ll=1:numel(starts{trial})
    plot([starts{trial}(ll) starts{trial}(ll)],[-150 150],'k','LineWidth',2)
end
for jj=1:5
    for ii=1:5
        plot(squeeze(Macaque_data(trial).LFP(ii,jj,:)),'Color',[.6 .6 .8],'LineWidth',1);
        hold on
    end
end
ylabel('LFP'); xlabel('Time [ms]'); xlim([0 duration])


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
            LFP_waves_sign = Macaque_data(trial).LFP(:,:,wave_start:wave_stop);  
             % Consistent color scaling across frames
            minv = min(LFP_waves_sign,[],'all');
            maxv = max(LFP_waves_sign,[],'all');
            % How many frames per row
            Twave = size(LFP_waves_sign,3);
            ncols = 10;
            nrows = ceil(Twave / ncols);      
            figure('Name',sprintf('Wave %d (trial %d): LFP frames',kk,trial));
            for jj = 1:Twave
                subplot(nrows, ncols, jj);
                imagesc(LFP_waves_sign(:,:,jj),'Interpolation', 'bilinear')
                set(gca,'YDir','normal');
                title(sprintf('t = %d ms', wave_start + jj - 1), 'FontSize', 8);
                xticklabels([]); yticklabels([]);
                caxis([minv, maxv]);
            end        
            % Colormap: try brewermap if available, otherwise fall back
            if exist('brewermap','file')
                colormap(brewermap(256,'*RdBu'));
            else
                colormap(parula);
            end
            han = axes('Position',[0 0 1 1],'Visible','off');
            cb = colorbar('Position',[0.92 0.11 0.02 0.77]); %#ok<NASGU>
             
        end %end loop on selected waves
%% --- Wave velocities / directions / wavelengths  ---
wave_selected_speed=[];
wave_velocities   = [];
wave_directions   = [];
wave_wavelengths  = [];
wave_wavelengths2 = [];

for trial = 1:number_trials
    for wave = 1:length(starts{trial})
        wave_start = starts{trial}(wave);
        wave_stop  = stops{trial}(wave);

        % ---- average speed for this wave ----
        wave_speeds      = Macaque_data(trial).speed(wave_start:wave_stop);
        avg_wave_speed   = mean(wave_speeds,'omitnan');
        wave_velocities  = [wave_velocities; avg_wave_speed];
        wave_selected_speed = [wave_selected_speed; wave_speeds(:)];
    

        % ---- direction (use gradient with dx so units match) ----
        phase_event   = Macaque_data(trial).Ph(:,:,wave_start:wave_stop);
        n_time_event  = size(phase_event,3);
        [ny,nx,~] = size(phase_event);
         % Preallocate arrays to hold the spatial gradients at each time point
        FX_event = zeros(n_ch, n_ch, n_time_event);
        FY_event = zeros(n_ch, n_ch, n_time_event);
        
        % Loop over the time frames of the wave event and compute spatial gradients.
        % (We include the physical spacing dx in the gradient call. 
        % For computing the angle, the scale factor cancels.)
        for tt = 1:n_time_event
            [fx, fy] = gradient(phase_event(:,:,tt), dx);
            FX_event(:,:,tt) = fx;
            FY_event(:,:,tt) = fy;
        end
        
        % Mean temporal phase slope (rad/s)
        dphidt_wave = Macaque_data(trial).Ph(:,:,wave_start:wave_stop);
        dphidt_wave = angle(exp(1i*dphidt_wave(:,:,2:end)) .* ...
                            conj(exp(1i*dphidt_wave(:,:,1:end-1))));
        dphidt_mean = fs * mean(dphidt_wave(:),'omitnan');   %  
        
        % Oriented direction:
        % Average spatial gradients (rad/cm)
        % Attention to axis conventions in matlab!!!
        FX_avg = mean(FX_event,3);   % FX = d phi / dy (rows)
        FY_avg = mean(FY_event,3);   % FY = d phi/dx (cols)
        
        % Map to physical axes: x = columns, y = rows
        kx = mean(FY_avg(:));        % d phi/dx
        ky = mean(FX_avg(:));        % d phi/dy
        
        % Unsigned gradient direction (k direction)
        theta_grad = atan2(ky, kx);

%  
        if dphidt_mean > 0
            % 
            theta_oriented = wrapToPi(theta_grad + pi);
        else
            theta_oriented = theta_grad;
        end
        wave_directions = [wave_directions; theta_oriented];

        % Wavelength from gradient
        avg_grad_norm = hypot(kx, ky);   % rad/cm
        if avg_grad_norm < 1e-3
            % too flat to estimate a reliable wavelength
            lambda_gradient = NaN;
        else
            lambda_gradient = 2*pi / avg_grad_norm;  % cm
        end
        wave_wavelengths = [wave_wavelengths; lambda_gradient];

        % ---- optional: wavelength from speed/frequency (quick zero-crossing est.) ----
        phase_signal = reshape(phase_event, n_ch*n_ch, n_time_event);
        mean_phase_signal = mean(phase_signal,1);
        zx = find(diff(sign(mean_phase_signal))~=0);
        if ~isempty(zx)
            estimated_period = 2*mean(diff(zx));  % samples
            nu_wave = fs / estimated_period;      % Hz
        else
            nu_wave = NaN;
        end
        lambda_speed_freq = avg_wave_speed / nu_wave;     % cm
        wave_wavelengths2 = [wave_wavelengths2; lambda_speed_freq];        
    end
end

% Pack results (same names as A)
Macaque_final.wave_directions = wave_directions;
Macaque_final.wave_wavelengths = wave_wavelengths;
Macaque_final.wave_wavelengths2 = wave_wavelengths2;
Macaque_final.all_speeds = wave_selected_speed;

titlefile=sprintf('Macaque_R_trials_%d_DWfilter_%d_thPGD_%.2f_the_%d.mat',number_trials,width_filter,threshold_input,th_e);
save(titlefile,'Macaque_data','Macaque_final');


    %% --- Helper: phase-scramble surrogate   ---
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
