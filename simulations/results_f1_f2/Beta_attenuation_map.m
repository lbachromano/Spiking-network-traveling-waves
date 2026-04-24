% Generates spatial maps of Beta attenuation latencies via Hilbert envelopes.
% Identifies dominant Beta peaks, filters LFP, and computes threshold-crossing latencies.

clear

addpath(fullfile(pwd, '..', '..', 'functions'));


folder_name='examples_paper';   
trial_number=8; %used as example %1 11


name = fullfile(folder_name, sprintf('LFP_absSum_i%d_J1.mat', trial_number));
CR_x = load(fullfile(folder_name,'Rec_x.mat'));  Rec_x = CR_x.data(:);
CR_y = load(fullfile(folder_name,'Rec_y.mat')); Rec_y = CR_y.data(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold     = 0.8;      % envelope crossing threshold
win_ms   = -400 : 600; %window around tref_ms where to look for the beta attenuation time
order= 4;
fs = 1000;     % sampling rate (Hz)
skip_initial_window=50; %skip the first 100 ms containing artifacts
preparation=1000; %ms; the time when the external rate starts rising
tref_ms=1200; %point of mid raise of ext input
width_filter=5; %for filter
beta_range=[13,30];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LFP_from_brian=load(name);
lfp_temp=(LFP_from_brian.data)';
lfp_temp(:,1:skip_initial_window)=[]; 
T=length(lfp_temp); 
n_channel=size(lfp_temp,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Fix 2: Better Segment Handling ---
lfp_preparation = lfp_temp(:, skip_initial_window:preparation);
duration = size(lfp_preparation, 2); 
%Use a smaller nperseg if duration is short to allow more averaging
nperseg  = min(duration, 512); 
noverlap = floor(nperseg/2);
nfft     = 2^nextpow2(max(1024, nperseg));

% 1. Detrend and compute PSD
x_mean = detrend(mean(lfp_preparation, 1, 'omitnan'));
[Pm, f_psd] = pwelch(x_mean, hamming(nperseg), noverlap, nfft, fs);

% 2. Log-transform or Flatten (to handle 1/f)
Pm_log = 10*log10(Pm); 

% 3. Define Beta Search Window
beta_idx = (f_psd >= beta_range(1)) & (f_psd <= beta_range(2));
f_beta = f_psd(beta_idx);
p_beta = Pm_log(beta_idx);

% 4. Find the local peak within that window
[pks, locs] = findpeaks(p_beta, f_beta, 'SortStr', 'descend');

if ~isempty(pks)
    domFreq_LFP = locs(1); % The highest local peak in the Beta band
else
    % Fallback: if no local peak, take the maximum in range
    [~, k] = max(p_beta);
    domFreq_LFP = f_beta(k);
end
freqrange = [max(0.5, domFreq_LFP - width_filter), min(fs/2 - 1e-6, domFreq_LFP + width_filter)];
fprintf('LFP-only band centered at %.2f Hz.\n', domFreq_LFP);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BETA ENVELOPES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

envelopes=zeros(n_channel,T);
for u=1:n_channel

    filtered_lfp=(buttfilt(lfp_temp(u,:),freqrange,fs,'bandpass',order));
    x =filtered_lfp;
    t = (0:1:length(x)-1)./1000;
    envelopes(u,:)=abs(hilbert(x));

end
normalized_full_envelopes = envelopes ./ max(envelopes, [], 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- load recording-site coordinates (meters) ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xu = round(Rec_x, 9);   % meters 
yu = round(Rec_y, 9);

% --- get sorted unique coordinates along each axis ---
xs = sort(unique(xu));               % all distinct x positions
ys = sort(unique(yu));               % all distinct y positions
edge_l = numel(xs);                  % should equal numel(ys) == grid size

 
[~, x_R] = ismember(xu, xs);         % 1..edge_l (column index)
[~, y_R] = ismember(yu, ys);         % 1..edge_l (row index)

% Optional sanity check: channel labels at their (x,y) locations (mm)
figure; hold on
for u = 1:numel(xu)
    text(1000*Rec_x(u), 1000*Rec_y(u), num2str(u), 'HorizontalAlignment','center');
end
xlabel('x (mm)'); ylabel('y (mm)'); axis equal; box on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND BETA ATTENUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

win_samp = tref_ms +  win_ms;  % vector of sample indices

% 3) for each channel, find the local peak?s index within that window
mean_envelope=mean(envelopes,1); %average envelope across channels

[max_values, rel_idx]    =  findpeaks( mean_envelope(win_samp)' );  %find the time of the peak of averaged enveloped

%selected_index_pk_envelope=max(rel_idx); %to pick up the last max
% selected_index_pk_envelope=min(rel_idx); %to pick the first max
% 
[~,highest_peak] =max(max_values);
 selected_index_pk_envelope=rel_idx(highest_peak); %to pick the highest peak

event_sample =win_samp(selected_index_pk_envelope);
 

figure
plot(normalized_full_envelopes')
hold on
plot([event_sample event_sample], [0 1.2],'k','linewidth',2)
xlabel('time')
ylabel('beta envelope')
set(gca,'fontsize',18)




start_ms      = -10;     % window start, in ms relative to event
end_ms        = +200;     % window end,   in ms relative to event

% 2) Convert ms to sample offsets
start_offset = round( start_ms * fs/1000 );   % negative number
end_offset   = round( end_ms   * fs/1000 );   % positive number

% 3) Build your sample window (e.g. 1241:1591 if event_sample=1400)
win_samples  = (event_sample + start_offset) : (event_sample + end_offset);

% 4) Automatically build the time?axis in ms
time_axis_ms = (win_samples ) * (1000/fs);

% 5) Extract & normalize the envelopes in that window
zoom_env = envelopes(:, win_samples);
norm_env = zoom_env ./ max(zoom_env,[],2);

% 6) Find each channel?s crossing index & convert to ms
n_chan = size(envelopes,1);
latency_ms = nan(n_chan,1);
for u=1:n_chan
    [~, idx]      = min( (norm_env(u,:) - threshold).^2 );
    latency_ms(u) = time_axis_ms(idx);
end



%%

 
LFP_latency_map = nan(edge_l, edge_l);
for u = 1:numel(x_R)
    r = y_R(u); c = x_R(u);
    if r>0 && c>0
        LFP_latency_map(r, c) = latency_ms(u);
    end
end

% figure;
% imagesc(LFP_latency_map);
% set(gca,'YDir','normal'); colorbar; axis image
% xticks([]); yticks([]);
% title('\beta-envelope latency (ms)')


 
%%
figure
plot(time_axis_ms,norm_env','linewidth',1.2)
xlabel('time')
ylabel('beta envelope')
set(gca,'fontsize',18)
axis square

% Mask out outliers

% Suppose latency_ms is your [n_chan×1] vector of latencies
outlier_mask = isoutlier(latency_ms,'median');  

LFP_latency_masked = LFP_latency_map;   % edge_l×edge_l
for u = 1:n_channel
    if outlier_mask(u)
        LFP_latency_masked(y_R(u),x_R(u) ) = NaN;
    end
end


figure
h = imagesc(LFP_latency_masked);          % plot the matrix as-is
h.AlphaData = ~isnan(LFP_latency_masked); % keep your transparency mask

set(gca,'YDir','normal')      % or equivalently:  axis xy
                          %   (undoes the automatic flip)
colorbar
%title('Masked \beta Envelope Latencies (outliers hidden)')
xticks([]); yticks([])


% Outliers masked out and substituted with the moving mean of surrounding values
windowSize = [3 3]; 
% 2. Fill NaNs using the moving mean of surrounding values
LFP_filled = fillmissing(LFP_latency_masked, 'movmean', windowSize);
% 3. Plot the smoothed/filled result
figure
h = imagesc(LFP_filled);
set(gca,'YDir','normal')      % or equivalently:  axis xy
                             %   (undoes the automatic flip)
colorbar
%title('Masked \beta Envelope Latencies (outliers hidden)')
xticks([]); yticks([])
