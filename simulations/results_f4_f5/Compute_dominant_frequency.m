function f0_ave_LFP=Compute_dominant_frequency(prep_LFP)

n_channel=size(prep_LFP,1);

% prep_LFP is [nChan x T] 
fs      = 1000;
nperseg = 1000; 
noverlap = nperseg/2;
fmin = 5; fmax = 150;

nChan = size(prep_LFP,1);
f0_ch = nan(nChan,1);              % per-channel peak (mean over segments)
sigma_within = nan(nChan,1);       % optional: within-channel segment SD

 
 % Channel-averaged PSD, then one peak:
nperseg = 2000; noverlap = nperseg/2; fs = 1000; fmin=5; fmax=150;
[nChan, ~] = size(prep_LFP);
% compute PSD per channel once
for c=1:nChan
    sig = prep_LFP(c,:);
    [Pxx(:,c), f] = pwelch(sig-mean(sig), hamming(nperseg), noverlap, [], fs);
end
band = f>=fmin & f<=fmax;

% FREQUENCY FROM median: parabolic correction around max
%Pmean = mean(Pxx(band,:),2);  fb = f(band);
Pmean = median(Pxx(band,:),2);  fb = f(band);
[~,k] = max(Pmean); dF = fb(2)-fb(1);
if k>1 && k<numel(Pmean)
    a=Pmean(k-1); b=Pmean(k); g=Pmean(k+1);
    p = 0.5*(a-g)/(a-2*b+g);
    f0_LFP = fb(k) + p*dF;
else
    f0_LFP = fb(k);
end


f0_LFP_sing=nan(size(Pxx,2),1);
for j=1:size(Pxx,2)    
    sing_pxx=Pxx(band,j);
    [~,k] = max(sing_pxx);  
    if k>1 && k<numel(sing_pxx)
        a=sing_pxx(k-1); b=sing_pxx(k); g=sing_pxx(k+1);
        p = 0.5*(a-g)/(a-2*b+g);
        f0_LFP_sing(j) = fb(k) + p*dF;
    else
        f0_LFP_sing(j) = fb(k);
    end
end

f0_ave_LFP=mean(f0_LFP_sing);
sd_LFP=std(f0_LFP_sing);

end



function f0_segments = jackknife_peak_freq(sig, fs, nperseg, noverlap, fmin, fmax)
    step = nperseg - noverlap;
    idx = 1:step:(numel(sig)-nperseg+1);
    f0_segments = nan(1,numel(idx));
    for ii=1:numel(idx)
        seg = sig(idx(ii):idx(ii)+nperseg-1);
        [f0, ~] = peak_features_from_psd(seg, fs, fmin, fmax);
        f0_segments(ii) = f0;
    end
end

function [f0, fwhm, f_centroid, sig_centroid] = peak_features_from_psd(sig, fs, fmin, fmax)
    nperseg = 1000; noverlap = nperseg/2;
    [PSD,f] = pwelch(sig - mean(sig), hamming(nperseg), noverlap, [], fs);

    band = (f>=fmin) & (f<=fmax);
    fb = f(band); Pb = PSD(band);

    % --- centroid (for reference only; often too large) ---
    Ptot = sum(Pb);
    f_centroid = sum(fb.*Pb) / Ptot;
    sig_centroid = sqrt(sum(Pb.*(fb - f_centroid).^2)/Ptot);

    % --- peak index and parabolic interpolation of peak frequency ---
    [Ppk, k] = max(Pb);
    dF = fb(2) - fb(1);
    if k>1 && k<numel(Pb)
        alpha=Pb(k-1); beta=Pb(k); gamma=Pb(k+1);
        p = 0.5*(alpha - gamma)/(alpha - 2*beta + gamma);  % vertex offset in bins
        f0 = fb(k) + p*dF;
        P0 = beta - 0.25*(alpha - gamma)*p;                 % refined peak power (optional)
    else
        f0 = fb(k);
        P0 = Ppk;
    end

    % --- FWHM (–3 dB width) with linear interpolation on each side ---
    halfP = P0*0.5;

    % left side
    iL = k; 
    while iL>1 && Pb(iL)>halfP, iL=iL-1; end
    if iL==1
        fL = fb(1);
    else
        % interpolate between (iL, iL+1)
        t = (halfP - Pb(iL)) / (Pb(iL+1)-Pb(iL));
        fL = fb(iL) + t*(fb(iL+1)-fb(iL));
    end

    % right side
    iR = k;
    while iR<numel(Pb) && Pb(iR)>halfP, iR=iR+1; end
    if iR==numel(Pb)
        fR = fb(end);
    else
        % interpolate between (iR-1, iR)
        t = (halfP - Pb(iR-1)) / (Pb(iR)-Pb(iR-1));
        fR = fb(iR-1) + t*(fb(iR)-fb(iR-1));
    end

    fwhm = max(0, fR - fL);
end
