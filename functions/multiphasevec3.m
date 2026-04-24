function [phase,pow] = multiphasevec3(f,S,Fs,width)
% FUNCTION [phase,pow] = multiphasevec3(f,S,Fs,width)
%
% Returns the phase and power as a function of time for a range of
% frequencies (f).
%
% INPUT ARGS:
%   f = [2 4 8];   % Frequencies of interest
%   S = dat;       % Signal to process (samples x chan)
%   Fs = 256;      % Sampling frequency
%   width = 6;     % Width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   phase- Phase data [freqs,time]
%   power- Power data [freqs,time]

nF = length(f);
nS = size(S);
dt = 1/Fs;
st = 1./(2*pi*(f/width));

% Preallocate cell array for wavelets
curWaves = cell(1, nF);

% Generate wavelets for each frequency
for i = 1:nF
    % Create time vector manually instead of using linspace
    time_vec = -3.5*st(i):dt:3.5*st(i);
    curWaves{i} = morlet(f(i), time_vec, width);
end

nCurWaves = cellfun(@length, curWaves);
Lys = nS(2) + nCurWaves - 1;    % length of convolution of S and curWaves{i}

% Compute next power of 2 for each convolution length
Ly2s = arrayfun(@(x) pow2(nextpow2(x)), Lys);

ind1 = ceil(nCurWaves/2);       % start index of signal after convolution
pow = zeros(nS(1), nF, nS(2));
phase = zeros(nS(1), nF, nS(2));

for i = 1:nF
    Ly2 = Ly2s(i);
    
    % Perform convolution using FFT
    y1 = ifft(bsxfun(@times, fft(S,Ly2,2), fft(curWaves{i},Ly2)));
    y1 = y1(:, ind1(i):(ind1(i)+nS(2)-1));
    
    % Compute power and phase
    pow(:,i,:) = abs(y1).^2;
    phase(:,i,:) = angle(y1);
end
end

% Morlet wavelet function (you might want to modify this as well)
function y = morlet(f, t, width)
% Basic Morlet wavelet construction
sigma = width / (2*pi*f);
y = exp(-t.^2 / (2*sigma^2)) .* exp(1i*2*pi*f*t);
y = y ./ sqrt(sum(abs(y).^2));
end