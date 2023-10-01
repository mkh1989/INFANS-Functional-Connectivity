function [pxy, f] = infans_cpsd(y, x, fs)
% =========================================================================
%
% This function calculates cross power spectrum density (cpsd).
%
% INPUTS:
%   x   : the first signal/channel
%   y   : the second signal/channel
%   fs  : sampling frequency
%
% OUTPUT:
%   pxy : a vector contains the values of cross power specteral density
%   f   : a vector contains the frequencies at which cpsd are calculated
%
% NOTE:
%   we use Welch’s averaged, modified periodogram method to estimate the spectrum,
%   as we are dealing with finite data. Both the windowing of the data and
%   the use of Welch’s averaged periodogram reduce the frequency resolution of the coherence.
%   Welch’s periodogram,by default, uses segments of 2/9 times the window length,
%   thus reducing the frequency resolution to approximately a fifth of its value
%
% =========================================================================

% set parameters
samples = length(x);
step    = floor(samples / 9);
wlen    = step * 2;
win     = hamming(wlen);
nfft    = max(256, pow2(nextpow2(wlen)));

% claculate the number of segments
segments = floor((samples - wlen) / step) + 1;

% segment the data
win        = win / norm(win) / sqrt(pi);
windowed_x = zeros(wlen, segments);
windowed_y = zeros(wlen, segments);

for s = 1: segments
    windowed_x(:, s) = win .* x((s - 1) * step + (1: wlen));
    windowed_y(:, s) = win .* y((s - 1) * step + (1: wlen));
end

% get Fourier transform of the segmented data
FTx             = fft(windowed_x, nfft, 1);
FTx             = FTx(1: ceil((nfft + 1) / 2), :);
FTx([1 end], :) = FTx([1 end], :) / sqrt(2);

FTy             = fft(windowed_y, nfft, 1);
FTy             = FTy(1: ceil((nfft + 1) / 2), :);
FTy([1 end], :) = FTy([1 end], :) / sqrt(2);

% calculate the cpsd
pxy = mean(conj(FTx) .* FTy, 2);

% return the frequencies
f = linspace(0, pi, nfft) * (fs/2);
end