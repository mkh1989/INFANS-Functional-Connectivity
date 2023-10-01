function [coh, icoh, lcoh] = infans_compute_coh_icoh_lcoh(data, fs)
% =========================================================================
%
% This function calculates coherence, imaginary part of coherence, and lagged coherence.
%
% INPUTS:
%   data    : a matrix containing EEG recordings (samples * channels)
%   fs      : the sampling frequency of the data
%
% OUTPUT:
%   coh     : a matrix containing the max of magnitude squared coherence among all frequencies (channels * channels) 
%   icoh    : a matrix containing the max of imaginary part of coherence among all frequencies (channels * channels) 
%   lcoh    : a matrix containing the max of lagged coherence among all frequencies (channels * channels)
%     
% NOTE:
%   cxy:  cross-spectral density between x and y
%   cxx:  autospectral density of x
%   cyy:  autospectral density of y
%   Coherence function (c)            : cxy/sqrt(cxx*cyy)
%   Magnitude-squared Coherence (coh) : |c|^2 = |cxy|^2/(cxx*cyy) 
%   Imaginary Coherence (icoh)        : abs(imag(c))               
%   Lagged Coherence (lcoh)           : abs(imag(c))/sqrt(1-real(c)^2)
%   
% =========================================================================

% gets the size of the data
[samples, channels] = size(data);

% setting parameters
olength = fix(samples / 10);
wlength = 2 * olength;
win     = hamming(wlength);
nfft    = max(512, pow2(nextpow2(wlength)));

% allocates memory for outputs
coh  = zeros(channels);
icoh = zeros(channels);
lcoh = zeros(channels);

% calculates the indexes for each pair of sensors
for ch1 = 1: channels - 1
    for ch2 = ch1 + 1: channels
        [pxy,~]        = cpsd(data(: , ch1), data(: , ch2), win, olength, nfft, fs);
        [pxx,~]        = cpsd(data(: , ch1), data(: , ch1), win, olength, nfft, fs);
        [pyy,~]        = cpsd(data(: , ch2), data(: , ch2), win, olength, nfft, fs);
        c              = pxy ./ sqrt(pxx .* pyy);
        msc            = max(abs(c) .^ 2);
        imc            = max(abs(imag(c)));
        lagc           = max(abs(imag(c)) ./ sqrt(1-real(c).^2));
        coh (ch1, ch2) = msc;
        coh (ch2, ch1) = msc;
        icoh(ch1, ch2) = imc;
        icoh(ch2, ch1) = imc;
        lcoh(ch1, ch2) = lagc;
        lcoh(ch2, ch1) = lagc;
    end
end

end


