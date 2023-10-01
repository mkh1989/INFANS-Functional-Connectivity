function [plv, pli, wpli] = infans_compute_plv_pli_wpli(data, fs)
% =========================================================================
%
% This function calculates Phase Synchronization Indexes.
%
% INPUTS:
%   data : a matrix containing EEG recordings (samples * channels)
%   fs      : the sampling frequency of the data (just for wpli)

% OUTPUT:
%   plv  : a matrix containing phase locking values between all possible pairs of channels (channels * channels)
%   pli  : a matrix containing phase lag indexes between all possible pairs of channels (channels * channels)
%   wpli : a matrix containing wighted phase lag indexes between all possible pairs of channels (channels * channels)
% 
% =========================================================================

% gets the size of the data
[samples, channels] = size(data);

% setting parameters for wPLI
olength = fix(samples / 10);
wlength = 2 * olength;
win     = hamming(wlength);
nfft    = max(512, pow2(nextpow2(wlength)));

% allocates memory for outputs
plv  = zeros(channels);
pli  = zeros(channels);
wpli = zeros(channels);

% computes analytic signal of each channel
analytics = hilbert(data);

% computes the phase of analytic signals ([-pi, pi])
phases = angle(analytics);

% calculates the indexes for each pair of sensors
for ch1 = 1: channels - 1
    for ch2 = (ch1 + 1): channels     
        % computes relative phase between two channels
        relativePhase = phases(:, ch1) - phases(:, ch2); % ([-2*pi, 2*pi])
               
        % computes PLV
        plv(ch1, ch2) = abs(mean(exp(1i * wrapTo2Pi(relativePhase))));
        plv(ch2, ch1) = plv(ch1, ch2);
        
        % computes PLI
        pli(ch1, ch2) = abs(mean(sign(wrapToPi(relativePhase))));
        pli(ch2, ch1) = pli(ch1, ch2);
        
        % computes wPLI
        [pxy,~]         = cpsd(data(: , ch1), data(: , ch2), win, olength, nfft, fs);
        ipxy            = imag(pxy);
        wpli(ch1, ch2)  = abs(mean(ipxy)) ./ mean(abs(ipxy));
        wpli(ch2, ch1)  = wpli(ch1, ch2);
        
    end
end

end