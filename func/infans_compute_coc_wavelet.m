function [ccc, coc] = infans_compute_coc_wavelet(data)
% =========================================================================
%
% This function calculates the circular omega complexity. 
% The different EEG frequency bands are extracted by wavelet decomposition
%
% INPUTS:
%   data : a matrix containing EEG recordings (samples * channels)
%
% OUTPUT:
%   ccc : a 3D matrix containing circular correlation coeffcient matrices
%   relating to different sub-bands (sub-band * channel * channel)
%   coc : a number indicationg the circular omega complexity corresponding to the data
% 
% =========================================================================

% gets the size of the data
[samples, channels] = size(data);

% setting parameters
nsubbands  = 4;  % the number of frequency bands
decomposed = zeros(channels, nsubbands, samples);
ip         = zeros(nsubbands, channels, samples);

% decomposes each channel to 4 sub bands 
for ch = 1:channels
    % decomposes using SWT
    [approx, details]  = swt(data(:, ch), 6, 'sym2'); 
    
    % selects 4 conventional EEG rhythms
    y = [approx(6,:); details(6,:); details(5,:); details(4,:)]; 
    decomposed(ch, :, :) = y;
end
    
% calculates instantaneous phases
for ch = 1:channels
    for subband = 1:nsubbands
        sig                = squeeze(decomposed(ch, subband, :));
        z                  = hilbert(sig);
        temp               = unwrap(angle(z));
        ip(subband, ch, :) = temp(1:end);
    end
end
    
% calculates circular correlation coefficient
cirMatrix   = zeros(channels); % circular correlation matrix
ccc         = zeros(nsubbands, channels, channels);
subband_coc = zeros(1, nsubbands);
for subband = 1:nsubbands
    seg = squeeze(ip(subband, :, :))';   
    for ch1 = 1:channels   
        for ch2 = 1:channels
            chan1 = seg(:, ch1);
            chan2 = seg(:, ch2);
            
            % computes weighted sum of cos and sin of angles
            chan1_bar = angle(sum(ones(size(chan1)) .* exp(1i*chan1), 1));
            chan2_bar = angle(sum(ones(size(chan2)) .* exp(1i*chan2), 1));
            
            % computes circular correlation coeffcient
            num = sum(sin(chan1 - chan1_bar) .* sin(chan2 - chan2_bar));
            den = sqrt(sum(sin(chan1 - chan1_bar).^2) .* sum(sin(chan2 - chan2_bar).^2));          
            cirMatrix(ch1, ch2) = num / den;
        end   
    end 
    
    % saves the circular correlation coeffcient matrix
    ccc(subband, :, :) = cirMatrix;
    
    % computes COC measure value
    eigs  = eig(cirMatrix);
    eigs  = eigs/sum(eigs);
    omega = sum(eigs .* log2(eigs));
    tempcoc   = 1 + omega/log2(channels);   
    subband_coc(subband) = abs(tempcoc);
end 

% average of coc across all sub bands
coc = mean(abs(subband_coc)); 

end