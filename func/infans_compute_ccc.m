function ccc = infans_compute_ccc(data)
% =========================================================================
%
% This function calculates the circular correlation coeffcient matrix.
%
% INPUTS:
%   data : a matrix containing EEG recordings (samples * channels)
%
% OUTPUT:
%   ccc : a matrix containing circular correlation coeffcients (channel * channel)
% 
% REF: 
%   Payam Shahsavari Baboukani, et. al,
%   "A novel multivariate phase synchrony measure: Application to multichannel newborn EEG analysis"
%   https://doi.org/10.1016/j.dsp.2018.08.019.
% =========================================================================

% gets the size of the data
[samples, channels] = size(data);

% setting parameters
ip = zeros(channels, samples);

% decomposes each channel to 4 sub bands 
% for ch = 1:channels
%     % decomposes using SWT
%     [approx, details]  = swt(data(:, ch), 6, 'sym2'); 
%     
%     % selects 4 conventional EEG rhythms
%     y = [approx(6,:); details(6,:); details(5,:); details(4,:)]; 
%     decomposed(ch, :, :) = y;
% end
    
% calculates instantaneous phases
for ch = 1:channels 
    sig       = data(:,ch);
    z         = hilbert(sig);
    temp      = unwrap(angle(z));
    ip(ch, :) = temp(1:end);
end
    
% calculates circular correlation coefficient
ccc = zeros(channels);  % memory alocation for circular correlation coeffcient matrix
seg = ip';   
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
        ccc(ch1, ch2) = num / den;
    end   
end 

end