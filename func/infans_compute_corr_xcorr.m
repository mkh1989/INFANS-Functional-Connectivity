function xcor = infans_compute_corr_xcorr(data, maxlag)
% =========================================================================
%
% This function calculates cross-correlation and Pearson correlation (maxlag = 0).
%
% INPUTS:
%   data   : a matrix containing EEG recordings (samples * channels)
%   maxlag : the maximum value for the lag between two channels
%
% OUTPUT:
%   If maxlag = 0 (Pearson Correlation)
%   index  : a matrix containing corr values (channels * channels)
%
%   If maxlag ~= 0 (Auto Correlation & Cross Correlation)
%   index  : a 3D array containing cross correlation values (channels * channels * lag)
%   lag = -maxlag : maxlag
%   
%   Example:
%   maxlag = 1 ---> lag = -1, 0, 1 ---> size(index) = [channels, channels, 3];
%   index(:,:,1)  for lag = -1
%   index(:,:,2)  for lag = 0
%   index(:,:,3)  for lag = -1
%
% =========================================================================

% gets the size of the data
[~, channels] = size(data);

% computes correlation and cross correlation
[c, lags] = xcorr(data, maxlag, 'coeff');

if maxlag == 0 
    % Pearson Correlation (removes self correlation)
    xcor = reshape(c, channels, channels);
else 
    % cross-correlation
    xcor = zeros(channels, channels, length(lags));
    for i = 1:length(lags)
        xcor(:,:,i) = reshape(c(i,:), channels, channels);
    end
end

end

