function surrogate = infans_surrogate(data, option)
% =========================================================================
%
% This function makes surrogate data for evaluating the statistical
% significance of a connectivity index.
%
% INPUTS:
%   data   : a matrix containing EEG recordings (samples * channels)
%   type   : the type of surrogate data (shuffle or phase)
%
% OUTPUT:
%   surrogate : a matrix surrogates (samples * channels)
%
% NOTE: 
%   1- for COR, xCOR, and COH option = 'shuffle' 
%   2- for PLV, PLI, and WPLI option = 'phase'
%
% =========================================================================

% converts the data to single precision float
data = single(data);
    
% makes surrogates
if strcmp(option, 'phase')
    % calculates FFT of each column
    ft_data = fft(data);
    
    % randomizes the phase of the data
    semiphase = rand(floor((size(ft_data, 1) - 1) / 2), size(data, 2), 'single' );
    
    % creates a symmetric conjugate phase
    phase = zeros(size(ft_data));   
    phase(2: size(semiphase, 1) + 1, :) = semiphase;
    phase(end - size(semiphase, 1) + 1: end, :) = flip(semiphase, 1);
    
    % addes of the new phase in the signal.
    ft_surrogate = abs(ft_data) .* exp(1i * 2 * pi * phase);
    
    % calculates inverse FFT
    surrogate = ifft(ft_surrogate, [], 1, 'symmetric');
    
% If the option is to perform a suffling
else
    surrogate = zeros(size(data), 'single');
    
    % shuffles of each signal.
    for channels = 1:size(data, 2)
        surrogate(:, channels) = data(randperm(size( data, 1 )), channels);
    end
end

end