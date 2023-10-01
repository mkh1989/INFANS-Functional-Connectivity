function epochs = infans_segment_data(data, W, O, fs)
% Inputs:
% data         : channel * number of time samples
% label        : a row vector containing labels for each time point (0: nonseizure, 1:seizure)
% W            : length of each epoch in second
% O            : length of overlap in second
% fs           : sampling frequency

% Outputs:
% epochs:      : a 3D matrix containing eeg epochs (number of epochs * number of channels * number of time samples)
% epoch_labels : a column vector containing the label of each epoch (0: nonseizure, 1:seizure) 

    % initialize
    W = W * fs;
    O = O * fs;
    segments = fix((size(data,2) - O) / (W - O));
    
    % segment data
    epochs = zeros(segments, size(data,1), W);
    for s = 1: segments
        epochs(s, :, :)   = data(:, (s - 1) * (W - O) + (1: W));
    end

end