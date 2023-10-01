function EEGOUT = infans_change_montage(EEG, montage)
% This function chnage the montage of EEG channels
% Four montages are available: as recorded, laplace, double-banan, common average

% INPUTS: 
%   EEG    : EEGLAB structure
%   montage : a string indicating the type of montage

% OUTPUT:
%   EEGOUT : EEGLAB structure containing signals in the new montage (channel * time)
%   Note: in the double-banana montage the number of channels is less by one
%   Note: the default montage for the output signals is "common average"

EEGOUT = EEG;

switch montage
    case 'recorded'
        EEGOUT.newmontage = double(EEG.data);
    case 'laplace'
        EEGOUT.newmontage = double(laplacian_montage(EEG.data, [EEG.chanlocs.X], [EEG.chanlocs.Y], [EEG.chanlocs.Z]));
    case 'banana'
        EEGOUT.newmontage = double(double_banana_montage(EEG.data));
    case 'average'
        EEGOUT.newmontage = common_average(EEG.data);
    otherwise
        EEGOUT.newmontage = common_average(EEG.data);
end

end