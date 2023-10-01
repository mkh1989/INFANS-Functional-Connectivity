%% initialization
clear; close all; clc;

dir               = '.\SET Files\';
montage           = 'average'; % 'recorded', 'laplace', 'average', 'banana'
windowLength      = 5;   % in second
overlapLength     = 4;   % in second
lagRange          = 26;  % max lag for cross correlation (in sample)
zscoreFlag        = 1;   % set 1 for zscore normalization otherwise zero
Fs                = 256; % sampling frequency
hcutoff           = 0.5; % high cut-off frequency
lcutoff           = 12;  % low cut-off frequency
filterOrder       = 7;   % the order of high- and low-pass filters
% welchWindowLegth  = 512; % window length for welch method for cross-spectrum analysis
% welchOverlapLegth = 500; % overlap for welch method for cross-spectrum analysis
measures = {'xCOR', 'COR', 'COH', 'iCOH', 'lCOH', 'PLV', 'PLI', 'wPLI'};
outputFileName = [montage '_' num2str(windowLength) 's_' num2str(overlapLength) 's_' num2str(hcutoff) 'Hz_' num2str(lcutoff) 'Hz.mat'];
lpButter = designfilt('lowpassiir', 'FilterOrder', filterOrder, ...
'HalfPowerFrequency', lcutoff, 'SampleRate', Fs, 'DesignMethod', 'butter');     
hpButter = designfilt('highpassiir', 'FilterOrder', filterOrder, ...
'HalfPowerFrequency', hcutoff, 'SampleRate', Fs, 'DesignMethod', 'butter');

%% loop over segments for calculation of connectivity measures
for dataID = 1:50
    % loads the data
    eegFileName = [dir 'epoch_' num2str(dataID) '.set'];
    EEGOUT      = pop_loadset(eegFileName);
    
    fprintf(['working on the epoch ' num2str(dataID) ' ...\n'])
    % changes the montage
    EEGOUT = infans_change_montage(EEGOUT, montage);

    % filters the data
    EEG = filtfilt(hpButter,filtfilt(lpButter,EEGOUT.newmontage'))';
    
    % segments the data
    segmentedEEG = infans_segment_data(EEG, windowLength, overlapLength, Fs);
    
    % loop over segments
    [segments, ~, ~] = size(segmentedEEG);
    upd_connectivity = textprogressbar(segments, 'startmsg', 'Connectivity', 'endmsg', 'Done');
    for seg = 1:segments  
        segment(:,:)  = squeeze(segmentedEEG(seg,:,:))'; % (samples * channels)
        if zscoreFlag
            segment = zscore(segment); % (samples * channels)
        end

        % COR and xCOR
        c = infans_compute_corr_xcorr(segment, lagRange);
        FCs.xCOR{dataID, 1}(:,:,seg) = max(abs(c), [], 3) - eye(size(c,1));
        FCs.COR{dataID, 1}(:,:,seg)  = abs(c(:, :, lagRange + 1)) - eye(size(c,1));

        % COH, iCOH, lCOH
        [coh, icoh, lcoh] = infans_compute_coh_icoh_lcoh(segment, Fs);
        FCs.COH{dataID, 1}(:,:,seg)  = coh;
        FCs.iCOH{dataID, 1}(:,:,seg) = icoh;
        FCs.lCOH{dataID, 1}(:,:,seg) = lcoh;
        
        % PLV, PLI, wPLI
        [plv, pli, wpli] = infans_compute_plv_pli_wpli(segment, Fs);
        FCs.PLV{dataID, 1}(:,:,seg)  = plv;
        FCs.PLI{dataID, 1}(:,:,seg)  = pli;
        FCs.wPLI{dataID, 1}(:,:,seg) = wpli;
        
        % CCC
        % ccc = infans_compute_ccc(segment);
        % Index.CCC{dataID, 1}(:,:,seg) = ccc;
        
        % updates the progress bar
        upd_connectivity(seg);
    end
    fprintf('-------------------------------\n')
    fprintf('\n')
end
save(outputFileName, 'FCs')

%% thresholding 
 mean_measure = zeros(1, length(measures));
 std_measure  = zeros(1, length(measures));
 thresh       = zeros(1, length(measures));
 for meas = 1:length(measures)
    concatenatedMatrix = cat(3, FCs.(measures{meas}){:});  
    mean_measure(meas) = mean(concatenatedMatrix, 'all');
    std_measure(meas)  = std(concatenatedMatrix, 0, 'all');
    thresh(meas)       = mean_measure(meas) + std_measure(meas);

    clear concatenatedMatrix
    
    for dataID = 1:50
        Binarized.(measures{meas}){dataID, 1} = zeros(size(FCs.xCOR{1, 1}));       
        Binarized.(measures{meas}){dataID, 1} = double(FCs.(measures{meas}){dataID, 1} > thresh(meas));
    end
 end
save(outputFileName, 'Binarized', '-append')
 
%% graph analysis
% channels = size(Index.xCOR.Raw{1, 1}, 1);
for dataID = 1:50
    fprintf(['working on the epoch ' num2str(dataID) ' ...\n'])
    for meas = 1:length(measures)
        fprintf(['Measure ' cell2mat(measures(meas)) ' ...\n'])
        segments  = size(Binarized.(measures{meas}){dataID, 1}, 3);
        upd_graph = textprogressbar(segments, 'startmsg', 'Graph Analysis', 'endmsg', 'Done');
        
        for seg = 1:segments
            
            tempGraph(:,:) = Binarized.(measures{meas}){dataID, 1}(:, :, seg);
            
            % calculates global efficiency
            GE = efficiency_bin(tempGraph);
            
            % calculates local efficiency
            LE = mean(efficiency_bin(tempGraph, 1));
            
            % calculates modularity
            if sum(nnz(tempGraph)) ~= 0
                [Ci, Q] = modularity_und(tempGraph);
                MOD = Q;
            else
                MOD = NaN;
            end
            
            % claculates mean clustering coefficient
            MCCOEFF = mean(clustering_coef_bu(tempGraph));
            
            % calculates mean closeness centrality
            G = graph(tempGraph);
            MCLCEN = mean(centrality(G, 'closeness'));
            
            % calculates average degree
            AD = mean(degrees_und(tempGraph));
            
            % saves the metrics in corresponding location
            Metrics.(measures{meas}){dataID, 1}.GlobalEfficiency(seg, 1) = GE;
            Metrics.(measures{meas}){dataID, 1}.LocalEfficiency(seg, 1)  = LE;
            Metrics.(measures{meas}){dataID, 1}.Modularity(seg, 1)       = MOD;
            Metrics.(measures{meas}){dataID, 1}.MeanClustCoeff(seg, 1)   = MCCOEFF;
            Metrics.(measures{meas}){dataID, 1}.MeanCloseCent(seg, 1)    = MCLCEN;
            Metrics.(measures{meas}){dataID, 1}.AverageDegree(seg, 1)    = AD;
            
            % updates the progress bar
            upd_graph(seg);
        end
    end
    
    fprintf('-------------------------------\n')
    fprintf('\n')
end
save(outputFileName, 'Metrics', '-append')