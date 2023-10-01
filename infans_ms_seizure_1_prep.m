clear all
close all
clc

% filter designing (0.50-30 Hz) and notch (50 Hz)
lpButter = designfilt('lowpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 30, 'SampleRate', 256, 'DesignMethod', 'butter');
     
hpButter = designfilt('highpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 0.5, 'SampleRate', 256, 'DesignMethod', 'butter');

bsFilt = designfilt('bandstopiir','FilterOrder', 8, ...
         'HalfPowerFrequency1', 49.5, 'HalfPowerFrequency2',50.5, ...
         'SampleRate', 256);

fileDir = '.\120-second\Raw';
saveDir = '.\120-second\Preprocessed';

for i = 1:50
    % loading files
    fileName1 = [fileDir '\epoch_' num2str(i) '.set'];
    EEG = pop_loadset(fileName1);
    
    % z-score normalization
    % EEG.data = zscore(EEG.data')';
    
    % clearing some stuff
    EEG.filepath    = []; 
    EEG.event       = [];
    EEG.etc         = [];
    EEG.urevent     = [];
    EEG.icaact      = [];
    EEG.icawinv     = [];
    EEG.icasphere   = []; 
    EEG.icaweights  = [];
    EEG.icachansind = []; 
    
    % bandpass and notch filtering
    EEG.data = filtfilt(bsFilt,double(EEG.data)')';
    EEG.data = filtfilt(hpButter,filtfilt(lpButter,EEG.data'))';
    
    % saving the data
    fileName2 = ['\epoch_' num2str(i) '.set'];
    pop_saveset(EEG, 'filename', fileName2, 'filepath', saveDir, 'savemode', 'onefile');
    clear EEG
end