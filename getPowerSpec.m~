
%%% Input required: 
% - segmented data (e.g., 30s epoch) — output of import_EEG_preproc_GUI
% - scoring labels
% Ouput can be plotted using plotPowerSpec

clear;
close all;

%% Set paths

path = '/Users/nico/Documents/ICM/Iceberg/Data/healthy';
path_data = [path filesep 'filt_reref_seg'];
path_labels = [path filesep 'Hypno'];
files=dir([path_data filesep '*.mat']);
ft_defaults;

%% loop on subjects

for nF=1:length(files)
 %%   
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    separator = strfind(SubID,'_');
    Sub = SubID(1:separator(1)-1);

    fprintf('... working on Subject %s (%g/%g)\n',string(Sub),nF,length(files))        
        
    % Import the data
    load([folder_name filesep SubID])
    scoringlabels = readcell([path_labels filesep sprintf('%sS.TXT',Sub)]);
    % scoringlabels(end+1) = {'W'};

    data_pow = [];
    max_signal = [];
    epoch_score = [];
    
    for nCh=1:length(data.label)
        for nTr=1:length(data.trial)-1    % exclude the last epoch, which may be shorter than 30s
            fprintf('channel %1.0f/%1.0f, trial %1.0f/%1.0f\n',nCh,length(data.label),nTr,length(data.trial)-1)
            w_window                        = 6*data.fsample;
            w_overlap                       = w_window/2;
            df                              = 0.2;
            freqV                           = 1:0.2:60;
            signal                          = data.trial{nTr}(nCh,:);
            [pow,freqs]                     = pwelch(signal,w_window,w_overlap,freqV,data.fsample,'psd');
            data_pow(nTr,nCh,:)             = 10*log10(pow);
            max_signal(nTr,nCh)             = max(abs(signal));
            epoch_score{nTr}                = scoringlabels{nTr};
        end
    end  

    data_labels=data.label;
    save([path filesep 'Pow_f_ft' filesep [Sub '_Pow_f_ft']],'data_pow','max_signal','epoch_score','data_labels','freqs');

    
end