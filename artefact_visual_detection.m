%% Initialise paths and toolboxes
clear;
close all;
set(0,'DefaultUIControlFontSize',16);

if exist('ft_read_data.m')==0
    warning('You need to add fiedltrip to your path!')
    fprintf('>>> Select the fieldtrip main folder\n')
    ft_folder = uigetdir('','Select the fieldtrip main folder');
    addpath(ft_folder)
    ft_defaults;
end


% Path to EDF files: select folder containing the EDF files
fprintf('>>> Select the folder containing the EEG files\n')
subfolder = uigetdir('','Select the folder containing the EDF files');

%% Select EEG files to plot
% Return the subject IDs from the data folder
filelist    = dir([subfolder filesep '**' filesep '*.set']);
pick        = listdlg('ListString',{filelist.name},'PromptString','Select the EDF file to check');
filelist    = filelist(pick);
fprintf('>>> You have selected %g files\n',length(filelist))

%% Import .SET (EEGlab) data into fieldtrip structure

cfg=[];
cfg.trialfun                = 'eventlist_trialfun';
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = [201, 202];
cfg.trialdef.prestim        = 0.5;
cfg.trialdef.poststim       = 1;

cfg.dataset                 = [filelist.folder filesep filelist.name];
cfg                         = ft_definetrial(cfg);
% Preprocessing
cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
cfg.dftfreq        = [50 100 150];
cfg.demean                  = 'yes';
data_preproc                = ft_preprocessing(cfg);

%% reject trials
cfg          = [];
cfg.method   = 'summary';
cfg.alim     = 5e-5;
data         = ft_rejectvisual(cfg,data_preproc);

