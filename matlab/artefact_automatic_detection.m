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

%% EOG
cfg            = [];
cfg.datafile   = [filelist.folder filesep filelist.name];
cfg.continuous = 'yes';

cfg.trialfun                = 'eventlist_trialfun';
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = [201, 202];
cfg.trialdef.prestim        = 0.5;
cfg.trialdef.poststim       = 1;

cfg.dataset                 = [filelist.folder filesep filelist.name];
cfg                         = ft_definetrial(cfg);


% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = '*EOG*';
cfg.artfctdef.zvalue.cutoff      = 6;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [2 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% feedback
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_EOG] = ft_artifact_zvalue(cfg);

% %% Jump
% cfg.artfctdef.zvalue=rmfield(cfg.artfctdef.zvalue,'artifact');
% % channel selection, cutoff and padding
% cfg.artfctdef.zvalue.channel = 1:30;
% cfg.artfctdef.zvalue.cutoff = 20;
% cfg.artfctdef.zvalue.trlpadding = 0;
% cfg.artfctdef.zvalue.artpadding = 0;
% cfg.artfctdef.zvalue.fltpadding = 0;
% 
% % algorithmic parameters
% cfg.artfctdef.zvalue.cumulative = 'yes';
% cfg.artfctdef.zvalue.medianfilter = 'yes';
% cfg.artfctdef.zvalue.medianfiltord = 9;
% cfg.artfctdef.zvalue.absdiff = 'yes';
% 
% % make the process interactive
% cfg.artfctdef.zvalue.interactive = 'yes';
% 
% [cfg, artifact_jump] = ft_artifact_zvalue(cfg);
% 

%% muscle
cfg.artfctdef.zvalue=rmfield(cfg.artfctdef.zvalue,'artifact');
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel      = 1:30;
cfg.artfctdef.zvalue.cutoff       = 10;
cfg.artfctdef.zvalue.trlpadding   = 0;
cfg.artfctdef.zvalue.fltpadding   = 0;
cfg.artfctdef.zvalue.artpadding   = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter     = 'yes';
cfg.artfctdef.zvalue.bpfreq       = [110 140];
cfg.artfctdef.zvalue.bpfiltord    = 9;
cfg.artfctdef.zvalue.bpfilttype   = 'but';
cfg.artfctdef.zvalue.hilbert      = 'yes';
cfg.artfctdef.zvalue.boxcar       = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_muscle] = ft_artifact_zvalue(cfg);
