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
fprintf('>>> Select the folder containing the EDF files\n')
subfolder = uigetdir('','Select the folder containing the EDF files');

%% Select EDFs to plot
% Return the subject IDs from the data folder
filelist = dir([subfolder filesep '**' filesep '*.edf']);
pick=listdlg('ListString',{filelist.name},'PromptString','Select the EDF file to check');
filelist = filelist(pick);
fprintf('>>> You have selected %g EDF files\n',length(filelist))

%% Select channels to plot
hdr = ft_read_header([filelist.folder filesep filelist.name]);
pick_channels = listdlg('ListString',hdr.label,'PromptString','Select the channels to check');
all_channels  = hdr.label(pick_channels);
    
%% display data without preprocessing
cfg=[];
cfg.dataset         = [filelist.folder filesep filelist.name];
cfg.channel         = all_channels;
cfg.continuous      = 'yes';
cfg.allowoverlap    = 'true';
cfg.viewmode        = 'vertical';
cfg.blocksize       = 30;
cfg.ylim            = 'maxmin';
cfg                 = ft_databrowser(cfg);

