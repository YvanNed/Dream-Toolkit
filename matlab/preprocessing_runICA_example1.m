%%% Test the influence of high-pass filtering on ERP

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


if ~isempty(findstr(pwd,'thandrillon'))
    subjfolder = '/Volumes/GoogleDrive/My Drive/EEG Club - ICM/2021-2022/Session4: Artefact rejection/N170 Raw Data and Scripts Only/';
else
    % Path to EEG files: select folder containing the EEG files
    fprintf('>>> Select the folder containing the EEG files\n')
    subjfolder = uigetdir('','Select the folder containing the EDF files');
end

%% Select EEG files to preprocess
filelist    = dir([subjfolder filesep '**' filesep '*N170.set']);
fprintf('>>> %g .SET files found\n',length(filelist))

%% Clean, filter and re-reference
nF=1; %:length(filelist)

% retrieve metadata (header)
hdr = ft_read_header([filelist(nF).folder filesep filelist(nF).name]);

% initialise configuration (cfg) structure
cfg = [];

% retrieve name of file (.set)
cfg.dataset   = [filelist(nF).folder filesep filelist(nF).name];

% indicate type of file (continuous or not)
cfg.continuous = 'yes';

% demean and filter the data
cfg.demean         = 'yes';

cfg.lpfilter       = 'yes';        % enable low-pass filtering
cfg.lpfilttype     = 'but';
cfg.lpfiltord      = 4;
cfg.lpfreq         = 30;

cfg.hpfilter       = 'yes';        % enable low-pass filtering
cfg.hpfilttype     = 'but';
cfg.hpfiltord      = 4;
cfg.hpfreq         = 0.1;

cfg.dftfilter      = 'no';        % enable notch filtering to eliminate power line noise
cfg.dftfreq        = [50 100 150];


%% Re-referencing (over all electrodes, linked mastoids or bipolar)
cfg.channel        = hdr.label(1:30);
cfg.reref          = 'yes';
cfg.refchannel     = 'all';
data_cont          = ft_preprocessing(cfg);

%% Run ICA
rankICA = rank(data_cont.trial{1,1});
cfg        = [];
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
cfg.numcomponent = rankICA;
comp_cont = ft_componentanalysis(cfg, data_cont);

%% Construct layout from EEG labels and biosemi layout
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=data_cont.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

%% Plot components
figure;
cfg = [];
cfg.component = 1:length(comp_cont.label);       % specify the component(s) that should be plotted
cfg.layout    = layout; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
%     cfg.marker = 'labels';
ft_topoplotIC(cfg, comp_cont)
set(gcf,'Position',[1           1        1871         984]);

cfg=[];
cfg.continuous='yes';
cfg.allowoverlap='true';
cfg.viewmode='vertical';
cfg.blocksize=20;
cfg = ft_databrowser(cfg, comp_cont);

%% Reject components
cfg = [];
cfg.component = [1];
data_clean1 = ft_rejectcomponent(cfg, comp_cont, data_cont);
 
cfg.component = [1 6 5 8];
data_clean2 = ft_rejectcomponent(cfg, comp_cont, data_cont);

%% Plot comparison
figure('Position',[5         378        1391         2*419]);
subplot(2,1,1);
plot(data_cont.time{1},data_cont.trial{1}(find(ismember(data_cont.label,{'O1','O2'})),:)','k');
hold on;
plot(data_clean1.time{1},data_clean1.trial{1}(find(ismember(data_clean1.label,{'O1','O2'})),:)','r');
xlim([217.7374  235.8553])
ylim([-111.4003  246.7167])
set(gcf,'Color','w')
set(gca,'FontSize',18,'FontWeight','bold')

subplot(2,1,2);
plot(data_cont.time{1},data_cont.trial{1}(find(ismember(data_cont.label,{'O1','O2'})),:)','k');
hold on;
plot(data_clean2.time{1},data_clean2.trial{1}(find(ismember(data_clean2.label,{'O1','O2'})),:)','r');
xlim([217.7374  235.8553])
ylim([-111.4003  246.7167])
set(gcf,'Color','w')
set(gca,'FontSize',18,'FontWeight','bold')
