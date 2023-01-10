%%
clear all
close all

Folder_Name='/Users/thandrillon/Data/AC_EDF_LUSPY/';
File_Name='AC_Pilote.edf';


hdr=ft_read_header([Folder_Name filesep File_Name]);
% data=ft_read_data([Folder_Name filesep File_Name]);

evt=LUSPY_GetEvents([Folder_Name filesep File_Name(1:end-4) '.csv']);
% start recording
temp=hdr.orig.T0;
hours=temp(4);
minutes=temp(5);
secondes=temp(6);
millisecondes=0;
start_time=3600*(hours)+60*(minutes)+(secondes)+(millisecondes)/1000;
    
before_events=evt(find_trials(evt.Evenement,'f_D'),:);
for k=1:length(before_events.HeurePrecise)
    temp=before_events.HeurePrecise{k};
    temp(temp==' ')=[];
    hours=temp(1:findstr(temp,'h')-1);
    minutes=temp(findstr(temp,'h')+1:findstr(temp,'min')-1);
    sep=findstr(temp,'s');
    secondes=temp(findstr(temp,'min')+3:sep(1)-1);
    millisecondes=temp(sep(1)+1:findstr(temp,'ms')-1);
    before_events.Time(k)=3600*str2num(hours)+60*str2num(minutes)+str2num(secondes)+str2num(millisecondes)/1000-start_time;
end

after_events=evt(find_trials(evt.Evenement,'f_F'),:);
for k=1:length(after_events.HeurePrecise)
    temp=after_events.HeurePrecise{k};
    temp(temp==' ')=[];
    hours=temp(1:findstr(temp,'h')-1);
    minutes=temp(findstr(temp,'h')+1:findstr(temp,'min')-1);
    sep=findstr(temp,'s');
    secondes=temp(findstr(temp,'min')+3:sep(1)-1);
    millisecondes=temp(sep(1)+1:findstr(temp,'ms')-1);
    after_events.Time(k)=3600*str2num(hours)+60*str2num(minutes)+str2num(secondes)+str2num(millisecondes)/1000-start_time;
end

%% Check locking

%%
%%% Import data, select for breaks and cut
cfg=[];
cfg.trialfun          = 'LUSPY_Epoching';
cfg.dataset           = [Folder_Name filesep  File_Name];
cfg.demean            = 'no';
cfg.evt               = before_events;
cfg.channel           = hdr.label(ismember(hdr.chantype,'unknown') & ~strcmp(hdr.label,'Mic')); %hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG'))) & (cellfun(@isempty,regexp(hdr.label,'EMG'))) & (cellfun(@isempty,regexp(hdr.label,'ECG'))) & (cellfun(@isempty,regexp(hdr.label,'Mic')))));
cfg.reref             = 'yes';
cfg.refchannel        = {'A1','A2'};

cfg.trialdef.prestim  = -5; % in seconds
cfg.trialdef.poststim = 0; % in seconds
cfg                   = ft_definetrial(cfg);
data_before           = ft_preprocessing(cfg); % read raw data

cfg.evt               = after_events;
cfg.trialdef.prestim  = 0; % in seconds
cfg.trialdef.poststim = 5; % in seconds
cfg                   = ft_definetrial(cfg);
data_after            = ft_preprocessing(cfg); % read raw data

%%% Extract power
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          =  0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
cfg.toi         =  5;                         % analysis 2 to 30 Hz in steps of .2 Hz
cfg.t_ftimwin    =  ones(length(cfg.foi),1).*5;   % length of time window =6 sec
cfg.keeptrials   = 'yes';
cfg.pad           = 'nextpow2';
TFdata_after           = ft_freqanalysis(cfg, data_after);
TFdata_after.logpow=squeeze(mean(log10(TFdata_after.powspctrm),1));

TFdata_before           = ft_freqanalysis(cfg, data_before);
TFdata_before.logpow=squeeze(mean(log10(TFdata_before.powspctrm),1));


% cfg.keeptrials   = 'no';
% cfg.output        = 'fooof_aperiodic';
% cfg.foi          =  0.5:0.2:20;                         % analysis 2 to 30 Hz in steps of .2 Hz
% fractal = ft_freqanalysis(cfg, data);
