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
    subjfolder = '/Volumes/GoogleDrive/My Drive/EEG Club - Materials/2021-2022/Session4: Artefact rejection/N170 Raw Data and Scripts Only/';
else
    % Path to EEG files: select folder containing the EEG files
    fprintf('>>> Select the folder containing the EEG files\n')
    subjfolder = uigetdir('','Select the folder containing the EDF files');
end

%% Select EEG files to preprocess
filelist    = dir([subjfolder filesep '**' filesep '*N170.set']);
fprintf('>>> %g .SET files found\n',length(filelist))

%% Clean, filter and re-reference
for nF=1:length(filelist)
    
    % retrieve metadata (header)
    hdr = ft_read_header([filelist(nF).folder filesep filelist(nF).name]);
    
    % initialise configuration (cfg) structure
    cfg = [];
    
    % retrieve name of file (.set)
    cfg.dataset   = [filelist(nF).folder filesep filelist(nF).name];
    
    % indicate type of file (continuous or not)
    cfg.continuous = 'yes';
    
    % define trials to use for epoching
    cfg.trialfun                = 'eventlist_trialfun';
    cfg.trialdef.eventtype      = 'trigger';
    cfg.trialdef.eventvalue     = [201, 202];
    cfg.trialdef.prestim        = 0.5;
    cfg.trialdef.poststim       = 1;
    cfg.baselinewindow          = [-0.1 0];
    cfg                         = ft_definetrial(cfg);
    cfg.trl(:,3) = -cfg.trl(:,3);
    cfg.trl = cfg.trl([1:80, 100:180],:); %(cfg.trl(:,2)>hdr.nSamples,:)=[];
%     cfg(cfg.trl(:,1)<0,:)=[];
    
    % demean and filter the data
    cfg.demean         = 'yes';

    cfg.lpfilter       = 'yes';        % enable low-pass filtering
    cfg.lpfilttype     = 'but';
    cfg.lpfiltord      = 4;
    cfg.lpfreq         = 30;

    cfg.dftfilter      = 'no';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq        = [50 100 150];
    
%     cfg.bsfilter       = 'no';        % enable low-pass filtering
%     cfg.bsfilttype     = 'but';
%     cfg.bsfiltord      = 5;
%     cfg.bsfreq         = [45 55];
        
    % Re-referencing (over all electrodes, linked mastoids or bipolar)
    cfg.channel        = hdr.label(1:30);
    cfg.reref          = 'yes';
    cfg.refchannel     = 'all';
    
    % Preprocessing
    data                = ft_preprocessing(cfg);
    
    %%% CASE 1: no high-pass
    cfg.hpfilter       = 'no';
%     cfg.hpfilttype     = 'but';
%     cfg.hpfiltord      = 4;
%     cfg.hpfreq         = 0.1;
    data1              = ft_preprocessing(cfg,data);

    %%% CASE 2: 0.1Hz high-pass
    cfg.hpfilter       = 'yes';
    cfg.hpfilttype     = 'but';
    cfg.hpfiltord      = 4;
    cfg.hpfreq         = 0.1;
    data2              = ft_preprocessing(cfg,data);
    
    %%% CASE 3: 0.5Hz high-pass
    cfg.hpfilter       = 'yes';
    cfg.hpfilttype     = 'but';
    cfg.hpfiltord      = 4;
    cfg.hpfreq         = 0.5;
    data3              = ft_preprocessing(cfg,data);
    
    %%% CASE 4: 1Hz high-pass
    cfg.hpfilter       = 'yes';
    cfg.hpfilttype     = 'but';
    cfg.hpfiltord      = 4;
    cfg.hpfreq         = 1;
    data4              = ft_preprocessing(cfg,data);
    
    %%% CASE 5: 0.1Hz high-pass
    cfg.hpfilter       = 'yes';
    cfg.hpfilttype     = 'but';
    cfg.hpfiltord      = 4;
    cfg.hpfreq         = 0.1;
    cfg.hpfiltdir      = 'onepass';
    data5              = ft_preprocessing(cfg,data);
    
    % Average across trials
    cfg2            = [];
    cfg2.keeptrials = 'yes';
    cfg2.latency    = [-0.15 0.75];
    av_data1         = ft_timelockanalysis(cfg2, data1);
    av_data2         = ft_timelockanalysis(cfg2, data2);
    av_data3         = ft_timelockanalysis(cfg2, data3);
    av_data4         = ft_timelockanalysis(cfg2, data4);
    av_data5         = ft_timelockanalysis(cfg2, data5);
    
    % Remove trials with 20% of channels with absolute amplitude over 100uV
    clean_erp1(nF,:,:) = squeeze(mean(av_data1.trial(mean((max(abs(av_data1.trial),[],3)>100)')<0.2,:,:),1));
    clean_erp2(nF,:,:) = squeeze(mean(av_data2.trial(mean((max(abs(av_data2.trial),[],3)>100)')<0.2,:,:),1));
    clean_erp3(nF,:,:) = squeeze(mean(av_data3.trial(mean((max(abs(av_data3.trial),[],3)>100)')<0.2,:,:),1));
    clean_erp4(nF,:,:) = squeeze(mean(av_data4.trial(mean((max(abs(av_data4.trial),[],3)>100)')<0.2,:,:),1));
    clean_erp5(nF,:,:) = squeeze(mean(av_data5.trial(mean((max(abs(av_data5.trial),[],3)>100)')<0.2,:,:),1));
end


%%
figure;
hp=[];
hp(1)=plot(av_data1.time,squeeze(mean(clean_erp1(:,match_str(data.label,'PO8'),:),1)),'Color',[0 0 1],'LineWidth',2);
hold on;
hp(2)=plot(av_data1.time,squeeze(mean(clean_erp2(:,match_str(data.label,'PO8'),:),1)),'Color',[1 0 0],'LineWidth',2);
hp(3)=plot(av_data1.time,squeeze(mean(clean_erp3(:,match_str(data.label,'PO8'),:),1)),'Color',[0.8 0 0],'LineWidth',2);
hp(4)=plot(av_data1.time,squeeze(mean(clean_erp4(:,match_str(data.label,'PO8'),:),1)),'Color',[0.6 0 0],'LineWidth',2);

hp(5)=plot(av_data1.time,squeeze(mean(clean_erp5(:,match_str(data.label,'PO8'),:),1)),'Color',[0.8 0 0],'LineWidth',2,'LineStyle','--');

set(gcf,'Color','w')
set(gca,'FontSize',18,'FontWeight','bold')

xlim([-0.1 0.7])
xlabel('Time from Stimulus Onset (s)')
ylabel('Voltage (\muV)')

line([0 0],ylim,'Color',[1 1 1]*0.5,'LineStyle',':')
line(xlim,[0 0],'Color',[1 1 1]*0.5,'LineStyle',':')

legend(hp,{'no HP','0.1Hz','0.5Hz','1Hz','0.1Hz*'},'Location','EastOutside')
