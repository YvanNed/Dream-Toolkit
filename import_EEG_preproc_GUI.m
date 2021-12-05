
%%%%%
% Segmenting data in 30-s epochs and adding scoring labels 

%% Initialise paths and toolboxes

clear;
close all;
set(0,'DefaultUIControlFontSize',16);

% Path to EDF files: select folder containing the EDF files
fprintf('>>> Select the folder containing the EDF files\n')
datadir = uigetdir('','Select the folder containing the EDF files');

%% Select EDFs to check

% Return the subject IDs from the data folder
filelist = dir([datadir filesep '**' filesep '*.edf']);
pick=listdlg('ListString',{filelist.name},'PromptString','Select the EDF file to check');
filelist = filelist(pick);
fprintf('>>> You have selected %g EDF files\n',length(filelist))

ft_defaults;

%% loop on subjects

for nF=1:length(filelist)
    
    % Which subject
    file_name = filelist(nF).name;
    folder_name = filelist(nF).folder;
    seps=findstr(file_name,'.');
    SubID=file_name(1:seps(1)-1);
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(filelist))
    hdr=ft_read_header([folder_name filesep file_name]);

    % Define epochs
    cfg=[];
    cfg.trialfun                = 'contseg_trialfun';
    cfg.SubID                   = SubID;
    cfg.dataset                 = [folder_name filesep file_name];
    cfg.trialdef.lengthSegments = 30; % in seconds;
    cfg = ft_definetrial(cfg);

    % Selecting channels
    pick_channels = listdlg('ListString',hdr.label,'PromptString','Select the channels to check');
    all_channels  = hdr.label(pick_channels);
    Num_ch = numel(all_channels);
    
    % Filtering
    cfg.channel        = cellstr(all_channels);
    cfg.demean         = 'yes';
    cfg.lpfilter       = 'yes';        % enable low-pass filtering
    cfg.lpfilttype     = 'but';
    cfg.lpfiltord      = 4;
    cfg.lpfreq         = 40;
    cfg.hpfilter       = 'yes';        % enable high-pass filtering
    cfg.hpfilttype     = 'but';
    cfg.hpfiltord      = 4;
    cfg.hpfreq         = 0.1;
    cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq        = [50 100 150]; % set up the frequencies for notch filtering

    % Re-referencing (over all electrodes, linked mastoids or bipolar)
    while true
        
        cfg.reref = 'yes';
        
        msg = "Choose the referencing method";
        opts = ["Average over all channels" "Linked mastoids" "Bipolar referencing"];
        choice = menu(msg,opts);
        
        switch choice
            case 1
                Refmethod = {'Average over all channels'};
                cfg.refmethod = 'avg';
                cfg.refchannel = hdr.label;  % all channels, incl those not picked earlier
                break
            case 2
                Refmethod = {'Linked mastoids'};
                cfg.refmethod = 'avg';
                if numel(find(contains(cfg.channel,{'A1', 'A2'}))) == 2
                    cfg.refchannel = {'A1', 'A2'}; 
                elseif numel(find(contains(cfg.channel,{'M1', 'M2'}))) == 2
                    cfg.refchannel = {'M1', 'M2'}; 
                else
                    warning('No mastoids found amongst the channels. Please choose another referencing method')
                    continue
                end
                break
            case 3
                Refmethod = {'Bipolar'};
                cfg.reref = 'yes';
                cfg.refmethod = 'bipolar';
                cfg.refchannel = hdr.label;  % all channels, incl those not picked earlier
                break
        end
    
    end

    % Display parameters
    Table_param = table(cellstr(strjoin(cfg.channel)),cfg.hpfreq,cfg.lpfreq,cfg.dftfreq,Refmethod,'VariableNames',{'Channels','high-pass frequency','low-pass frequency','Notch filtering','Referencing method'});
    disp(Table_param);
    
    % Preprocessing
    dat                = ft_preprocessing(cfg); 
    
    % Resampling
    cfgbs=[];
    cfgbs.resamplefs   = 256;
    cfgbs.detrend      = 'no';
    cfgbs.demean       = 'yes';
    data               = ft_resampledata(cfgbs,dat); 
       
    % Save the structure
    % save([path_data filesep 'f_ft_' file_name(1:end-4)],'data');
    
    if length(filelist)>1
        fprintf(1,'Press SPACE to continue to next EDF file.\n')
        pause;
    end
    
end





