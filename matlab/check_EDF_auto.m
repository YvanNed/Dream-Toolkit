

%% EDF format checks

% There are 3 potential issues related to the conversion of raw signals to EDF:

% 1. Signal clipping:    signal cut once it exceeds an amplitude threshold
%                        (the min-max range set before EDF conversion was
%                        too narrow)
% 2. Bit depth:          signal shows a stair-like progression (the min-max
%                        range set before EDF conversion was too wide)
% 3. Inverted polarity:  signal multiplied by -1

% The present script is a fully automatic version of the GUI: loops across 
% all datasets and stores a table that contains information about signal 
% clipping and bit depth. Returns the channels whose signal may be clipped 
% or appear in low resolution, and offers the possibility to visually 
% inspect them.

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

% Return the subject IDs from the data folder
filelist = dir([subfolder filesep '**' filesep '*.EDF']);
% filelist = dir([subfolder filesep '*.edf']);
fprintf('>>> %s EDF files found\n',string(numel(filelist)))

%% Loop across subjects

summary_table=array2table(zeros(2000,6),'VariableNames',{'File','Channel','Unit','Min','Max','BinGap'});   % 2000 is an arbitrary number for preallocation
summary_table.File=categorical(summary_table.File);
summary_table.Channel=categorical(summary_table.Channel);
summary_table.Unit=categorical(summary_table.Unit);
nc=0;

for S = 1:length(filelist)

    % Parameters subject
    subID = filelist(S).name;
    Sub = subID(1:end-4);
    subfolder = filelist(S).folder;

    % Import the data
    cfg = [];
    cfg.dataset = [subfolder filesep subID];
    
    fprintf(1,'>>> >>> Importing data from Subject %s... (%g/%g),\n',Sub,S,size(filelist,1))
    data = ft_read_data(cfg.dataset);
    hdr = ft_read_header(cfg.dataset);
    all_channels  = hdr.label;
    
    % Store data; will be retrieved for visual inspection later 
    KeepData(S,:) = [{subID} {data}];
    
    %% Is the data preprocessed already? Plot the power spectrum
    
%     % Get scoring labels
%     fprintf('>>> Select the folder of the excel file containing the scoring labels\n')
%     filedir = uigetdir('','Select the folder of the text file containing the scoring labels');
%     scoringlabels = readcell([filedir filesep sprintf('%sS.TXT',Sub)]);
%   
%     % Define epochs 
%     cfg.trialfun                = 'contseg_trialfun';
%     cfg.trialdef.lengthSegments = 30; % in seconds;
%     cfg = ft_definetrial(cfg);
%     Dat = ft_preprocessing(cfg); 
%     
%     data_pow = [];
%     max_signal = [];
%     epoch_score = [];
% 
%     % Compute spectral power
%     for nCh=1:length(all_channels)
%         for nTr=1:length(Dat.trial)-1    % exclude the last epoch, which may be shorter than 30s
%             fprintf('channel %1.0f/%1.0f, trial %1.0f/%1.0f\n',nCh,length(Dat.label),nTr,length(Dat.trial)-1)
%             w_window                        = 6*Dat.fsample;
%             w_overlap                       = w_window/2;
%             df                              = 0.2;
%             freqV                           = 1:0.2:60;
%             signal                          = Dat.trial{nTr}(nCh,:);
%             [pow,freqs]                     = pwelch(signal,w_window,w_overlap,freqV,Dat.fsample,'psd');
%             data_pow(nTr,nCh,:)             = 10*log10(pow);
%             max_signal(nTr,nCh)             = max(abs(signal));
%             epoch_score{nTr}                = scoringlabels{nTr};
%         end
%     end 
%     
%     % Plot and save power spectrum
%     plotPowerSpec(subfolder,Sub,data_pow,max_signal,epoch_score,freqs,all_channels);
        
    %% Check for signal clipping and bit depth issue
       
    for i = 1:length(all_channels)  
        
        nc = nc+1;
        
        % calculates the difference in amplitude between neighoring data points
        Data = data(i,:);
        delta_ampl = abs(diff(Data));    
        delta_ampl(delta_ampl==0) = [];   % ignore the peak at 0 uV
        
        if isempty(delta_ampl)  % For channels with no signal, this line allows to keep running the code
            delta_ampl = 0;
        end
        
        % Prepare the summary table
        table_mat=[min(Data) max(Data) min(delta_ampl)];
        summary_table.File(nc) = subID;
        summary_table.Channel(nc) = all_channels{i};
        summary_table.Unit(nc) = hdr.orig.PhysDim(hdr.orig.chansel(i),:);
        summary_table.Min(nc) = table_mat(1);
        summary_table.Max(nc) = table_mat(2);
        summary_table.BinGap(nc) = table_mat(3);
        
        % Get the first 2 values of the normal distribution (will be used
        % to identify signal clipping peaks)
        HIST = histcounts(data((i),:),'Normalization','count');
        Min1 = HIST(1); 
        if delta_ampl==0  % For channels with no signal, this line allows to keep running the code
            Min2 = 0;
        else
            Min2 = HIST(2);
        end
        MIN(nc,:) = [Min1 Min2]; 
                
    end
end

emptyrows = find(summary_table.File == '0');
summary_table(emptyrows,:) = [];

% Separate channels in mv / uV (threshold values change accordingly)
mVchan = find(summary_table.Unit == 'mV');
uVchan = find(summary_table.Unit == 'uV');

if isempty(mVchan) && isempty(uVchan)
    warning('Channel unit unknown, therefore files could not be checked. Please indicate the unit for each channel (''mV'' or ''uV'') in hdr.orig.PhysDim') 
    return
end

%% Store and inspect problematic files: Signal clipping and Bit Depth files

% Signal Clipping (SC): if distribution within -500 - +500 uV and if peak
% is 20 times higher than the neighboring data point.
SC_files_idx = sort([uVchan(summary_table.Min(uVchan)>-500 | summary_table.Max(uVchan)<500); ...
                     mVchan(summary_table.Min(mVchan)>-0.50 | summary_table.Max(mVchan)<0.50)]); 
nopeak = find(MIN(SC_files_idx,1)<(MIN(SC_files_idx,2)+1)*20);   % If the first value of the distribution is not 20 times bigger than the second value, don't consider this a peak --> no clipping.
SC_files_idx(nopeak) = [];            
SC_files = summary_table(SC_files_idx,1:5);  

if ~isempty(SC_files)
    fprintf('>>> %s channels (from %s EDF files) might have a signal clipping issue.\n\n',string(size(SC_files,1)),string(numel(unique(SC_files(:,1)))));
    disp(SC_files)
else
    fprintf('>>> No signal clipping issue was detected.\n');
end

% Bit Depth (BD): if the bin difference is greater than 0.1 uV 
BD_files_idx = sort([uVchan(summary_table.BinGap(uVchan)>0.1);mVchan(summary_table.BinGap(mVchan)>0.001)]); 
BD_files = [summary_table(BD_files_idx,1:3) summary_table(BD_files_idx,6)];

if ~isempty(BD_files)
     fprintf('>>> %s channels (from %s EDF files) might have a bit depth issue.\n\n',string(size(BD_files,1)),string(numel(unique(BD_files(:,1)))))
     disp(BD_files)
else
    fprintf('>>> No resolution issue was detected.\n');
end

%% Save summary table, SC and BD files

path_summary = [subfolder filesep 'SummaryTable'];
if exist(path_summary,'file')==0
    mkdir(path_summary)
end
save([path_summary filesep 'SummaryTables'],'summary_table','SC_files','BD_files');

%% Visual inspection of problematic files: plot histograms

if ~isempty(~BD_files_idx) || ~isempty(SC_files_idx)
    reply = input('Do you want to visually inspect these channels? Press ''y'' (YES) or ''n'' (NO) > ','s');
    if reply == 'y'
        fprintf('>>> Ok, plotting histograms for these %s channels...\n',string(size(SC_files,1)+size(BD_files,1)));

        % Plot SC
        for S = 1:size(SC_files,1)
            cfg.dataset = [subfolder filesep char(SC_files.File(S))];
            hdr = ft_read_header(cfg.dataset);
            whichdata = find(contains(string(char(KeepData(:,1))),string(SC_files.File(S))));
            dat = KeepData{whichdata,2};
            all_channels = hdr.label;
            the_channel = char(SC_files.Channel(S));
            % Plot figure
            f=figure;ax=gca;
            chan_idx = find(ismember(all_channels, string(the_channel)));
            Data = dat(chan_idx,:);
            if hdr.orig.PhysDim(hdr.orig.chansel(chan_idx),1:2) == 'uV'
                H = histogram(Data,-max(abs(Data)):0.005:max(abs(Data)),'EdgeColor','#1167b1'); 
            else
                H = histogram(Data,-max(abs(Data)):0.000005:max(abs(Data)),'EdgeColor','#1167b1'); 
            end
            xlim([-1.05 1.05]*max(abs(Data)))
            Max2 = sort(H.Values(2:end-1), 'descend');
            if Max2~=0
                ylim([0 (Max2(2)*1.3)])  % Take a bin next to 0 as the ylimit, add 30% above
            end
            t = title({sprintf('%s (%s)',char(SC_files.File(S)),the_channel);'Amplitude distribution'},'Interpreter','none');
            t.FontWeight = 'normal';
            xlabel(sprintf('Amplitude (%s)',hdr.orig.PhysDim(hdr.orig.chansel(chan_idx),1:2))); ylabel('Data points distribution')
            ax.FontSize = 14;
        end

        % Plot BD
        for S = 1:size(BD_files,1)
            cfg.dataset = [subfolder filesep char(BD_files.File(S))];
            hdr = ft_read_header(cfg.dataset);
            whichdata = find(contains(string(char(KeepData(:,1))),string(BD_files.File(S))));
            dat = KeepData{whichdata,2};
            all_channels = hdr.label;
            the_channel = char(BD_files.Channel(S));
            % Plot figure
            f=figure; ax=gca;
            chan_idx = find(contains(all_channels, string(the_channel)));
            Data = data(chan_idx,:);
            delta_ampl = abs(diff(Data));
            H2=histogram(delta_ampl,0:0.01:50,'EdgeColor','#1167b1');
            xlabel(sprintf('Delta amplitude (%s)',hdr.orig.PhysDim(hdr.orig.chansel(chan_idx),1:2))); ylabel('Data points distribution') 
            t = title({sprintf('%s (%s)',char(BD_files.File(S)),the_channel);'Absolute difference in amplitude between neighboring data points'},'Interpreter','none');
            t.FontWeight = 'bold';
            ax.FontSize = 14;
            ylim([0 1]*max(H2.Values(H2.BinEdges(1:end-1)>=0.01)))
        end
    end
end




