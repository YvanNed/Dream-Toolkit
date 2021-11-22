

%% EDF format checks

% There are 3 potential issues related to the conversion of raw signals to EDF:

% 1. Signal clipping:    signal cut once it exceeds an amplitude threshold
%                        (the min-max range set before EDF conversion was
%                        too narrow)
% 2. Bit depth:          signal shows a stair-like progression (the min-max
%                        range set before EDF conversion was too wide)
% 3. Inverted polarity:  signal multiplied by -1

% The present script is an fully automatic version of the semi-automatic 
% code: loops across all datasets and generates a table that summarises the
% issues (if any) detected in each dataset and channel.

%% Initialise paths and toolboxes

clear;
close all;
set(0,'DefaultUIControlFontSize',16);

% Path to EDF files: select folder containing the EDF files
fprintf('>>> Select the folder containing the EDF files\n')
subfolder = uigetdir('','Select the folder containing the EDF files');

% Return the subject IDs from the data folder
filelist = dir([subfolder filesep '**' filesep '*.edf']);
fprintf('>>> %s EDF files found\n',string(numel(filelist)))

%% Loop across subjects

summary_table=array2table(zeros(0,6),'VariableNames',{'File','Channel','Unit','Min','Max','BinGap'});
summary_table.File=categorical(summary_table.File);
summary_table.Channel=categorical(summary_table.Channel);
summary_table.Unit=categorical(summary_table.Unit);
nc=0;

%for S = 1:length(filelist)
for S = 1:25

    % Parameters subject
    subID = filelist(S).name;
    Sub = subID(1:end-4);
    subfolder = filelist(S).folder;
    
    % Import the data
    cfg = [];
    cfg.dataset = [subfolder filesep subID];
    
    fprintf(1,'>>> >>> Importing data from Subject %s...\n',Sub)
    data = ft_read_data(cfg.dataset);
    hdr = ft_read_header(cfg.dataset);
    all_channels  = hdr.label;
        
    %% Check for signal clipping and bit depth issue
    
    for i = 1:length(all_channels)
        
        nc = nc+1;
        Data = data((i),:);
        
        delta_ampl = abs(diff(Data));    % calculates the difference in amplitude between neighoring data points
        delta_ampl(delta_ampl==0) = [];  % ignore the 0 peak
        
        table_mat=[min(Data) max(Data) min(delta_ampl)];
        summary_table.File(nc) = subID;
        summary_table.Channel(nc) = all_channels{i};
        summary_table.Unit(nc) = hdr.orig.PhysDim(i,:);
        summary_table.Min(nc) = table_mat(1);
        summary_table.Max(nc) = table_mat(2);
        summary_table.BinGap(nc) = table_mat(3);
        
    end
end

%% Store problematic files

mVchan = find(summary_table.Unit == 'mV');
uVchan = find(summary_table.Unit == 'uV');

% Signal clipping (SC) files
SC_files_idx = sort([uVchan(summary_table.Min(uVchan)>-500 | summary_table.Max(uVchan)<500) ...
                       mVchan(summary_table.Min(mVchan)>-0.050 | summary_table.Max(mVchan)<0.050)]); 
SignalClipping_files = [summary_table(SC_files_idx,1) summary_table(SC_files_idx,2)];    

if ~isempty(SignalClipping_files)
    fprintf('>>> %s EDF files might have a signal clipping issue. Please check the following channels\n\n',string(numel(unique(SignalClipping_files(:,1)))));
    disp(SignalClipping_files)
    fprintf('>>> Do you want to visually inspect the data? Press ''y'' for YES or ''n'' for NO');
    if strcmp(reply,'y')
        fprintf('>>> Plotting histograms for the %s channels where signal clipping might be present\n',string(size(SignalClipping_files,1)));
        for S = 1:size(SignalClipping_files,1)
            cfg.dataset = [subfolder filesep char(SignalClipping_files.File(S))];
            data = ft_read_data(cfg.dataset);
            hdr = ft_read_header(cfg.dataset);
            all_channels  = hdr.label;
            the_channel  = char(SignalClipping_files.Channel(S));
            % Plot figure
            f=figure;ax=gca;
            chan_idx = find(contains(all_channels, string(the_channel)));
            Data = data(chan_idx,:);
            H = histogram(Data,-max(abs(Data)):0.05:max(abs(Data)),'EdgeColor','#1167b1'); 
            xlim([-1.05 1.05]*max(abs(Data)))
            ylim([0 1]*max(H.Values(H.BinEdges(1:end-1)<=-0.05 | H.BinEdges(1:end-1)>=0.05)))
            t = title({sprintf('%s (%s)',char(SignalClipping_files.File(S)),the_channel);'Amplitude distribution'},'Interpreter','none');
            t.FontWeight = 'normal';
            xlabel(sprintf('Amplitude (%s)',hdr.orig.PhysDim(chan_idx,1:2))); ylabel('Data points distribution')
            ax.FontSize = 14;
        end
    else
        return
    end
else
    fprintf('>>> No signal clipping issue was detected.\n');
end

% Bit Depth (BD) files
BD_files_idx = sort([uVchan(summary_table.BinGap(uVchan)>0.1) mVchan(summary_table.BinGap(mVchan)>0.00001)]); 
BitDepth_files = summary_table.File(BD_files_idx);
if ~isempty(BitDepth_files)
    fprintf('>>> %s EDF files might have a bit depth issue. Please check the following channels\n\n',string(numel(unique(DitDepth_files(:,1)))));
    disp(BitDepth_files)
else
    fprintf('>>> No resolution issue was detected.\n');
end













