

%% EDF format checks

% There are 3 potential issues related to the conversion of raw signals to EDF:

% 1. Signal clipping:    signal cut once it exceeds an amplitude threshold
%                        (the min-max range set before EDF conversion was
%                        too narrow)
% 2. Bit depth:          signal shows a stair-like progression (the min-max
%                        range set before EDF conversion was too wide)
% 3. Inverted polarity:  signal multiplied by -1

%% Initialise paths and toolboxes

clear;
close all;
set(0,'DefaultUIControlFontSize',16);

% Path to EDF files: select folder containing the EDF files
fprintf('>>> Select the folder containing the EDF files\n')
subfolder = uigetdir('','Select the folder containing the EDF files');

%% List qll EDFs to check
% Return the subject IDs from the data folder
filelist = dir([subfolder filesep '**' filesep '*.edf']);

%% Loop across subjects
summary_table=array2table(zeros(0,6),'VariableNames',{'File','Channel','Unit','Min','Max','BinDepth'});
summary_table.File=categorical(summary_table.File);
summary_table.Channel=categorical(summary_table.Channel);
summary_table.Unit=categorical(summary_table.Unit);

for S = 1:length(filelist)
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
    
    % Pick the channels to check
    all_channels  = hdr.label;
    
    %%%  Visualise the data
    
    %% Check for signal clipping and bit depth issue
    nc=0;
    for i = 1:length(all_channels)
        nc=nc+1;
        Data = data((i),:);
        
        delta_ampl = abs(diff(Data));
        delta_ampl(delta_ampl==0)=[];
        table_mat=[min(Data) max(Data) min(delta_ampl)];
        summary_table.File(nc)=subID;
        summary_table.Channel(nc)=all_channels{i};
        summary_table.Unit(nc)=hdr.orig.PhysDim(i,:);
        summary_table.Min(nc)=table_mat(1);
        summary_table.Max(nc)=table_mat(2);
        summary_table.BinDepth(nc)=table_mat(3);
        
    end
end





