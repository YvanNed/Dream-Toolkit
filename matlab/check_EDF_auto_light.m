

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

summary_table=array2table(zeros(2000,7),'VariableNames',{'File','Channel','Unit','Min','Max','BinGap','SamplingRate'});   % 2000 is an arbitrary number for preallocation
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
        summary_table.SamplingRate(nc) = hdr.orig.SampleRate(hdr.orig.chansel(i));
        
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


