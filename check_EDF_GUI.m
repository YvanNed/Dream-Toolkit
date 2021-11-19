

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

%% Select EDFs to check

% Return the subject IDs from the data folder
filelist = dir([subfolder filesep '**' filesep '*.edf']);
pick=listdlg('ListString',{filelist.name},'PromptString','Select the EDF file to check');
filelist = filelist(pick);
fprintf('>>> You have selected %g EDF files\n',length(filelist))

%% Loop across subjects

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
    pick_channels = listdlg('ListString',hdr.label,'PromptString','Select the channels to check');
    all_channels  = hdr.label(pick_channels);
    Num_ch = numel(all_channels);
     
%%%  Visualise the data
%    cfg             = [];
%    cfg.dataset     = [rootdir filesep subID];
%    cfg.channel      = cellstr(Channels);
%    Preproc_data    = ft_preprocessing(cfg);
%    cfg.blocksize   = 30; % in sec
%    cfg.channel     = cellstr(Channels); 
%    cfg.viewmode    = 'vertical';
%    ft_databrowser(cfg, Preproc_data);
    
    %% Check for signal clipping and bit depth issue
        
    fprintf(1,'>>> >>> >>> Plotting histograms for each of the %s channels...\n',string(Num_ch))

    for i = 1:Num_ch
        
        Data = data(pick_channels(i),:);

        f=figure;

        % Check for signal clipping issue
        subplot(1,2,1); ax=gca;

        H = histogram(Data,-max(abs(Data)):0.05:max(abs(Data)),'EdgeColor','#1167b1'); 
        xlim([-1.05 1.05]*max(abs(Data)))
        ylim([0 1]*max(H.Values(H.BinEdges(1:end-1)<=-0.05 | H.BinEdges(1:end-1)>=0.05)))
        t = title('Amplitude distribution');
        t.FontWeight = 'normal';
        xlabel('Amplitude (μV)'); ylabel('Data points distribution')
        ax.FontSize = 14;

        % Check for bit depth issue 
        subplot(1,2,2); ax=gca;

        % Plot the absolute difference in amplitude between neihboring data points 
        delta_ampl = abs(diff(Data));

        H2=histogram(delta_ampl,0:0.01:50,'EdgeColor','#1167b1');
        xlabel('Delta amplitude (μV)'); ylabel('Data points distribution') 
        t = title({'Absolute difference in amplitude between';'neighboring data points'});
        t.FontWeight = 'normal';
        ax.FontSize = 14;
        
        T = sgtitle({sprintf('Subject %s',Sub);sprintf('Channel %s',all_channels{i})}); 
        T.FontWeight = 'bold';
        ylim([0 1]*max(H2.Values(H2.BinEdges(1:end-1)>=0.01)))
        f.Position = [-96,1387,906,420];
        set(gcf,'color','w');

    end

    %% Check for polarity issue
    fprintf(1,'>>> >>> >>> Checking potential polarity issues...\n')

    % Plot the correlation matrix between the selected channels.
    [r, pV] = corr(data(pick_channels,:)');
    
    g=figure; ax=gca;
    
    imagesc(r); 
    colorbar
    caxis([-1 1])
    xticks(1:Num_ch)
    yticks(1:Num_ch)
    xticklabels(all_channels)
    yticklabels(all_channels)
    ax.XAxis.TickLength = [0 0];
    ax.YAxis.TickLength = [0 0];
    xtickangle(ax,45)
    t = title({sprintf('Subject %s',Sub);'Correlation between all channels'});
    t.FontWeight = 'normal';
    
    % Add text for negative correlation (r < -0.2)
    t = cell(Num_ch, Num_ch);
 
    for i=1:Num_ch
        for j=1:Num_ch
            t(i, j) = cellstr(num2str(round(r(i,j), 2)));
        end
    end
    
    t = cellfun(@str2num,t);
    [x,y] = find(t<=-0.2);
    for i = 1:numel(x)
        text(x(i), y(i), string(t(x(i),y(i))), 'HorizontalAlignment', 'Center', 'FontSize', 12);
    end

    ax.FontSize = 14;
    g.Position = [811,1185,773,623];
    set(gcf,'color','w');
    
    if length(filelist)>1
        fprintf(1,'Press SPACE to continue to next EDF file.\n')
        pause;
        close all;
    end
    
end





