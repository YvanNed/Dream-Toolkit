

%% EDF format checks

% There are 3 potential issues related to the conversion of raw signals to
% EDF formats:

% 1. Signal clipping: signal cut once it exceeds an amplitude threshold
% (the min-max range set before conversion was too narrow)

% 2. Bit depth: signal shows a stair-like progression (the min-max
% range set before conversion was too large)

% 3. Inverted polarity: signal multiplied by -1 

%% Load the data

clear;
close all;

rootdir = '/Users/nico/Documents/ICM/Iceberg/Data/healthy';
addpath '/Users/nico/Documents/ICM/Edison/Scripts'
run localdef.m
addpath((path_fieldtrip));
addpath((path_fieldtrip_adv));
addpath((path_fieldtrip_adv2));

% Loop across subjects
Subs_healthy = {'73','75','91','94','98','99','104','106','107','109'};
Subs_parks = {'79','95','165','185','193','222','246','269','282','284'};

for S = 1:2

    % Parameters subject
    subID = cell2mat(strcat(Subs_healthy(S),'.edf'));
    separator = strfind(subID,'.edf'); 
    Sub = subID(1:separator(1)-1);

    % Import the data
    cfg = [];
    cfg.dataset = [rootdir filesep subID];
    data = ft_read_data(cfg.dataset);
    hdr = ft_read_header(cfg.dataset);

    % A messy code to return the correct EEG channel indices
    all_channels = string(hdr.label);
    eeg_channels_orig_idx = find(contains(string(hdr.orig.Transducer),'EEG'));
    
    for y = 1:numel(eeg_channels_orig_idx)
        full_label = hdr.orig.Transducer(eeg_channels_orig_idx(y),:);
        separator2 = strfind(full_label,' '); 
        short_label = full_label(separator2(1)+1:separator2(2)-1);
        
        EEG_channels_idx(y) = find(contains(string(hdr.label),short_label));
    end
    
    EEG_channels = all_channels(EEG_channels_idx);
    
%     Visualise the data
%     cfg = [];
%     Preproc_data = ft_preprocessing(cfg);
%     cfg.blocksize = 30; % in sec
%     cfg.channel = cellstr(EEG_channels); 
%     cfg.viewmode = 'vertical';
%     ft_databrowser(cfg, Preproc_data);
%     
    
    %% Check for signal clipping and bit depth issue
    
    % For each channel, check for signal clipping / bith depth issue
    for i = 1:numel(EEG_channels)

        Data = data(EEG_channels_idx(i),:);

        f=figure;

        % Signal clipping: outstanding values (ie, min/max) in the histogram?
        subplot(1,2,1); ax=gca;

        histogram(Data,-max(abs(Data)):0.05:max(abs(Data)),'EdgeColor','#1167b1') 
        xlim([-1.05 1.05]*max(abs(Data)))
        t = title('Data points distribution');
        t.FontWeight = 'normal';

        xlabel('Amplitude'); ylabel('Distribution')
        ax.FontSize = 14;

        % Plot the absolute difference in amplitude between neihboring data points 
        % --> Signal clipping if two values peak out
        % --> Bit depth issue if gaps observed between evenly distributed values
        subplot(1,2,2); ax=gca;

        delta_ampl = abs(diff(Data));

        histogram(delta_ampl,0:0.01:50,'EdgeColor','#1167b1')
        t = title({'Absolute difference in amplitude between';'neighboring data points'});
        t.FontWeight = 'normal';

        xlabel('Delta amplitude'); ylabel('Distribution') 
        ax.FontSize = 14;
        f.Position = [459,1143,906,420];

        % Artificial bit depth issue: stair-like evolution
        % figure; histogram(abs(diff(round(data(EEG_channel,:)/10)*10)),0:0.01:50)

        T = sgtitle({sprintf('Subject %s',Sub);sprintf('Channel %s',EEG_channels(i))}); 
        T.FontWeight = 'bold';

        f.Position = [-60,1399,906,420];

    end

    %% Check for polarity issue

    % Compute the correlation between the 3 EEG channels. Positive correlation
    % means that all channels have the same polarity (but doesn't rule out the
    % possibilty that all have an inverse polarity——needs to be manually
    % checked by comparing with raw signals on Compumedics)
    [r, pV] = corr(data(EEG_channels_idx,:)');

    % Plot the correlation matrix
    g=figure; ax=gca;
    imagesc(r); 

    colorbar
    caxis([-1 1])
    xticklabels(EEG_channels)
    yticklabels(EEG_channels)
    ax.XAxis.TickLength = [0 0];
    ax.YAxis.TickLength = [0 0];
    t = title('Correlation between EEG channels');
    t.FontWeight = 'normal';
    ax.FontSize = 14;

    g.Position = [861,1398,560,420];

end



