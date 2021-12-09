

clear;
close all;

%% Set paths

path = '/Users/nico/Documents/ICM/Iceberg/Data/healthy';
path_data = [path filesep 'Pow_f_ft'];
files=dir([path_data filesep '*.mat']);
ft_defaults;

%% loop on subjects: retrieve and squeeze power for each sleep stage

nc1=0;

for nF=1:length(files)
% for nF=2
   
    % Parameters
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    separator = strfind(SubID,'_');
    Sub = SubID(1:separator(1)-1);
    fprintf('... working on Subject %s (%g/%g)\n',string(Sub),nF,length(files))        
    
    chan_labels={'Fp1','C3','Cz'};
    nc1=nc1+1; 
    
    % Import data and extract power
    load([path_data filesep SubID])
    
    for nCh=1:length(chan_labels)
        
        temp_pow_W=data_pow(intersect(match_str(string(epoch_score),{'W'}),find(max_signal(:,match_str(data_labels,chan_labels{nCh}))<250)),match_str(data_labels,chan_labels{nCh}),:);
        Pow_Spec_W(nc1,nCh,:)=squeeze(nanmean(temp_pow_W,1));
        
        temp_pow_N1=data_pow(intersect(match_str(string(epoch_score),{'1'}),find(max_signal(:,match_str(data_labels,chan_labels{nCh}))<250)),match_str(data_labels,chan_labels{nCh}),:);
        Pow_Spec_N1(nc1,nCh,:)=squeeze(nanmean(temp_pow_N1,1));
        
        temp_pow_N2=data_pow(intersect(match_str(string(epoch_score),{'2'}),find(max_signal(:,match_str(data_labels,chan_labels{nCh}))<250)),match_str(data_labels,chan_labels{nCh}),:);
        Pow_Spec_N2(nc1,nCh,:)=squeeze(nanmean(temp_pow_N2,1));
        
        temp_pow_N3=data_pow(intersect(match_str(string(epoch_score),{'3','4'}),find(max_signal(:,match_str(data_labels,chan_labels{nCh}))<250)),match_str(data_labels,chan_labels{nCh}),:);
        Pow_Spec_N3(nc1,nCh,:)=squeeze(nanmean(temp_pow_N3,1));

        temp_pow_R=data_pow(intersect(match_str(string(epoch_score),{'R'}),find(max_signal(:,match_str(data_labels,chan_labels{nCh}))<250)),match_str(data_labels,chan_labels{nCh}),:);
        Pow_Spec_R(nc1,nCh,:)=squeeze(nanmean(temp_pow_R,1));
        
    end
    
%     figure('Position',[1,83,1440,714]);
% 
%     for nCh=1:3
% 
%         subplot(1,3,nCh);
%         plot(freqs,squeeze(nanmean(Pow_Spec_W(nF,nCh,:),1))','LineWidth',2)
%         hold on; plot(freqs,squeeze(nanmean(Pow_Spec_N1(nF,nCh,:),1))','LineWidth',2)
%         hold on; plot(freqs,squeeze(nanmean(Pow_Spec_N2(nF,nCh,:),1))','LineWidth',2)
%         hold on; plot(freqs,squeeze(nanmean(Pow_Spec_N3(nF,nCh,:),1))','LineWidth',2)
%         hold on; plot(freqs,squeeze(nanmean(Pow_Spec_R(nF,nCh,:),1))','LineWidth',2)
% 
%         title([chan_labels{nCh}]);
%         xlabel('Freq (Hz)')
%         ylabel('Log Power')
%         set(gcf,'Color','w')
%         set(gca,'FontSize',18,'FontWeight','bold')
%         legend({'Wake','N1','N2','N3','REM'});
% 
%     end
% 
%     S = sgtitle(sprintf('Sub %s',Sub));
%     S.FontSize = 24;
%     S.FontWeight = 'bold';
% 
%     fpath = [folder_name filesep 'Plots'];
%     % saveas(gca,fullfile(fpath,sprintf('%s',Sub)),'jpg')

end

%% Plot power across sleep stages, on 3 channels

figure('Position',[1 377 1084 420]);

for nCh=1:3
    
    subplot(1,3,nCh);
    plot(freqs,squeeze(nanmean(Pow_Spec_W(:,nCh,:),1))','LineWidth',2)
    hold on; plot(freqs,squeeze(nanmean(Pow_Spec_N1(:,nCh,:),1))','LineWidth',2)
    hold on; plot(freqs,squeeze(nanmean(Pow_Spec_N2(:,nCh,:),1))','LineWidth',2)
    hold on; plot(freqs,squeeze(nanmean(Pow_Spec_N3(:,nCh,:),1))','LineWidth',2)
    hold on; plot(freqs,squeeze(nanmean(Pow_Spec_R(:,nCh,:),1))','LineWidth',2)

    title([chan_labels{nCh}]);
    xlabel('Freq (Hz)')
    ylabel('Log Power')
    set(gcf,'Color','w')
    set(gca,'FontSize',18,'FontWeight','bold')
    legend({'Wake','N1','N2','N3','REM'});

end

%% Save

fpath = [folder_name filesep 'Plots'];
% saveas(gca,fullfile(fpath,'mean_sub'),'jpg')

