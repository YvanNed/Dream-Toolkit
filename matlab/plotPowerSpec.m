
function plotPowerSpec(subfolder,Sub,data_pow,max_signal,epoch_score,freqs,all_channels)

nc1=0;
chan_labels={'Fp1','C3','Cz'};
nc1=nc1+1; 

% Import data and extract power
for nCh=1:length(chan_labels)

    temp_pow_W=data_pow(intersect(match_str(string(epoch_score),{'W'}),find(max_signal(:,match_str(all_channels,chan_labels{nCh}))<250)),match_str(all_channels,chan_labels{nCh}),:);
    Pow_Spec_W(nc1,nCh,:)=squeeze(nanmean(temp_pow_W,1));

    temp_pow_N1=data_pow(intersect(match_str(string(epoch_score),{'1'}),find(max_signal(:,match_str(all_channels,chan_labels{nCh}))<250)),match_str(all_channels,chan_labels{nCh}),:);
    Pow_Spec_N1(nc1,nCh,:)=squeeze(nanmean(temp_pow_N1,1));

    temp_pow_N2=data_pow(intersect(match_str(string(epoch_score),{'2'}),find(max_signal(:,match_str(all_channels,chan_labels{nCh}))<250)),match_str(all_channels,chan_labels{nCh}),:);
    Pow_Spec_N2(nc1,nCh,:)=squeeze(nanmean(temp_pow_N2,1));

    temp_pow_N3=data_pow(intersect(match_str(string(epoch_score),{'3','4'}),find(max_signal(:,match_str(all_channels,chan_labels{nCh}))<250)),match_str(all_channels,chan_labels{nCh}),:);
    Pow_Spec_N3(nc1,nCh,:)=squeeze(nanmean(temp_pow_N3,1));

    temp_pow_R=data_pow(intersect(match_str(string(epoch_score),{'R'}),find(max_signal(:,match_str(all_channels,chan_labels{nCh}))<250)),match_str(all_channels,chan_labels{nCh}),:);
    Pow_Spec_R(nc1,nCh,:)=squeeze(nanmean(temp_pow_R,1));

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

fpath = [subfolder filesep 'PowSpec'];
if exist(fpath,'file')==0
    mkdir(fpath)
end
saveas(gca,fullfile(fpath,sprintf('PowerSpectrum_%s',Sub)),'jpg')
save([fpath filesep sprintf('Power_%s',Sub)],'data_pow','freqs','Pow_Spec_W','Pow_Spec_N1','Pow_Spec_N2','Pow_Spec_N3','Pow_Spec_R');

end