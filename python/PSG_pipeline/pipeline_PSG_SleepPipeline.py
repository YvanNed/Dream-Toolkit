#!/usr/bin/env python
# coding: utf-8

# Start-up instructions
# (1) Install Conda
# Download at https://www.anaconda.com/
# Then everything is on the terminal

# (2) Install MNE with Conda
# conda update --name=base conda
# conda --version
# conda create --channel=conda-forge --strict-channel-priority --name=mne mne
# conda activate mne

# (3) Install YASA with Conda
# conda config --add channels conda-forge
# conda config --set channel_priority strict
# conda install yasa

# (4) Install Spyder with Conda
# conda install spyder

# In[1]:    
# importing packages
import os
import numpy as np
import matplotlib.pyplot as plt
import mne
import yasa
import pandas as pd
import glob
import seaborn as sns
import sklearn

# set up paths
if 'thandrillon' in os.getcwd():
    path_data='/Users/thandrillon/Data/ICEBERG/'
else:
    path_data='your_path'

if os.path.exists(path_data+"reports")==False:
    os.makedirs(path_data+"reports")
    
files = glob.glob(path_data + 'PSG' + os.sep + '**' + os.sep + '*.edf')
redo=0

# In[2]: Loop on the EDF
for file in files:
    
    # In[2.1] LOAD RAW DATA
    file_name=file.split(os.sep)[-1]
    path_info=file.split(os.sep)
    group_info=path_info[-2].split(' ')[-1]
    sub_ID=file_name.split('.')[0]
    report_prefix=path_data+"reports"+os.sep + sub_ID + '_' + group_info + '_'
    
    if os.path.exists(report_prefix+"report.html") and redo==0:
        continue

    # Define the channels you want to load
    selected_channels = ["Fp1", "C3", "O1", "A2"]
    raw = mne.io.read_raw_edf(file, preload=True, encoding="latin-1", include=selected_channels)

    # Exploring metadata
    n_time_samps = raw.n_times
    time_secs = raw.times
    ch_names = raw.ch_names
    n_chan = len(ch_names)  # note: there is no raw.n_channels attribute
    sf = raw.info['sfreq']
    
    if np.size(raw.ch_names)==0:
        continue
    
    print(raw)
    print('the data object has {} time samples and {} channels.'
          ''.format(n_time_samps, n_chan))
    print('The last time sample is at {} seconds.'.format(time_secs[-1]))
    print('The first few channel names are {}.'.format(', '.join(ch_names[:3])))
    print()  # insert a blank line in the output
    print('The sampling rate is {} Hz.'.format(raw.info['sfreq']))
    
    # adding to report
    report = mne.Report(title=file_name)
    report.add_raw(raw=raw, title="Raw: "+file_name, psd=False) 

    # In[2.2]: Plot a spectrogram
    # Load expert hypnogramn
    if group_info=='iRBD':
        hypno_name=path_data+os.sep+'Hypno'+os.sep+'HYP '+group_info+os.sep+sub_ID+"RBD.TXT"
    elif group_info=='PD':
        hypno_name=path_data+os.sep+'Hypno'+os.sep+'HYP '+group_info+os.sep+sub_ID+"P.TXT"
    elif group_info=='VS':
        hypno_name=path_data+os.sep+'Hypno'+os.sep+'HYP '+group_info+os.sep+sub_ID+"S.TXT"
        
    if os.path.exists(hypno_name):
        expert_hypno = np.loadtxt(hypno_name, dtype=str)
    else:
        continue
    
    if np.floor(raw.n_times/sf/30)!=np.size(expert_hypno):
        continue

        
    if np.size(np.where(np.isin(expert_hypno, ['N1','N2','N3'])))==0:
        expert_hypno = np.char.replace(expert_hypno, '1', 'N1')
        expert_hypno = np.char.replace(expert_hypno, '2', 'N2')
        expert_hypno = np.char.replace(expert_hypno, '3', 'N3')
    hypno_vec = np.full(np.size(expert_hypno), np.nan)  
    hypno_vec[np.where(np.isin(expert_hypno, 'W'))[0]]  = 0
    hypno_vec[np.where(np.isin(expert_hypno, 'N1'))[0]] = 1
    hypno_vec[np.where(np.isin(expert_hypno, 'N2'))[0]] = 2
    hypno_vec[np.where(np.isin(expert_hypno, 'N3'))[0]] = 3
    hypno_vec[np.where(np.isin(expert_hypno, 'R'))[0]]  = 4
    
    # C3
    channel_name = 'C3'  # Load channel C3
    data_c3, times = raw[channel_name, :]
    data_c3 = data_c3[0,:]
    
    ref_name = 'A2'  # Load channel A2
    data_a2, times = raw[ref_name, :]
    data_a2 = data_a2[0,:]

    data_c3_reref = data_c3 - data_a2 # Rereference C3 with A2
    hypno_fullvec = yasa.hypno_upsample_to_data(hypno=hypno_vec, sf_hypno=(1/30), data=data_c3_reref, sf_data=sf)
    fig = yasa.plot_spectrogram(data_c3_reref, sf, hypno_fullvec, fmin=0.5, fmax=40)
    
    report.add_figure(
        fig=fig,
        title="Hypnogram and Spectrogram on "+channel_name+"-"+ref_name,
        image_format="PNG",
        )
    del data_c3 
    del data_c3_reref
    
    # O1
    sf = raw.info['sfreq']
    channel_name = 'O1'  # Load channel C3
    data_o1, times = raw[channel_name, :]
    data_o1 = data_o1[0,:]
    data_o1_reref = data_o1 - data_a2 # Rereference C3 with A2
    hypno_fullvec = yasa.hypno_upsample_to_data(hypno=hypno_vec, sf_hypno=(1/30), data=data_o1_reref, sf_data=sf)
    fig = yasa.plot_spectrogram(data_o1_reref, sf, hypno_fullvec, fmin=0.5, fmax=40)
    
    report.add_figure(
        fig=fig,
        title="Hypnogram and Spectrogram on "+channel_name+"-"+ref_name,
        image_format="PNG",
        )
    
    del data_o1 
    del data_a2
    del data_o1_reref
    plt.close('all')

    # In[2.3]: Detect spindles
    channel_names = ['Fp1','C3','O1','A2'] # Select a few EEG channels
    sub_raw = raw.copy().pick(channel_names)
    sub_raw.load_data()
    sub_raw.set_eeg_reference(ref_channels=['A2']) # Rereference to A2 but could be A1 and A2
    sub_raw.drop_channels(['A2'])
        
    sub_raw = sub_raw.copy().filter(l_freq=0.4, h_freq=30)
    
    sp_det = yasa.spindles_detect( # Defauts parameters can be changed
            data = sub_raw,
            freq_sp = (12,15),
            freq_broad = (1,30),
            duration = (0.5, 2.5),
            min_distance = 500,
            thresh = {'corr': 0.65, 'rel_pow': None, 'rms': 1.5}, #Can be modified
            multi_only = False,
            remove_outliers = False,
            verbose = False,
            hypno=hypno_fullvec, 
            include=(2, 3)
            )
    if sp_det is not None:
        sp_det.summary(grp_chan=True)
        spindle_table = sp_det.summary()

        sp_det.plot_average(center='Peak', time_before=1, time_after=1);
        fig3=plt.gcf();
        report.add_figure(
            fig=fig3,
            title="Average spinldes",
            image_format="PNG",
            )
    
        sp_summary = sp_det.summary(grp_chan=True, grp_stage=True, aggfunc='mean')
        plt.close(fig3)
        
        spindle_table.to_csv(report_prefix+"automated_spindle_detection.csv")
        sp_summary.to_csv(report_prefix+"summary_spindle_detection.csv")
        
    
    # In[2.3]: Detect Slow Waves
    channel_names = ['Fp1','C3','O1'] # Select a few EEG channels
    sw = yasa.sw_detect(sub_raw, sf, ch_names=channel_names, hypno=hypno_fullvec, include=(2, 3))
    
    if sw is not None:
        sw_table = sw.summary()
        sw_summary = sw.summary(grp_chan=True, grp_stage=True, aggfunc='mean')
    
        # First, we need to get a long-format dataframe of peak-locked event.
        df_sync = sw.get_sync_events(center="NegPeak", time_before=0.4, time_after=0.8)
        # Then we can use seaborn.lineplot to create the desired plot:
        sns.lineplot(data=df_sync, x="Time", y="Amplitude", hue="Channel", 
                      style="Stage", palette="Set1", dashes=True)
    
        sns.despine()
        fig4=plt.gcf();
        report.add_figure(
            fig=fig4,
            title="Average Slow Waves per electrode and stage",
            image_format="PNG",
            )
        plt.close(fig4)
    
        sw_table.to_csv(report_prefix+"automated_sw_detection.csv")
        sw_summary.to_csv(report_prefix+"summary_sw_detection.csv")
    
    # In[2.4]: Automated sleep staging
    # From Vallat & Walker eLife 2021
    # For each PSG night, we extracted a single central EEG, left EOG, and chin EMG.
    # We chose a central EEG (e.g., C4-M1 or C4-Fpz depending on the dataset) since the American Academy of Sleep Medicine (AASM)
    # recommends that a central EEG should be included in all PSG recordings, and it is therefore more likely to be present
    # in a variety of PSG recordings.
    # These signals were then downsampled to 100 Hz to speed up computation time,
    # and bandpass-filtered between 0.40 Hz and 30 Hz.
    # No artifact removal was applied to the PSG data before running the sleep-staging algorithm.
    channel_names = ['C3', 'Fp1','A2']
    sub_raw = raw.copy().pick(channel_names)
    sub_raw.set_channel_types(dict(zip(channel_names, ['eeg','eeg','eeg'])))
    
    sub_raw.load_data()
    psg_reref = sub_raw.copy().resample(
        sfreq=100)
    
    eeg_channels = mne.pick_types(psg_reref.info, eeg=True)
    
    psg_reref.set_eeg_reference(ref_channels=['A2'])
    psg_reref_bp = psg_reref.copy().filter(l_freq=0.4, h_freq=30, picks=eeg_channels)
    
    sls = yasa.SleepStaging(psg_reref_bp, eeg_name="C3", eog_name="Fp1")
    pred_hypno = sls.predict()
    sls.plot_predict_proba();
    fig2=plt.gcf();
    

    
    report.add_figure(
        fig=fig2,
        title="Hypnodensity",
        image_format="PNG",
        )
    
    hypnodensities = sls.predict_proba()
    confidence = sls.predict_proba().max(1)
    df_pred = pd.DataFrame({'Expert': expert_hypno, 'Predicted': pred_hypno, 'Confidence': confidence})
    df_pred.reset_index(drop=True, inplace=True)
    hypnodensities.reset_index(drop=True, inplace=True)
    table_scoring = pd.concat([df_pred,hypnodensities], axis=1)   # FIX! NA VALUES  
    
    # compare with expert hypnogram
    accuracy = (expert_hypno == pred_hypno).sum() / pred_hypno.size
    print("The overall agreement is %.3f" % accuracy)
    # confusion matrix
    labels_hypno=['W','N1','N2','N3','R']
    cm = sklearn.metrics.confusion_matrix(expert_hypno, pred_hypno,labels=labels_hypno, normalize='true')
    bac = sklearn.metrics.balanced_accuracy_score(expert_hypno, pred_hypno)
        
    # Plot confusion matrix
    fig, ax = plt.subplots(1, layout="constrained")
    im = ax.imshow(cm, interpolation="nearest", cmap=plt.cm.Blues)
    ax.set(title="Norm. Confusion matrix")
    fig.colorbar(im)
    tick_marks = np.arange(len(labels_hypno))
    plt.xticks(tick_marks, labels_hypno, rotation=45)
    plt.yticks(tick_marks, labels_hypno)
    ax.set(ylabel="True label", xlabel="Predicted label")
    
    report.add_figure(
        fig=fig,
        title="Confusion Matrix",
        image_format="PNG",
        caption="Acc: %.3f | BAcc: %.3f" % (accuracy, bac),
       )
    
    table_scoring.to_csv(report_prefix+"automated_scoring.csv")

    # In[2.4]: Saving the report and other files
    report.save(report_prefix+"report.html", overwrite=True, open_browser=False)
    plt.close('all')
