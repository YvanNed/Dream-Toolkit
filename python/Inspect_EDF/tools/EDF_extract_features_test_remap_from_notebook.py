#!/usr/bin/env python
# coding: utf-8

# In[1]:
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mne
import yasa
import pandas as pd
import glob
import seaborn as sns
from scipy.stats import sem
import copy
import json
from pathlib import Path 


matplotlib.use('Agg') #Qt5Agg #Agg

# Définition du chemin des données
if 'thandrillon' in os.getcwd():
    path_data = '/Users/thandrillon/Data/Apomorphee'
elif 'yvan.nedelec' in os.getcwd():
    path_data = '/Users/yvan.nedelec/OneDrive - ICM/Documents/RE_yvan/projects/noemie_rescue'
else:
    path_data = '/Users/noemiewane/Library/Mobile Documents/com~apple~CloudDocs/Internat/Master 2/Stage/Données/Version 4/'

# Création du dossier "reports"
reports_dir = os.path.join(path_data, "reports_251006")
os.makedirs(reports_dir, exist_ok=True)

# Création du dossier "plots"
plots_path = os.path.join(path_data, "plots_251006")
os.makedirs(plots_path, exist_ok=True)

# Recherche des fichiers EDF
files = glob.glob(os.path.join(path_data+os.sep+'data', '*.edf'))
redo = 1

# load json file for remap channels and re-ref
config_dict = json.loads(Path("C:/Users/yvan.nedelec/OneDrive - ICM/Documents/RE_yvan/projects/noemie_rescue/data/config_param/remap_reref_persubject_test_voila1.json").read_text())

# In[2]:    
# Boucle sur les fichiers EDF
all_sp_waveform_av = pd.DataFrame()
all_sw_waveform_av = pd.DataFrame()
all_sw = pd.DataFrame()
all_sw_av = pd.DataFrame()
all_sp = pd.DataFrame()
all_sp_av = pd.DataFrame()

all_low_sw = pd.DataFrame()
all_low_sw_av = pd.DataFrame()

all_high_sw = pd.DataFrame()
all_high_sw_av = pd.DataFrame()

all_proba_hypno_av = pd.DataFrame()

psds_all_elec_N1_APO = []
psds_all_elec_N2_APO = []
psds_all_elec_N1_PLA = []
psds_all_elec_N2_PLA = []
# open table about subject identity and condition
subj_table = pd.read_csv(os.path.join(path_data, 'CondTable.csv'), sep=';')

# list to catch the files encountering a pb (erreur de lecture des données eeg, de lecture de l'hypnogramme, de mismatch taille de l'hypno et taille eeg...)
failed_files = []

#______________________________________________________________________________
# to select a specific file based on subject label (debug or test feature)_____
# sub_file = "10_N1"
# # create a pattern as a regex (regular expression) object that allows operation such as search, match etc. 
# # we add (?<!\d) before the subject name to avoid subject names with digits before the subject name of interest
# pattern = re.compile(rf"(?<!\d){re.escape(sub_file)}") # the re.escape allows to avoid errors if there are special characters (e.g. .,+,-,(,) etc) in the string of interest (you have to add a \ to avoid the function of the special character). if you write "(pdf)", re.escape("(pdf)") will return "\(pdf\)"
# # search for the pattern 
# mask_file = np.array([bool(pattern.search(f)) for f in files])
# file = files[np.where(mask_file)[0][0]]

# else you can loop through files to identify which one have the targeted name
# for t, test in enumerate(files):
#     if os.path.basename(test).split(".")[0] == sub_file:
#         file = files[t]
#______________________________________________________________________________

all_ERP_SW = []
report = mne.Report(title='Micro-Structure all subjects')
for file in files:
    file_name = os.path.basename(file)  # Ex: "15_N1.edf"
    file_ID = file_name.split('.')[0]  # Extrait "15_N1"
    sub_ID = file_ID.split("_")[0]
    night_ID = file_ID.split("_")[-1]
    drug_ID = subj_table.loc[subj_table['SubID'] == int(sub_ID), night_ID].values[0]
    group_ID = 'G_N1_'+subj_table.loc[subj_table['SubID'] == int(sub_ID), 'N1'].values[0]
    report_prefix = os.path.join(reports_dir, f"{file_ID}_report")
    
    sub_config = config_dict[file_ID]

    header_html = f"<h2>Subject: {sub_ID} — Condition: {drug_ID} — Session: {night_ID}</h2>"
    report.add_html(
        html=header_html,
        title=f"{sub_ID}_{drug_ID}_{night_ID}",
        section=None  # or use a specific section name if grouping
    )
    
    selected_channels = list(sub_config["remap"].keys())
    
    try:
        raw = mne.io.read_raw_edf(file, preload=True, encoding="latin-1", include=selected_channels)
    except Exception as e:
        print(f"⚠️ Erreur lecture fichier {file_name} : {e}")
        failed_files.append((file, 'eeg loading'))
        continue
    
    # remap channels (based on the generated dict from jupy notebook):
    raw.rename_channels(sub_config["remap"])
    
    # set eeg reference
    if  sub_config["ref_channels"] != 'average':
        raw.set_eeg_reference(ref_channels = sub_config["ref_channels"])
        # get rid of the ref
        raw.drop_channels(sub_config["ref_channels"])
    else:
        raw.set_eeg_reference(ref_channels = sub_config["ref_channels"])
    
    # get sf
    sf = raw.info["sfreq"]
    
    print(f"📊 {file_name}: {raw.n_times} échantillons, {len(raw.ch_names)} canaux à {sf} Hz.")
    
    # Chargement de l'hypnogramme expert
    hypno_name = os.path.join(path_data+os.sep+'data', file_ID+"_Hypnogram_remapped.txt")
    if not os.path.exists(hypno_name):
        print(f"⚠️ Hypnogramme introuvable pour {sub_ID} : {hypno_name}")
        failed_files.append((file, 'no hypno'))
        continue

    try:
        expert_hypno = np.loadtxt(hypno_name, dtype=str)
        expert_hypno = expert_hypno.astype('<U10')
    except Exception as e:
        print(f"⚠️ Erreur chargement hypnogramme pour {sub_ID} : {e}")
        failed_files.append((file, 'hypno loading'))
        continue
    
    
    # Vérification de la compatibilité
    expected_epochs = int(np.floor(raw.n_times / sf / 30))
    if expected_epochs != len(expert_hypno):
        print(f"⚠️ Incompatibilité hypnogramme ({len(expert_hypno)}) vs EEG ({expected_epochs} epochs)")
        failed_files.append((file, 'mismatch hypno - eeg'))
        continue

    # Normalisation des labels
    if np.any(np.isin(expert_hypno, "?")):
        expert_hypno = np.strings.replace(expert_hypno, '?', 'W').astype('<U10')
    
    if not np.any(np.isin(expert_hypno, ['N1', 'N2', 'N3'])):
        expert_hypno = np.strings.replace(expert_hypno, '1', 'N1')
        expert_hypno = np.strings.replace(expert_hypno, '2', 'N2')
        expert_hypno = np.strings.replace(expert_hypno, '3', 'N3')
        expert_hypno = np.strings.replace(expert_hypno, '4', 'N3')
    else:
        try:
            expert_hypno = np.loadtxt(hypno_name, dtype=str)
            print(f"✅ Hypnogramme chargé, nombre d'entrées : {np.size(expert_hypno)}")
        except Exception as e:
            print(f"⚠️ Erreur lors du chargement de l'hypnogramme : {e}")
            continue
    if expert_hypno.ndim>1:
        expert_hypno=expert_hypno[:,2]
        
    # Conversion en vecteur numérique
    hypno_vec = np.full(len(expert_hypno), np.nan)
    mapping = {'W': 0, 'N1': 1, 'N2': 3, 'N3': 3, 'R': 4}
    for stage, value in mapping.items():
        hypno_vec[np.where(expert_hypno == stage)] = value
        
        
    new_sf = 200
    raw.resample(new_sf, npad="auto")
    
    # --- Génération des figures hypnogramme/spectrogramme ---
    # Pour le canal C4-M1
    channel_name = new_channels[1]
    data_c4, times = raw[channel_name, :]
    data_c4 = data_c4[0, :]
    hypno_fullvec = yasa.hypno_upsample_to_data(hypno=hypno_vec,
                                                sf_hypno=(1/30),
                                                data=data_c4,
                                                sf_data=new_sf)    

 
    # --- Détection des spindles ---
    sub_raw = raw.copy().filter(l_freq=0.1, h_freq=30)
    sp_det = yasa.spindles_detect(
        data=sub_raw,
        freq_sp=(11, 16),
        hypno=hypno_fullvec,
        include=(3)
    )
    
    data_sp_avg = pd.DataFrame()
    if sp_det is not None:

        det_sp = sp_det.summary()
        
        data_sp=sp_det.get_sync_events(center='Peak', time_before=0.75, time_after=0.75)
        
        data_sp_avg = data_sp.groupby(['Time', 'Channel'], as_index=False)['Amplitude'].mean()
        data_sp_avg["SubID"]=sub_ID
        data_sp_avg["Night"]=night_ID
        data_sp_avg["Drug"]=drug_ID
        data_sp_avg["GroupID"]=group_ID
        
        fig3 = plt.figure(figsize=(10, 6))
        sns.lineplot(data=data_sp_avg, x='Time', y='Amplitude', hue='Channel', palette='tab10')
        
        plt.title('SP Waveforme Across Time for Each Electrode')
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude (µV)')
        plt.legend(title='Electrode')
        plt.tight_layout()
        
        report.add_figure(fig=fig3, title="Moyenne des fuseaux", image_format="PNG")
        plt.close(fig3)
        
        sp_det.summary().to_csv(report_prefix + "_spindle_detection_NREM.csv")

    else:
        det_sp = pd.DataFrame(columns=['Start', 'Peak', 'End', 'Duration', 'Amplitude', 'RMS', 'AbsPower',
               'RelPower', 'Frequency', 'Oscillations', 'Symmetry', 'Stage', 'Channel',
               'IdxChannel', 'SubID'])

    det_sp_av = det_sp.groupby(['Channel','Stage'], as_index=False).mean()
    det_sp_av = det_sp_av.drop(['Start', 'Peak', 'End','IdxChannel'], axis=1)
    row_count = det_sp.groupby(['Channel','Stage']).size().reset_index(name='count')
    det_sp_av = pd.merge(det_sp_av, row_count, on=['Channel','Stage'], how='left')
    unique, counts = np.unique(hypno_vec, return_counts=True)
    count_dict = dict(zip(unique.astype(int), counts))
    expected_keys = [3]    
    for key in expected_keys:
        if key not in count_dict:
            count_dict[key] = 0
    det_sp_av['density'] = np.nan

    # add NaN or 0 values if no sp detected in a certain stage
    channels = sub_raw.ch_names
    full_index = pd.MultiIndex.from_product([channels], names=['Channel'])
    # Step 2: Set index to ['Channel', 'Stage'] and reindex with full combinations
    det_sp_av['density'] = det_sp_av['density'].fillna(0)
    if det_sp_av.empty==False:
        det_sp_av['density']=det_sp_av['count']/(count_dict[3]/2)
        
        columns=det_sp_av.columns
        existing_channels = det_sp_av['Channel'].unique()
        # Find missing channels
        missing_channels = [ch for ch in channels if ch not in existing_channels]
        # Create rows for missing channels
        new_rows = []
        for ch in missing_channels:
            row = {
                'Channel': ch,
                'Stage': 3,
                'Duration': np.nan,
                'Amplitude': np.nan,
                'RMS': np.nan,
                'AbsPower': np.nan,
                'RelPower': np.nan,
                'Frequency': np.nan,
                'Oscillations': np.nan,
                'Symmetry': np.nan,
                'count': 0,
                'density': 0
            }
            new_rows.append(row)
        
        # Append to the DataFrame
        if new_rows:
            det_sp_av = pd.concat([det_sp_av, pd.DataFrame(new_rows)], ignore_index=True)
        
        # Ensure column order is preserved
        det_sp_av = det_sp_av[columns]
        
    else:
        columns=det_sp_av.columns
        for ch in channels:
            row = {
                'Channel': ch,
                'Stage': 3,
                'Duration': np.nan,
                'Amplitude': np.nan,
                'RMS': np.nan,
                'AbsPower': np.nan,
                'RelPower': np.nan,
                'Frequency': np.nan,
                'Oscillations': np.nan,
                'Symmetry': np.nan,
                'count': 0,
                'density': 0
            }
            det_sp_av = pd.concat([det_sp_av, pd.DataFrame([row])], ignore_index=True)
        det_sp_av = det_sp_av[columns]

    det_sp_av["SubID"]=sub_ID
    det_sp["SubID"]=sub_ID
    det_sp_av["Night"]=night_ID
    det_sp["Night"]=night_ID
    det_sp_av["Drug"]=drug_ID
    det_sp["Drug"]=drug_ID
    det_sp_av["GroupID"]=group_ID
    det_sp["GroupID"]=group_ID   
    
    # --- Détection des ondes lentes ---
    sw = yasa.sw_detect(sub_raw, sf, hypno=hypno_fullvec,include=(3),
                        freq_sw=(0.2, 4), dur_neg=(0.125, 2.5), dur_pos=(0.1, 2), amp_neg=(5, 200), amp_pos=(5, 150), amp_ptp=(10, 350), coupling=False, remove_outliers=False, verbose=False)
    
    
    det_sw = sw.summary()
    det_low_sw = det_sw[(det_sw['Frequency'] >= 0.2) & (det_sw['Frequency'] <= 1)]
    det_high_sw = det_sw[(det_sw['Frequency'] > 1) & (det_sw['Frequency'] <= 4)]

    # Step 1: Extract the negative peak times from det_sw
    events = (det_sw["NegPeak"] * raw.info["sfreq"]).astype(int)  # sample indices
    events = np.column_stack((events, np.zeros(len(events), dtype=int), np.ones(len(events), dtype=int)))
    
    
    # Parameters
    tmin, tmax = -1.0, 1.0
    baseline = (-1.0, -0.5)   # baseline window
    sfreq = raw.info["sfreq"]
    times = np.arange(int(tmin*sfreq), int(tmax*sfreq)) / sfreq
    b_start = int((baseline[0] - tmin) * sfreq)
    b_stop  = int((baseline[1] - tmin) * sfreq)
    
    # Define conditions
    conditions = {
        "All SW": det_sw,
        "0.2–1 Hz": det_sw[(det_sw["Frequency"] >= 0.2) & (det_sw["Frequency"] < 1)],
        "1–4 Hz": det_sw[(det_sw["Frequency"] >= 1) & (det_sw["Frequency"] <= 4)],
    }
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    
    for ax, (label, df) in zip(axes, conditions.items()):
        for ch in df["Channel"].unique():
            chan_events = df[df["Channel"] == ch]["NegPeak"].values
            
            epochs = []
            for ev in chan_events:
                ev_sample = int(ev * sfreq)
                start = ev_sample + int(tmin * sfreq)
                stop  = ev_sample + int(tmax * sfreq)
                if start >= 0 and stop < raw.n_times:
                    data, _ = raw[ch, start:stop]
                    # --- Baseline correction
                    base_mean = data[b_start:b_stop].mean()
                    data = data - base_mean
                    epochs.append(data[0])
                                
            if len(epochs) == 0:
                continue
            
            epochs = np.vstack(epochs)
            mean_wave = epochs.mean(axis=0)
            sem_wave = epochs.std(axis=0) / np.sqrt(len(epochs))
            
            ax.plot(times, mean_wave, label=ch)
            ax.fill_between(times, mean_wave - sem_wave, mean_wave + sem_wave, alpha=0.3)
            
            temp_mean_wave = pd.DataFrame({
                "SubID": sub_ID,
                "Drug": drug_ID,
                "Night": night_ID,
                "Channel": ch,
                "Condition": label,
                "Time": times,
                "Amplitude": mean_wave
            })
            all_ERP_SW.append(temp_mean_wave)
        
        ax.axvline(0, color="k", linestyle="--")
        ax.set_title(label)
        ax.set_xlabel("Time (s)")
    
    axes[0].set_ylabel("Amplitude (µV)")
    axes[-1].legend()
    fig.suptitle("Slow-wave ERPs by frequency band", fontsize=14)
    plt.tight_layout()
    plt.show()

    
    report.add_figure(fig=fig, title="Moyenne des ondes lentes", image_format="PNG")
    plt.close(fig)
    

    det_sw_av = det_sw.groupby(['Channel','Stage'], as_index=False).mean()
    det_sw_av = det_sw_av.drop(['Start', 'NegPeak', 'MidCrossing', 'PosPeak', 'End','IdxChannel'], axis=1)
    row_count = det_sw.groupby(['Channel','Stage']).size().reset_index(name='count')
    det_sw_av = pd.merge(det_sw_av, row_count, on=['Channel','Stage'], how='left')
    unique, counts = np.unique(hypno_vec, return_counts=True)
    count_dict = dict(zip(unique.astype(int), counts))
    expected_keys = [3]    
    for key in expected_keys:
        if key not in count_dict:
            count_dict[key] = 0
    det_sw_av['density'] = np.nan
        
    channels = sub_raw.ch_names
    full_index = pd.MultiIndex.from_product([channels], names=['Channel'])
    det_sw_av['density'] = det_sw_av['density'].fillna(0)
    if det_sw_av.empty==False:
        k=3
        det_sw_av['density']=det_sw_av['count']/(count_dict[k]/2)

    det_sw_av["SubID"]=sub_ID
    det_sw["SubID"]=sub_ID
    det_sw["GroupID"]=group_ID
    det_sw_av["GroupID"]=group_ID
    det_sw_av["Night"]=night_ID
    det_sw["Night"]=night_ID
    det_sw_av["Drug"]=drug_ID
    det_sw["Drug"]=drug_ID
   
    det_low_sw_av = det_low_sw.groupby(['Channel','Stage'], as_index=False).mean()
    det_low_sw_av = det_low_sw_av.drop(['Start', 'NegPeak', 'MidCrossing', 'PosPeak', 'End','IdxChannel'], axis=1)
    row_count = det_low_sw.groupby(['Channel','Stage']).size().reset_index(name='count')
    det_low_sw_av = pd.merge(det_low_sw_av, row_count, on=['Channel','Stage'], how='left')
    unique, counts = np.unique(hypno_vec, return_counts=True)
    count_dict = dict(zip(unique.astype(int), counts))
    expected_keys = [3]    
    for key in expected_keys:
        if key not in count_dict:
            count_dict[key] = 0
    det_low_sw_av['density'] = np.nan
        
    det_low_sw_av['density'] = det_low_sw_av['density'].fillna(0)
    if det_low_sw_av.empty==False:
        k=3
        det_low_sw_av['density']=det_low_sw_av['count']/(count_dict[k]/2)
        
        columns=det_low_sw_av.columns
        existing_channels = det_low_sw_av['Channel'].unique()
        # Find missing channels
        missing_channels = [ch for ch in channels if ch not in existing_channels]
        # Create rows for missing channels
        new_rows = []
        for ch in missing_channels:
            row = {
                'Channel': ch,
                'Stage': 3,
                'Duration': np.nan,
                'ValNegPeak': np.nan,
                'ValPosPeak': np.nan,
                'PTP': np.nan,
                'Slope': np.nan,
                'Frequency': np.nan,
                'count': 0,
                'density': 0
            }
            new_rows.append(row)
        
        # Append to the DataFrame
        if new_rows:
            det_low_sw_av = pd.concat([det_low_sw_av, pd.DataFrame(new_rows)], ignore_index=True)
        
        # Ensure column order is preserved
        det_low_sw_av = det_low_sw_av[columns]
        
    else:
        columns=det_low_sw_av.columns
        for ch in channels:
            row = {
                'Channel': ch,
                'Stage': 3,
                'Duration': np.nan,
                'ValNegPeak': np.nan,
                'ValPosPeak': np.nan,
                'PTP': np.nan,
                'Slope': np.nan,
                'Frequency': np.nan,
                'count': 0,
                'density': 0
            }
            det_low_sw_av = pd.concat([det_low_sw_av, pd.DataFrame([row])], ignore_index=True)
        det_low_sw_av = det_low_sw_av[columns]
        
    det_low_sw_av["SubID"]=sub_ID
    det_low_sw["SubID"]=sub_ID
    det_low_sw_av["GroupID"]=group_ID
    det_low_sw["GroupID"]=group_ID
    det_low_sw_av["Night"]=night_ID
    det_low_sw["Night"]=night_ID
    det_low_sw_av["Drug"]=drug_ID
    det_low_sw["Drug"]=drug_ID
      
    det_high_sw_av = det_high_sw.groupby(['Channel','Stage'], as_index=False).mean()
    det_high_sw_av = det_high_sw_av.drop(['Start', 'NegPeak', 'MidCrossing', 'PosPeak', 'End','IdxChannel'], axis=1)
    row_count = det_high_sw.groupby(['Channel','Stage']).size().reset_index(name='count')
    det_high_sw_av = pd.merge(det_high_sw_av, row_count, on=['Channel','Stage'], how='left')
    unique, counts = np.unique(hypno_vec, return_counts=True)
    count_dict = dict(zip(unique.astype(int), counts))
    expected_keys = [3]    
    for key in expected_keys:
        if key not in count_dict:
            count_dict[key] = 0
    det_high_sw_av['density'] = np.nan
        
    det_high_sw_av['density'] = det_high_sw_av['density'].fillna(0)
    if det_high_sw_av.empty==False:
        k=3
        det_high_sw_av['density']=det_high_sw_av['count']/(count_dict[k]/2)
        
        columns=det_high_sw_av.columns
        existing_channels = det_high_sw_av['Channel'].unique()
        # Find missing channels
        missing_channels = [ch for ch in channels if ch not in existing_channels]
        # Create rows for missing channels
        new_rows = []
        for ch in missing_channels:
            row = {
                'Channel': ch,
                'Stage': 3,
                'Duration': np.nan,
                'ValNegPeak': np.nan,
                'ValPosPeak': np.nan,
                'PTP': np.nan,
                'Slope': np.nan,
                'Frequency': np.nan,
                'count': 0,
                'density': 0
            }
            new_rows.append(row)
        
        # Append to the DataFrame
        if new_rows:
            det_high_sw_av = pd.concat([det_high_sw_av, pd.DataFrame(new_rows)], ignore_index=True)
        
        # Ensure column order is preserved
        det_high_sw_av = det_high_sw_av[columns]
    else:
        columns=det_high_sw_av.columns
        for ch in channels:
            row = {
                'Channel': ch,
                'Stage': 3,
                'Duration': np.nan,
                'ValNegPeak': np.nan,
                'ValPosPeak': np.nan,
                'PTP': np.nan,
                'Slope': np.nan,
                'Frequency': np.nan,
                'count': 0,
                'density': 0
            }
            det_high_sw_av = pd.concat([det_high_sw_av, pd.DataFrame([row])], ignore_index=True)
        det_high_sw_av = det_high_sw_av[columns]

    det_high_sw_av["SubID"]=sub_ID
    det_high_sw["SubID"]=sub_ID
    det_high_sw_av["GroupID"]=group_ID
    det_high_sw["GroupID"]=group_ID
    det_high_sw_av["Night"]=night_ID
    det_high_sw["Night"]=night_ID
    det_high_sw_av["Drug"]=drug_ID
    det_high_sw["Drug"]=drug_ID
    
    
    all_sw = pd.concat([all_sw, det_sw], ignore_index=True)
    all_sw_av = pd.concat([all_sw_av, det_sw_av], ignore_index=True)
    all_sp = pd.concat([all_sp, det_sp], ignore_index=True)
    all_sp_av = pd.concat([all_sp_av, det_sp_av], ignore_index=True)
    
    all_low_sw = pd.concat([all_low_sw, det_low_sw], ignore_index=True)
    all_low_sw_av = pd.concat([all_low_sw_av, det_low_sw_av], ignore_index=True)

    
    all_high_sw = pd.concat([all_high_sw, det_high_sw], ignore_index=True)
    all_high_sw_av = pd.concat([all_high_sw_av, det_high_sw_av], ignore_index=True)

    
    
    if not data_sp_avg.empty:
        all_sp_waveform_av = pd.concat([all_sp_waveform_av, data_sp_avg], ignore_index=True)
    # seems like data_sw_avg could corresponds to all_ERP_SW (but it was missing NightID) but I don't find data_sw_avg
    # if not data_sw_avg.empty:
    #     all_sw_waveform_av = pd.concat([all_sw_waveform_av, data_sw_avg], ignore_index=True)
    
    plt.close('all')
            
    # --- Sauvegarde du rapport ---
    report.save(reports_dir + os.sep + "allSubj_microstruct_report.html", overwrite=True, open_browser=False)

# print info of failed_files
print(f"{len(failed_files)} files encountered an error during the feature extraction")
df_failed = pd.DataFrame(failed_files, columns=['file', 'failure step'])
df_failed.to_csv(reports_dir+os.sep +"apomorphee_failed_files.tsv", sep='\t', index=True)

all_sw.to_csv(reports_dir+os.sep +"apomorphee_sw_relaxed_NREM.tsv", sep='\t', index=True)
all_sw_av.to_csv(reports_dir+os.sep +"apomorphee_sw_relaxed_NREM_bystage.tsv", sep='\t', index=True)
all_sp.to_csv(reports_dir+os.sep +"apomorphee_sp_NREM_all.tsv", sep='\t', index=True)
all_sp_av.to_csv(reports_dir+os.sep +"apomorphee_sp_NREM_bystage.tsv", sep='\t', index=True)

all_low_sw.to_csv(reports_dir+os.sep +"apomorphee_low_sw_relaxed_NREM.tsv", sep='\t', index=True)
all_low_sw_av.to_csv(reports_dir+os.sep +"apomorphee_low_sw_relaxed_NREM_bystage.tsv", sep='\t', index=True)
all_high_sw.to_csv(reports_dir+os.sep +"apomorphee_high_sw_relaxed_NREM.tsv", sep='\t', index=True)
all_high_sw_av.to_csv(reports_dir+os.sep +"apomorphee_high_sw_relaxed_NREM_bystage.tsv", sep='\t', index=True)

all_ERP_SW = pd.concat(all_ERP_SW, ignore_index=True)
all_ERP_SW.to_csv(reports_dir+os.sep +"apomorphee_erp_sw_relaxed_NREM_bystage.tsv", sep='\t', index=True)

# Example: aggregate across participants per Drug × Channel × Condition
grand_avg_ERP_SW = all_ERP_SW.groupby(["Drug", "Condition", "Channel", "Time"])["Amplitude"].mean().reset_index()


# In[2]: Plot slow waves waveforms
all_ERP_SW['Channel'] = pd.Categorical(all_ERP_SW['Channel'], categories=['F','C','O'], ordered=True)

SW_grouped = all_ERP_SW.groupby(['Time', 'Channel', 'Drug']).agg(
    Amplitude_mean=('Amplitude', 'mean'),
    Amplitude_sem=('Amplitude', 'sem')
).reset_index()

# Set Seaborn style
sns.set(style='whitegrid')

# Create FacetGrid: one row, 3 columns (Channels), hue=Drug or Drug_Night
g = sns.FacetGrid(
    data=SW_grouped,
    col='Channel',
    hue='Drug',
    sharey=True,
    sharex=True,
    margin_titles=True,
    height=4,
    aspect=1.5
)

# Plot function
def plot_mean_sem(data, color, **kwargs):
    x = data['Time']
    y = data['Amplitude_mean']
    x_sem = data['Amplitude_sem']
    plt.plot(x, y, color=color, **kwargs)
    plt.fill_between(x, y - x_sem, y + x_sem, color=color, alpha=0.3)

# Map plotting function
g.map_dataframe(plot_mean_sem)

# Format labels and titles
g.set_axis_labels("Time (s)", "Amplitude (µV)")
g.add_legend(title='Drug' if 'Drug' in SW_grouped.columns else 'Drug + Night')
g.set_titles(col_template="Channel: {col_name}")

plt.tight_layout()
plt.show()
this_path=plots_path
plt.savefig(this_path + os.sep + "av_SW_byTreament.png", dpi=300, bbox_inches='tight')

SW_grouped = all_ERP_SW.groupby(['Time', 'Channel', 'Drug', 'Night']).agg(
    Amplitude_mean=('Amplitude', 'mean'),
    Amplitude_sem=('Amplitude', 'sem')
).reset_index()

# Create FacetGrid: rows=Night, cols=Channel
g = sns.FacetGrid(
    data=SW_grouped,
    row='Night',
    col='Channel',
    hue='Drug',
    sharey=True,
    sharex=True,
    margin_titles=True,
    height=4,
    aspect=1.5,
)

# Map lineplot with confidence interval manually
def plot_mean_sem(data, color, **kwargs):
    x = data['Time']
    y = data['Amplitude_mean']
    x_sem = data['Amplitude_sem']
    plt.plot(x, y, color=color, **kwargs)
    plt.fill_between(x, y - x_sem, y + x_sem, color=color, alpha=0.3)

g.map_dataframe(plot_mean_sem)

# Add titles and labels
g.set_axis_labels("Time (s)", "Amplitude (µV)")
g.add_legend(title='Drug')
g.set_titles(row_template="Night: {row_name}", col_template="Channel: {col_name}")

plt.tight_layout()
plt.show()
plt.savefig(this_path + os.sep + "av_SW_byNight_byTreament.png", dpi=300, bbox_inches='tight')

# In[3]: Plot spindles waveforms
all_sp_waveform_av['Channel'] = pd.Categorical(all_sp_waveform_av['Channel'], categories=['F','C','O'], ordered=True)

# Step 1: Find combinations where abs(Amplitude) > 100
mask = all_sp_waveform_av['Amplitude'].abs() > 100
combinations_to_remove = all_sp_waveform_av.loc[mask, ['SubID', 'Night', 'Drug']].drop_duplicates()

# Step 2: Create a mask to identify all rows matching these combinations
def is_combination(row):
    return ((combinations_to_remove['SubID'] == row['SubID']) &
            (combinations_to_remove['Night'] == row['Night']) &
            (combinations_to_remove['Drug'] == row['Drug'])).any()

mask_to_nan = all_sp_waveform_av.apply(is_combination, axis=1)

# Step 3: Replace all columns for these rows with NaN
all_sp_waveform_clean_av = all_sp_waveform_av
all_sp_waveform_clean_av.loc[mask_to_nan, :] = np.nan

SP_grouped = all_sp_waveform_clean_av.groupby(['Time', 'Channel', 'Drug']).agg(
    Amplitude_mean=('Amplitude', 'mean'),
    Amplitude_sem=('Amplitude', 'sem')
).reset_index()

# Set Seaborn style
sns.set(style='whitegrid')

# Create FacetGrid: one row, 3 columns (Channels), hue=Drug or Drug_Night
g = sns.FacetGrid(
    data=SP_grouped,
    col='Channel',
    hue='Drug' if 'Drug' in SP_grouped.columns else 'Drug_Night',
    sharey=True,
    sharex=True,
    margin_titles=True,
    height=4,
    aspect=1.5
)

# Plot function
def plot_mean_sem(data, color, **kwargs):
    x = data['Time']
    y = data['Amplitude_mean']
    x_sem = data['Amplitude_sem']
    plt.plot(x, y, color=color, **kwargs)
    plt.fill_between(x, y - x_sem, y + x_sem, color=color, alpha=0.3)

# Map plotting function
g.map_dataframe(plot_mean_sem)

# Format labels and titles
g.set_axis_labels("Time (s)", "Amplitude (µV)")
g.add_legend(title='Drug' if 'Drug' in SP_grouped.columns else 'Drug + Night')
g.set_titles(col_template="Channel: {col_name}")

plt.tight_layout()
plt.show()
plt.savefig(this_path + os.sep + "av_SP_byTreament.png", dpi=300, bbox_inches='tight')

SP_grouped = all_sp_waveform_clean_av.groupby(['Time', 'Channel', 'Drug', 'Night']).agg(
    Amplitude_mean=('Amplitude', 'mean'),
    Amplitude_sem=('Amplitude', 'sem')
).reset_index()

# Create FacetGrid: rows=Night, cols=Channel
g = sns.FacetGrid(
    data=SP_grouped,
    row='Night',
    col='Channel',
    hue='Drug',
    sharey=True,
    sharex=True,
    margin_titles=True,
    height=4,
    aspect=1.5,
)

# Map lineplot with confidence interval manually
def plot_mean_sem(data, color, **kwargs):
    x = data['Time']
    y = data['Amplitude_mean']
    x_sem = data['Amplitude_sem']
    plt.plot(x, y, color=color, **kwargs)
    plt.fill_between(x, y - x_sem, y + x_sem, color=color, alpha=0.3)

g.map_dataframe(plot_mean_sem)

# Add titles and labels
g.set_axis_labels("Time (s)", "Amplitude (µV)")
g.add_legend(title='Drug')
g.set_titles(row_template="Night: {row_name}", col_template="Channel: {col_name}")

plt.tight_layout()
plt.show()
plt.savefig(this_path + os.sep + "av_SP_byNight_byTreament.png", dpi=300, bbox_inches='tight')

#%% same plot as the plot from EDF_extract_features_amp75_251006.py
# load data if you don't want to re-run the code
all_ERP_SW = pd.read_csv(reports_dir+os.sep +"apomorphee_erp_sw_relaxed_NREM_bystage.tsv", sep='\t')

all_ERP_SW['Channel'] = pd.Categorical(all_ERP_SW['Channel'], categories=['F','C','O'], ordered=True)

SW_grouped = all_ERP_SW.groupby(['Time', 'Channel', 'Drug', 'Condition']).agg(
    Amplitude_mean=('Amplitude', 'mean'),
    Amplitude_sem=('Amplitude', 'sem')
).reset_index()


# Create FacetGrid: rows=Night, cols=Channel
g = sns.FacetGrid(
    data=SW_grouped,
    row='Condition',
    col='Channel',
    hue='Drug',
    sharey=True,
    sharex=True,
    margin_titles=True,
    height=4,
    aspect=1.5,
)

# Map lineplot with confidence interval manually
def plot_mean_sem(data, color, **kwargs):
    x = data['Time']
    y = data['Amplitude_mean']
    x_sem = data['Amplitude_sem']
    plt.plot(x, y, color=color, **kwargs)
    plt.fill_between(x, y - x_sem, y + x_sem, color=color, alpha=0.3)

g.map_dataframe(plot_mean_sem)

# Add titles and labels
g.set_axis_labels("Time (s)", "Amplitude (µV)")
g.add_legend(title='Drug')
g.set_titles(row_template="Condition: {row_name}", col_template="Channel: {col_name}")

plt.tight_layout()
plt.show()
this_path = plots_path
plt.savefig(this_path + os.sep + "av_SW_byCond_byTreament.png", dpi=300, bbox_inches='tight')




