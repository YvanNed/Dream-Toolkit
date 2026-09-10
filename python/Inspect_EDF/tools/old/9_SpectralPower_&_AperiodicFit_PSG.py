# -*- coding: utf-8 -*-
"""
This script is an example to compute spectral power, aperiodic components, and extract bandpower for classic PSG recordings (~3 EEG channels).
It was developped on the apomorphee database.
It requires raw .edf EEG files, hypnograms, and output of "2_select&remap_channels_edf_voila.ipynb" to harmonize channel names if needed.
It works mainly at the epochs level. 

Steps:
    - reject epochs based on peak-to-peak amplitude
    - compute Power Spectrum Density (welch method, 4s sliding window, fmin 0.5 Hz, fmax 30 Hz)
    - compute aperiodic components (fit from 2-30Hz, then apply it to full spectrum to remove it)
    - extract frequency band power for linear (µV²/Hz), log (dB/Hz), and aperiodic-removed PSDs
    
Outputs (in a newly created folder reports_spectralpower):
    - individual reports with hypno/spectrogram, PSDs per sleep stage (dB/Hz), and aperiodic-removed PSDs per sleep stage (dB/Hz) plots
    - .tsv tables (can be open with excel) for frequency band power, aperiodic components, and periodic components (from SpecParam) at the epoch level
    - group-level plots of classic/aperiodic-removed PSDs per sleep stage
    - additional .tsv files to check faiuled files to load, % of rejected epochs, failed 1/f fit to PSD
 
@author: yvan.nedelec
"""

#%% Import package and define custom function
# packages
import os
import re
import mne
import yasa
import glob
import json
import sklearn
import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path 
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
from specparam import SpectralModel
from statsmodels.stats.multitest import fdrcorrection
from sklearn.metrics import accuracy_score, f1_score, cohen_kappa_score

# %matplotlib qt

# custom functions (required to apply the remapping from select&remap notebook)
def drop_suffix_duplicates(raw):
    """
    Repère les channels de type 'XXX-0', 'XXX-1', ...
    et ne garde que les '-0' lorsqu'ils coexistent.
    
    Exemple :
    - ['C4-0', 'C4-1', 'F4', 'O2'] -> on drop 'C4-1'
    - ['C4-1'] tout seul -> on ne drop rien.
    """
    ch_names = raw.ch_names
    # regroupement par 'base', ex: C4-0 -> base = C4
    groups = {}
    for ch in ch_names:
        if ch.endswith('-0') or ch.endswith('-1'):
            base = ch.rsplit('-', 1)[0]
            groups.setdefault(base, []).append(ch)

    to_drop = []
    for base, ch_list in groups.items():
        zeros = [c for c in ch_list if c.endswith('-0')]
        if zeros:
            # si un '-0' existe, on drop tous les autres suffixés (typiquement les '-1')
            for c in ch_list:
                if not c.endswith('-0'):
                    to_drop.append(c)

    if to_drop:
        print("🔎 Channels suffixés trouvés, on va drop : ", to_drop)
        raw.drop_channels(to_drop)
    else:
        print("✅ Aucun doublon '-0'/'-1' à nettoyer")

    return raw, to_drop

def adapt_remap_dict_to_suffixes(raw, remap_dict):
    """
    Adapte un remap_dict de la forme {"C4": "C", "F4": "F", ...}
    aux noms de channels réellement présents dans raw.ch_names,
    en remplaçant par exemple "C4" -> "C4-0" si "C4" n'existe pas
    mais "C4-0" oui.

    Retourne un nouveau dict prêt pour raw.rename_channels().
    """
    ch_set = set(raw.ch_names)
    new_remap = {}

    for base_label, target in remap_dict.items():
        if base_label in ch_set:
            # Cas simple : le label existe tel quel dans raw.ch_names
            new_remap[base_label] = target
        else:
            # Cas suffixé : essayer base_label-0
            candidate = f"{base_label}-0"
            if candidate in ch_set:
                new_remap[candidate] = target
            # sinon : on ignore ce label (aucun ch correspondant)

    return new_remap

# value used to clamp psd if = 0 or close to zero
eps = np.finfo(float).tiny

#%% define data path
if 'thandrillon' in os.getcwd():
    path_data = '/Users/thandrillon/Data/Apomorphee'
elif 'yvan.nedelec' in os.getcwd():
    path_data = '/Users/yvan.nedelec/OneDrive - ICM/Documents/RE_yvan/projects/noemie_rescue'
else:
    path_data = '/Users/noemiewane/Library/Mobile Documents/com~apple~CloudDocs/Internat/Master 2/Stage/Données/Version 4/'

# Création du dossier "reports"
reports_dir = os.path.join(path_data, "reports_spectralpower")
os.makedirs(reports_dir, exist_ok=True)

# Création du dossier "plots"
plots_path = os.path.join(path_data, "plots_spectralpower")
os.makedirs(plots_path, exist_ok=True)

# création du dossier derivatives (this will depend on the data structure i.e. BIDS vs not BIDS)
deriv_path = os.path.join(path_data, "derivatives")
os.makedirs(deriv_path, exist_ok=True)

# Recherche des fichiers EDF
files = glob.glob(os.path.join(path_data+os.sep+'data', '*.edf'))
redo = 1 # 0 = skip participants with report already computed (do not use yet because saved data tables will not contain those participants); 1 = redo participant even if they have a report (mandatory for now)

# load json file for remap channels and re-ref
json_path = "C:/Users/yvan.nedelec/OneDrive - ICM/Documents/RE_yvan/projects/noemie_rescue/data/config_param/remap_reref_persubject.json"
with open(json_path, 'r', encoding='utf-8') as f: config_dict = json.load(f)

# participant additional info
subj_table = pd.read_csv(os.path.join(path_data, 'CondTable.csv'), sep=';')

#______________________________________________________________________________
# to select a specific file based on subject label (debug or test feature)_____
sub_file = "10_N1"
# create a pattern as a regex (regular expression) object that allows operation such as search, match etc. 
# we add (?<!\d) before the subject name to avoid subject names with digits before the subject name of interest
pattern = re.compile(rf"(?<!\d){re.escape(sub_file)}") # the re.escape allows to avoid errors if there are special characters (e.g. .,+,-,(,) etc) in the string of interest (you have to add a \ to avoid the function of the special character). if you write "(pdf)", re.escape("(pdf)") will return "\(pdf\)"
# search for the pattern 
mask_file = np.array([bool(pattern.search(f)) for f in files])
file = files[np.where(mask_file)[0][0]]

# ==> then run the codes within the loop of participants 
#______________________________________________________________________________

#______________________________________________________________________________
# to compute the code on missing reports_______________________________________
# check files with report already computed
report_files = glob.glob(os.path.join(path_data+os.sep+'reports', '*.html'))
# get participants name from report already computed
report_names = {
    os.path.normcase(Path(rf).stem.replace("_report_report", ""))
    for rf in report_files
}

missing_files = [
    f
    for f in files
    if os.path.normcase(Path(f).stem) not in report_names
]
#______________________________________________________________________________

#%% loop over participants
failed_files = [] # list to catch file encountering an error
rejected_epo = pd.DataFrame() # df to catch the rejected epochs
all_bandpower_df = pd.DataFrame() # df to save bandpower, here we could use a method to load an existing file, but then we would have to be careful when we re-run the script over and over)
psds_list = [] # list to preserve individual psds data
ap_removed_psds_list = [] # list to preserve individual psds after removing the aperiodic comp. 
ap_offset_list = [] # list to preserve individual aperiodic offset
ap_exponent_list = [] # list to preserve individual aperiodic exponent
all_ap_df = pd.DataFrame() # df to save aperiodic components 
failed_fit = [] # list to catch when specparam failed to fit 1/f
all_periodic_df = pd.DataFrame() # df to save periodic components from SpecParam

for file in files: # change files to missing_files to compute only missing files (DON'T USE IT FOR NOW: .xlsx files will be overwritten and the group average plots will not work)
    
    #__________________________________________________________________________
    # ------------------------------ file info --------------------------------
    file_name = os.path.basename(file)  # Ex: "15_N1.edf"
    file_ID = file_name.split('.')[0]  # Extrait "15_N1"
    sub_ID = file_ID.split("_")[0]
    # night_ID = file_ID.split("_")[-1]
    # drug_ID = subj_table.loc[subj_table['SubID'] == int(sub_ID), night_ID].values[0]
    # group_ID = 'G_N1_'+subj_table.loc[subj_table['SubID'] == int(sub_ID), 'N1'].values[0]
    
    # check if report and create it 
    report_prefix = os.path.join(reports_dir, f"{file_ID}_report")
    if os.path.exists(report_prefix + ".html") and redo == 0:
        print(f"🚀 Rapport déjà généré pour {file_name}, passage au suivant...")
        continue
    report = mne.Report(title=file_name)
    # -------------------------------------------------------------------------
    #__________________________________________________________________________
    
    #__________________________________________________________________________
    # ------------------------------ load EEG ---------------------------------
    # get participant info from json
    sub_config = config_dict[file_ID]
    selected_channels = list(sub_config["remap"].keys())
    
    try:
        raw = mne.io.read_raw_edf(file, preload=True, encoding="latin-1", include=selected_channels) # load only selected channels from the notebook
    except Exception as e:
        print(f"⚠️ Erreur lecture fichier {file_name} : {e}")
        failed_files.append((file, 'eeg loading'))
        continue
    
    # remove duplicates channels (XX-0, XX-1 => will keep only XX-0)
    raw, dropped = drop_suffix_duplicates(raw)
    
    # remap channels (based on the generated dict from jupy notebook):
    # remap dict to take into account XX-0 created by mne for duplicates
    remap_dict_adapted = adapt_remap_dict_to_suffixes(raw, sub_config["remap"])
    raw.rename_channels(remap_dict_adapted)
    
    # set eeg reference
    if  sub_config["ref_channels"] != 'average':
        raw.set_eeg_reference(ref_channels = sub_config["ref_channels"])
        # get rid of the ref
        raw.drop_channels(sub_config["ref_channels"])
    else:
        raw.set_eeg_reference(ref_channels = sub_config["ref_channels"])
    
    # get sf
    sf = raw.info["sfreq"]
    
    # bandpass filtering
    # from 0.1 Hz to 40 Hz (two-ways FIR filter, Gramfort, 2013)  
    raw.filter(l_freq = 0.1, h_freq = 40, # 0.1 was chosen because we will compute delta band power from 0.5 to 4 Hz
               l_trans_bandwidth= 'auto', h_trans_bandwidth= 'auto',
               filter_length='auto', method = 'fir', phase = 'zero-double', 
               fir_window='hamming', fir_design='firwin')
    # -------------------------------------------------------------------------
    
    # --------------------------- load hypnogram ------------------------------
    hypno_name = os.path.join(path_data+os.sep+'data', file_ID+"_Hypnogram_remapped.txt") # those hypnogram have been harmonized with the check_hypno.py script 
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
    
    # Vérification de la compatibilité (longueur EEG vs longueur hypno)
    expected_epochs = int(np.floor(raw.n_times / sf / 30))
    if expected_epochs != len(expert_hypno):
        print(f"⚠️ Incompatibilité hypnogramme ({len(expert_hypno)}) vs EEG ({expected_epochs} epochs)")
        failed_files.append((file, 'mismatch hypno - eeg'))
        continue
    
    # Conversion in numerical vector 
    hypno_vec = np.full(len(expert_hypno), np.nan)
    mapping = {'W': 0, 'N1': 1, 'N2': 2, 'N3': 3, 'R': 4}
    for stage, value in mapping.items():
        hypno_vec[np.where(expert_hypno == stage)] = value
    # -------------------------------------------------------------------------
    
    # ---------------- Quick plot of the hypno/spectrogram --------------------
    for c, ch in enumerate(raw.ch_names):
        data, _ = raw[c, :]
        hypno_fullvec = yasa.hypno_upsample_to_data(
            hypno=hypno_vec, 
            sf_hypno=1/30, 
            data=data[0], 
            sf_data=sf
        )
        
        fig = yasa.plot_spectrogram(data[0], sf, hypno_fullvec, fmin=0.5, fmax=40)
        report.add_figure(fig=fig, title=f"Hypnogram and Spectrogram {ch}", image_format="PNG")
        plt.close(fig)
    # -------------------------------------------------------------------------
    #__________________________________________________________________________
    
    #__________________________________________________________________________
    # -------------- compute epochs and assign the sleep stage ----------------
    # quick look at epochs distribution
    # hyp_values, hyp_counts = np.unique(expert_hypno, return_counts=True)
    # for h, c in zip(hyp_values, hyp_counts):
    #     print(f"{h}: {c}")
    
    # define parameters (would be interesting to test 5s epochs like Pierre Champetier)
    epoch_dur = 30 # seconds
    ptp_criteria = dict(eeg=250e-6)  # criteria peak-to-peak amplitude of 200 µV to reject epochs. Could be specific for each sleep stage (e.g. go up to 300 µV for N3)
    flat_criteria = dict(eeg=1e-6) # criteria to reject epoch with peak-to-peak amp. < criteria (1 µV here)
    # compute epochs
    epochs = mne.make_fixed_length_epochs(raw, duration=epoch_dur, preload=True)
    
    # add sleep stage as event
    epochs.events[:, 2] = hypno_vec     # add sleep stage as events
    epochs.event_id = mapping
    
    # drop bad epochs and log the rejected epochs. We could also add a reject based on sample-to-sample amplitude (to identify pops and jumps)
    tot_epo = np.array(list({event: len(epochs[event]) for event in epochs.event_id.keys()}.values())) # extract the total number of epochs per sleep stage
    
    epochs.drop_bad(reject=ptp_criteria, flat=flat_criteria) # drop epochs with peak-to-peak amp > 200 µV (absolute difference between min and max of the epoch) and flat signal (if peak-to-peak amp. < 1 µV)
    
    good_epo = np.array(list({event: len(epochs[event]) for event in epochs.event_id.keys()}.values())) # extract the number of epochs after rejection per sleep stage
    bad_epo = tot_epo - good_epo 
    bad_prop = bad_epo / tot_epo
    stage_list = list(epochs.event_id.keys())
    rej_epo = pd.DataFrame({"FileID": file_ID, "Stage": stage_list, "bad epochs": bad_epo, "bad proportion": bad_prop, "total epochs": tot_epo}) # store bad epochs info as df
    rejected_epo = pd.concat([rejected_epo, rej_epo], ignore_index=True) # concatenate individual data to group data
    # -------------------------------------------------------------------------
    #__________________________________________________________________________
    
    #__________________________________________________________________________
    # ----------------------- compute spectral power --------------------------
    # ------------------------ (periodic activity) ----------------------------
    # define parameters
    method = 'welch'
    fmin = 0.5 # Hz
    fmax = 30 # Hz
    window_dur = 4 # 4 seconds sliding window lead to a frequency resolution of 0.25 Hz. Our lowest freq. is 0.5 and the rule of thumb is to have a resoluton that encompasses at least 2 cycle of the lowest freq. (so freq. res. of 0.25 Hz)
    n_per_seg = int(window_dur*sf) # taille de la FFT window en sample
    n_overlap = int(n_per_seg/2) # 50% overlap
    window = "hann" # hanning window, most commonly used in sleep research
    
    # compute psd (0.5-30 Hz with a freq res of 0.25 Hz should results in 119 frequency bins)
    psds = epochs.compute_psd(method = method, fmin=fmin, fmax=fmax, n_fft = n_per_seg, 
                              n_overlap = n_overlap, n_per_seg = n_per_seg, window = window)
    # preserve individual data for group plot
    psds_list.append(psds)
    
    # extract psd data and convert to µV. 
    psds_data = psds.get_data() * 1e12
    psdsdB_data = 10*np.log10(psds_data.copy()) # multiplying by 10 is to turn into dB
    freqs = psds.freqs
    
    # save psds (use mne.time_frequency.read_spectrum() to load back) (could be used depending on the next part of the pipeline)
    # psds.save(os.path.join(deriv_path,file_ID + '-psd.hdf5'))
    
    # get the code of sleep stages still in the data (after epoch rejection) and extract the array of the associated epochs (required for looping across epochs for 1/f fitting, and plotting)
    inv_event_id = {v: k for k, v in psds.event_id.items()} # get the code of event_id (corresponding to sleep stage)
    epoch_stages = np.array([inv_event_id[e] for e in psds.events[:, 2]]) # create array with the sleep stage of the remaining epochs (to loop over it)
    # -------------------------------------------------------------------------
    
    # ------------------ compute 1/f offset and exponent ----------------------
    # ------------------------ (aperiodic activity) ---------------------------
    
    # Here, we decided to fit the 1/f only on the 2-30 Hz range to minimize the influence of low-frequency artifacts and physiological non-stationarities (increase of slow waves in deep sleep)
    # Then, we extrapolate the 1/f to the full spectrum (0.5-30 Hz) to remove the aperiodic component to the full spectrum.
    # Therefore, interpretation of delta band should be careful, and sleep slow waves can be analyzed as events (in the temporal dimension)
    # if the spectrum goes to higher frequencies (>40 Hz), you can consider using aperiodic_mode = knee to take into account a 'bend' in the aperiodic component (but the method to remove the aperiodic comp. to spectrum needs to be changed) 
    
    # specparam algo parameters
    model_param = dict(
        peak_width_limits = [0.5, 20], # min = >2*freq. resolution (0.25 here); max = 20 Hz based on Sophia recommendation (default is 12 Hz)
        aperiodic_mode = 'fixed', # can use 'knee' if the loglog plot is not a straight line (the line have an angle at some point) or 'doublexp' if loglog is curvated see: https://specparam-tools.github.io/auto_tutorials/plot_05-AperiodicFitting.html#sphx-glr-auto-tutorials-plot-05-aperiodicfitting-py
        min_peak_height = 0.3 # to avoid capturing very small peaks and overfit (can be changed if too many/too few peaks are detected)
    )
    
    specparam_model = SpectralModel(**model_param) # initialyze model
    
    # prepare numpy arrays to store individual data
    ap_removed_psds_data = np.full(psds_data.shape, np.nan) # epochs x channels x freq. bin
    ap_offset_data = np.full((len(epoch_stages), len(psds.ch_names)), np.nan) # epochs x channels
    ap_exponent_data = np.full((len(epoch_stages), len(psds.ch_names)), np.nan) # epochs x channels
    rows_ap = []
    rows_periodic = []
    
    # loop across epochs and channels to fit the model
    for ei, stage in enumerate(epoch_stages):
        for ci, ch in enumerate(psds.ch_names):
            
            # select current psd (for the epoch and channel)
            cur_psd = psds_data[ei,ci,:]
            # reduced the fit only on freqs > 2 Hz
            mask = freqs >= 2
            freqs_reduced = freqs[mask]
            psd_reduced = cur_psd[mask]
            
            # fit
            try:
                # psd_reduced = np.maximum(psd_reduced, eps) # clamp psd to a minimal value if it contains zeros or negative values (not sure it is recommended and should not be necessary if you have clean data)
                specparam_model.fit(freqs_reduced, psd_reduced)
            except Exception as e:
                print(f"⚠️ Erreur Specparam fit pour {file_ID}, epoch {ei} : {e}")
                failed_fit.append((file, f'epoch {ei}, channel {ci}'))
                continue
            
            # get aperiodic comp.
            offset, exponent = specparam_model.get_params("aperiodic")
            
            # preserve individual data
            ap_offset_data[ei, ci] = offset
            ap_exponent_data[ei, ci] = exponent
            # extrapolate the fit to the range 0.5-30 Hz
            extrapolated_fit = offset - exponent*np.log10(freqs) # works only for aperiodic_mode = 'fixed', need to come up with another solution if 'knee'
            
            # get fitting metrics
            error = specparam_model.get_metrics('error', 'mae') # Mean absolute error of the model fit to the data.
            r_squared = specparam_model.get_metrics('gof', 'squared') # R-squared between the model fit and the data.
            
            # compute df
            rows_ap.append({
                "file_ID": file_ID,
                "epoch": int(ei),
                "stage": stage,
                "channel": ch,
                "offset": offset,
                "exponent": exponent,
                "r_squared": r_squared,
                "error": error,
                })
            
            # Remove the aperiodic comp. to psds_log 
            extrapolated_fit_dB = 10*extrapolated_fit # convert to dB
            ap_removed_psd = psdsdB_data[ei,ci,:] - extrapolated_fit_dB
            ap_removed_psds_data[ei, ci, :] = ap_removed_psd
            
            # get periodic info and store into another table
            periodic_peaks = specparam_model.results.get_params('periodic', version='converted')
            for pi, peak_info in enumerate(periodic_peaks):
                # print(pi, peak_info)
                rows_periodic.append({
                    "file_ID": file_ID,
                    "epoch": int(ei),
                    "stage": stage,
                    "channel": ch,
                    "peak_#": pi+1,
                    "total_peak_#": len(periodic_peaks),
                    "central_frequency": peak_info[0],
                    "peak_power": peak_info[1],
                    "peak_bandwidth": peak_info[2],
                    })
            
    
    # save aperiodic in dataframe
    df_ap = pd.DataFrame(rows_ap)
    all_ap_df = pd.concat([all_ap_df, df_ap], ignore_index=True)
    
    # save perdiodic info in dataframe
    df_periodic = pd.DataFrame(rows_periodic)
    all_periodic_df = pd.concat([all_periodic_df, df_periodic], ignore_index=True) # Sophia does not recommend to use the peak extraction from Specparam. Moreover, here it is at the epoch level, might be better to work at the sleep stage level (smoother spectrum)

    # store individual psds after removing aperiodic comp.
    ap_removed_psds = psds.copy()
    ap_removed_psds._data = ap_removed_psds_data # replace data in the mne.time_frequency.spectrum object 
    ap_removed_psds_list.append(ap_removed_psds)
    ap_offset_list.append(ap_offset_data) # as np array for now, could be considered to put in an mne object to use mne method
    ap_exponent_list.append(ap_exponent_data) # as np array for now, could be considered to put in an mne object to use mne method
    
    # save psds with removed aperiodic comp. (already in dB) (use mne.time_frequency.read_spectrum() to load back) (could be used depending on the next part of the pipeline)
    # ap_removed_psds.save(os.path.join(deriv_path,file_ID + '-psd.hdf5'))
    # -------------------------------------------------------------------------
    
    # ------------- compute band power (absolute and relative) ----------------
    # -------- from classical PSD and PSD with aperiodic comp. removed --------
    bands = {
        'delta': (0.5, 4),
        'theta': (4, 8),
        'alpha': (8, 12),
        'sigma': (12, 16),
        'beta': (16, 30)
    }
    
    # compute total power (to compute relative power per bands)
    total_power = psds_data.sum(axis=2) # sum of power across the whole PSD per epochs and channels, not focusing only in the defined band power, eventhough in our case the psd was computed only on our band power of interest (0.5 to 30 Hz)
    
    # compute band power in linear and log space at the epoch level
    bandpower_lin = {}
    bandpower_db = {}
    bandpower_ap_removed = {}
    for band, (fmin, fmax) in bands.items():
        idx = np.logical_and(freqs >= fmin, freqs <= fmax)
        bandpower_lin[band] = np.nanmean(psds_data[:, :, idx], axis=2)                      # µV²/Hz
        bandpower_db[band] = np.nanmean(psdsdB_data[:, :, idx], axis=2)                     # dB/Hz
        bandpower_ap_removed[band] = np.nanmean(ap_removed_psds_data[:, :, idx], axis=2)    # dB/Hz (with aperiodic comp. removed)
    
    # create empty list, that will be filled with dict containing relevant info and then turn to dataframe
    rows_bp = [] 
    
    for band in bands:
        pow_lin = bandpower_lin[band]                   # (n_epochs, n_ch)
        pow_db = bandpower_db[band]                     # (n_epochs, n_ch)
        pow_rel = pow_lin / total_power                 # relative per epoch
        pow_ap_removed = bandpower_ap_removed[band]     # (n_epochs, n_ch)
        
        for ei, stage in enumerate(epoch_stages):
            for ci, ch in enumerate(psds.ch_names):      
                rows_bp.append({
                    "file_ID": file_ID,
                    "epoch": int(ei),
                    "stage": stage,
                    "band": band,
                    "channel": ch,
                    "power_uV2": pow_lin[ei, ci],
                    "power_dB": pow_db[ei, ci],
                    "power_relative": pow_rel[ei, ci],
                    "power_aperiodic_removed_dB": pow_ap_removed[ei, ci],
                    })
                
    df_bandpower = pd.DataFrame(rows_bp)
    all_bandpower_df = pd.concat([all_bandpower_df, df_bandpower], ignore_index=True) 
    # -------------------------------------------------------------------------  
    #__________________________________________________________________________
    
    #__________________________________________________________________________
    # -------------------------- figures for report --------------------------- 
    # Create ONE figure with n_channels (should be 3) subplots horizontally aligned
    # for now I am keeping the figure averaging after passing in the log space since it reflects the data used for stats (epochs, already in log space)
    
    # ------- plot spectral power per channel for each sleep stage (ss) -------
    ss_colors = {'W':'#969696', 'N1':'#9e9ac8', 'N2':'#807dba', 'N3':'#6a51a3', 'R':'#c994c7'} # in order of W, N1, N2, N3, R
    ss_order = ['W', 'N1', 'N2', 'N3', 'R']# should correspond to event_id.keys(), allows to select which state you want to analyze
    ss_order = [s for s in ss_order if s in psds.event_id]  # only keep stages in the data
    
    n_channels = len(psds.ch_names)
    
    fig1, axes1 = plt.subplots(1, n_channels, figsize=(5 * n_channels, 5), sharey=True)

    # If only 1 channel, axes is not a list → enforce list
    if n_channels == 1:
        axes1 = [axes1]
    
    # loop through channel
    for ch_idx, ch_name in enumerate(psds.ch_names):
    
        ax1 = axes1[ch_idx]
        
        # loop through sleep stage (via events and event_id)
        for ss_idx, ss_name in enumerate(ss_order):
            ss_code = psds.event_id[ss_name] # get the sleep stage code corresponding to sleep stage name
            ss_color = ss_colors[ss_name]
    
            # mask to select psd epochs corresponding to this sleep stage
            mask = (psds.events[:,2] == ss_code)
            
            if not np.any(mask):
                continue  # aucune epoch pour ce stade
    
            # to plot from linear space     
            # # select PSD for this channel and this sleep stage: (n_epochs_stage, n_freqs)
            # data_stage = psds_data[mask, ch_idx, :]
    
            # mean_lin = data_stage.mean(axis=0)
            # sem_lin = data_stage.std(axis=0, ddof=1) / np.sqrt(data_stage.shape[0])
            # lower_lin = lower_lin = np.maximum(mean_lin - sem_lin, np.finfo(float).eps) # avoid values <= 0 (log does not like them)
            # upper_lin = mean_lin + sem_lin
            
            # compare when average across log (because stats will likely be done on )
            data_stage_log = psdsdB_data[mask, ch_idx, :]
            
            mean_log = data_stage_log.mean(axis=0)
            sem_log = data_stage_log.std(axis=0, ddof=1) / np.sqrt(data_stage_log.shape[0])
            lower_log = mean_log - sem_log
            upper_log = mean_log + sem_log
            
            ax1.plot(freqs, mean_log, label=f'LOG {ss_name} (n={data_stage_log.shape[0]})', color = ss_color) # change for mean_lin for linear space
    
            # Shading pour la SEM
            ax1.fill_between(freqs, lower_log, upper_log,
                             alpha=0.3, color = ss_color)       # change for lower_lin, upper_lin for linear space
        
        ax1.set_title(f'{ch_name}')
        ax1.set_xlabel('Frequency (Hz)', fontsize=12)
        ax1.set_ylabel('PSD (dB/Hz)', fontsize=12)   # change to µV²/Hz if linear space
        ax1.set_xlim(freqs[0], freqs[-1])
        
        # to plot in log-log (put x axis in log as well)
        # ax1.set_xscale('log')
        # ax1.set_xlim(freqs[0], freqs[-1])
        # ax1.set_xticks([0.5, 1, 2, 4, 8, 16, 30]) # change x ticks for better visualization in log
        # ax1.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        
        ax1.legend()
        
    fig1.tight_layout()
    
    report.add_figure(fig=fig1, title="PSD per sleep stage", image_format="PNG")
    plt.close(fig1)
    # -------------------------------------------------------------------------
    
    # ---------------------- plot psd after 1/f removal -----------------------
    # Create ONE figure with n_channels (should be 3) subplots horizontally aligne
    fig2, axes2 = plt.subplots(1, n_channels, figsize=(5 * n_channels, 5), sharey=True)

    # If only 1 channel, axes is not a list → enforce list
    if n_channels == 1:
        axes2 = [axes2]
    
    # loop through channel
    for ch_idx, ch_name in enumerate(ap_removed_psds.ch_names):
        ax = axes2[ch_idx]

        # loop through sleep stage (via events and event_id)
        for ss_idx, ss_name in enumerate(ss_order):
            ss_code = ap_removed_psds.event_id[ss_name] # get the sleep stage code corresponding to sleep stage name
            ss_color = ss_colors[ss_name]
    
            # mask to select psd epochs corresponding to this sleep stage
            mask = (ap_removed_psds.events[:,2] == ss_code)
            
            if not np.any(mask):
                continue  # aucune epoch pour ce stade
    
            # select PSD for this channel and this sleep stage: (n_epochs_stage, n_freqs), already in log
            data_stage = ap_removed_psds_data[mask, ch_idx, :]
            
            mean_log = data_stage.mean(axis=0)
            sem_log = data_stage.std(axis=0, ddof=1) / np.sqrt(data_stage.shape[0])
            lower_log = mean_log - sem_log
            upper_log = mean_log + sem_log
            
            ax.plot(freqs, mean_log, label=f'LOG {ss_name} (n={data_stage.shape[0]})', color = ss_color)
    
            # Shading pour la SEM
            ax.fill_between(freqs, lower_log, upper_log,
                             alpha=0.3, color = ss_color)

    
        ax.set_title(f'{ch_name}')
        ax.set_xlabel('Frequency (Hz)', fontsize=12)
        ax.set_ylabel('PSD (dB/Hz)', fontsize=12)   # change to µV²/Hz if linear space
        ax.set_xlim(freqs[0], freqs[-1])
        
        # to plot in log-log (put x axis in log as well)
        # ax.set_xscale('log')
        # ax.set_xlim(freqs[0], freqs[-1])
        # ax.set_xticks([0.5, 1, 2, 4, 8, 16, 30]) # change x ticks for better visualization in log
        # ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        
        ax.legend()
        
    fig2.tight_layout()
    
    report.add_figure(fig=fig2, title="PSD after removing aperiodic component per sleep stage", image_format="PNG")
    plt.close(fig2)
    # -------------------------------------------------------------------------
    #__________________________________________________________________________
    
    # add plot of the frequency band topo ?
    
    # ------------------------ Sauvegarde du rapport --------------------------
    report.save(report_prefix + ".html", overwrite=True, open_browser=False)
    # -------------------------------------------------------------------------

rejected_epo.to_csv(reports_dir+os.sep +"SpectralAnalysis_drop-log_epoch.tsv", sep='\t', index=True)
all_bandpower_df.to_csv(reports_dir+os.sep +"SpectralAnalysis_bandpower_per_epoch.tsv", sep='\t', index=True)
all_ap_df.to_csv(reports_dir+os.sep +"SpectralAnalysis_aperiodic_per_epoch.tsv", sep='\t', index=True)
df_failed_files = pd.DataFrame(failed_files)
df_failed_files.to_csv(reports_dir+os.sep +"SpectralAnalysis_failed_files_log.tsv", sep='\t', index=True)
df_failed_fit = pd.DataFrame(failed_fit)
df_failed_fit.to_csv(reports_dir+os.sep +"SpectralAnalysis_failed_fit_log.tsv", sep='\t', index=True)
all_periodic_df.to_csv(reports_dir+os.sep +"SpectralAnalysis_periodic_per_epoch.tsv", sep='\t', index=True)

#______________________________________________________________________________
# ------------------------- plot at the group level ---------------------------
if len(psds_list) == 0:
    print("⚠️ No PSD in psds_list, can't plot at the group level.")
else:
    # All PSDs should have the same frequencies and channels
    freqs = psds_list[0].freqs
    ch_names = psds_list[0].ch_names
    n_channels = len(ch_names)

    # optional check : same freqs / same ch_names for all
    for psds in psds_list[1:]:
        if not np.allclose(freqs, psds.freqs):
            raise ValueError("Frequencies between participants are not consistent, please harmonize.")
        if psds.ch_names != ch_names:
            raise ValueError("Channel names between participants are not consistent, please harmonize.")

    # define sleep stage name & colors
    ss_colors = {'W':'#969696', 'N1':'#9e9ac8', 'N2':'#807dba', 'N3':'#6a51a3', 'R':'#c994c7'} # in order of W, N1, N2, N3, R
    ss_order = ['W', 'N1', 'N2', 'N3', 'R']

    # Keep only sleep stage within at least one psds.event_id
    all_keys = set().union(*[set(p.event_id.keys()) for p in psds_list])
    ss_order = [s for s in ss_order if s in all_keys]

    # re-organize data as group_data[stage][ch_idx] = participant mean of PSD for each stage and channel
    group_data_log = {
        stage: {ch_idx: [] for ch_idx in range(n_channels)}
        for stage in ss_order
    }
    
    group_data_lin = {
        stage: {ch_idx: [] for ch_idx in range(n_channels)}
        for stage in ss_order
    }

    for pi, psds in enumerate(psds_list):
        # PSD en linéaire (µV²/Hz)
        psds_lin = psds.get_data() * 1e12  # shape: (n_epochs, n_channels, n_freqs)
        psds_lin = np.where(psds_lin > 0, psds_lin, np.nan) # put nan where psds <= 0
        psds_log = 10*np.log10(np.maximum(psds_lin.copy(), eps)) # clamp psd_lin with a minimal value when psd <= 0
        # psds_log = 10*np.log10(psds_lin)
        events = psds.events[:, 2]         # codes de stades (entiers)
        event_id = psds.event_id           # mapping nom -> code
        
        for stage in ss_order:
            if stage not in event_id:
                continue  # ce sujet n'a pas ce stade du tout dans ses event_id

            ss_code = event_id[stage]
            mask = (events == ss_code)

            if not np.any(mask):
                # ce sujet n'a pas d'epoch dans ce stade
                continue

            # moyenne par stade ET par canal (on garde la dimension freq)
            # résultat: (n_channels, n_freqs)
            subj_mean_stage_lin = np.nanmean(psds_lin[mask, :, :], axis=0)
            subj_mean_stage_log = np.nanmean(psds_log[mask, :, :], axis=0)

            # on stocke pour chaque canal
            for ch_idx in range(n_channels):
                group_data_log[stage][ch_idx].append(subj_mean_stage_log[ch_idx, :])
                group_data_lin[stage][ch_idx].append(subj_mean_stage_lin[ch_idx, :])

    # ------------------------------------------------------------------
    # 2) Plot group-level : une figure, 3 subplots (1 par canal)
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, n_channels, figsize=(5 * n_channels, 5), sharey=True)
    fig1, axes1 = plt.subplots(1, n_channels, figsize=(5 * n_channels, 5), sharey=True)

    if n_channels == 1:
        axes = [axes]
        axes1 = [axes1]

    for ch_idx, ch_name in enumerate(ch_names):
        ax = axes[ch_idx]
        ax1 = axes1[ch_idx]

        for ss_idx, stage in enumerate(ss_order):
            ss_color = ss_colors[stage]
            subj_arrays_log = group_data_log[stage][ch_idx]
            subj_arrays_lin = group_data_lin[stage][ch_idx]

            if len(subj_arrays_log) == 0:
                # Aucun sujet n'a de données pour ce stade et ce canal en log
                continue
            if len(subj_arrays_lin) == 0:
                # Aucun sujet n'a de données pour ce stade et ce canal en linéaire
                continue

            # empilement: shape (n_subjects_with_stage, n_freqs)
            data_log = np.vstack(subj_arrays_log)
            data_lin = np.vstack(subj_arrays_lin)

            # # moyenne & SEM en espace linéaire (si espace linéaire utilisé)
            mean_lin = np.nanmean(data_lin, axis=0)
            sem_lin = np.nanstd(data_lin, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(data_lin), axis=0))

            lower_lin = np.maximum(mean_lin - sem_lin, np.finfo(float).eps) # avoid values <= 0 (log does not like them)
            upper_lin = mean_lin + sem_lin

            # # passage en dB
            mean_db_lin = 10 * np.log10(mean_lin)
            lower_db_lin = 10 * np.log10(lower_lin)
            upper_db_lin = 10 * np.log10(upper_lin)
            
            # moyenne & SEM en espace log 
            # mean_db = data_log.mean(axis=0)
            # sem_db = data_log.std(axis=0, ddof=1) / np.sqrt(data_log.shape[0])
            mean_db = np.nanmean(data_log, axis=0)
            sem_db  = np.nanstd(data_log, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(data_log), axis=0))
            lower_db = mean_db - sem_db
            upper_db = mean_db + sem_db
            
            # plot in log
            ax.plot(
                freqs,
                mean_db,
                label=f'{stage} (n={data_log.shape[0]})',
                color=ss_color
            )

            ax.fill_between(
                freqs,
                lower_db,
                upper_db,
                alpha=0.3,
                color=ss_color
            )
            
            ax.set_title(f'{ch_name} - Group mean', fontsize=13)
            ax.set_xlabel('Frequency (Hz)', fontsize=12)
            if ch_idx == 0:
                ax.set_ylabel('PSD (dB/Hz)', fontsize=12)
            ax.set_xlim(freqs[0], freqs[-1])
            ax.tick_params(axis='both', labelsize=10)
            ax.legend(fontsize=9)
            
            # plot in lin
            ax1.plot(
                freqs,
                mean_db_lin,
                label=f'{stage} (n={data_lin.shape[0]})',
                color=ss_color
            )

            ax1.fill_between(
                freqs,
                lower_db_lin,
                upper_db_lin,
                alpha=0.3,
                color=ss_color
            )

            ax1.set_title(f'{ch_name} - Group mean', fontsize=13)
            ax1.set_xlabel('Frequency (Hz)', fontsize=12)
            if ch_idx == 0:
                ax1.set_ylabel('PSD (dB/Hz)', fontsize=12)
            ax1.set_xlim(freqs[0], freqs[-1])
            ax1.tick_params(axis='both', labelsize=10)
            ax1.legend(fontsize=9)

    fig.suptitle("Group-level PSD per sleep stage (average in log space)", fontsize=14)
    fig.tight_layout()
    
    fig1.suptitle("Group-level PSD per sleep stage (average in linear space)", fontsize=14)
    fig1.tight_layout()

    # Sauvegarde de la figure de groupe
    group_plot_path_log = os.path.join(reports_dir, "group_PSD_per_stage_fromlog.png")
    fig.savefig(group_plot_path_log, dpi=300)
    # plt.close(fig)
    group_plot_path_lin = os.path.join(reports_dir, "group_PSD_per_stage_fromlin.png")
    fig1.savefig(group_plot_path_lin, dpi=300)

    print(f"✅ Figure groupe sauvegardée : {group_plot_path_log}")
    print(f"✅ Figure groupe sauvegardée : {group_plot_path_lin}")

# ------------------------- plot at the group level ---------------------------
# ---------------- PSD after removing aperiodic comp. -------------------------
if len(ap_removed_psds_list) == 0:
    print("⚠️ No PSD in psds_list, can't plot at the group level.")
else:
    # All PSDs should have the same frequencies and channels
    freqs = ap_removed_psds_list[0].freqs
    ch_names = ap_removed_psds_list[0].ch_names
    n_channels = len(ch_names)

    # optional check : same freqs / same ch_names for all
    for psds in ap_removed_psds_list[1:]:
        if not np.allclose(freqs, psds.freqs):
            raise ValueError("Frequencies between participants are not consistent, please harmonize.")
        if psds.ch_names != ch_names:
            raise ValueError("Channel names between participants are not consistent, please harmonize.")

    # On garde seulement les stades effectivement présents dans au moins un psds.event_id
    # (par exemple si R n'est pas scoré dans certains sujets)
    all_keys = set().union(*[set(p.event_id.keys()) for p in ap_removed_psds_list])
    ss_order = [s for s in ss_order if s in all_keys]

    # ------------------------------------------------------------------
    # 1) On construit une structure group_data[stage][ch_idx]
    #    qui contient une liste d'array (n_freqs) = moyenne par SUJET
    #    des PSD (linéaires) pour ce stade et ce canal.
    
    # we can probably use mne.grand_average(psds_list) here 
    
    # ------------------------------------------------------------------
    group_data_ap_removed_log = {
        stage: {ch_idx: [] for ch_idx in range(n_channels)}
        for stage in ss_order
    }

    for psds in ap_removed_psds_list:
        # PSD are already in dB (with the aperiodic component removed)
        psds_lin = psds.get_data()  # shape: (n_epochs, n_channels, n_freqs)
        events = psds.events[:, 2]         # codes de stades (entiers)
        event_id = psds.event_id           # mapping nom -> code

        for stage in ss_order:
            if stage not in event_id:
                continue  # ce sujet n'a pas ce stade du tout dans ses event_id

            ss_code = event_id[stage]
            mask = (events == ss_code)

            if not np.any(mask):
                # ce sujet n'a pas d'epoch dans ce stade
                continue

            # moyenne par stade ET par canal (on garde la dimension freq)
            # résultat: (n_channels, n_freqs)
            subj_mean_stage_lin = np.nanmean(psds_lin[mask, :, :], axis=0)

            # on stocke pour chaque canal
            for ch_idx in range(n_channels):
                group_data_ap_removed_log[stage][ch_idx].append(subj_mean_stage_lin[ch_idx, :])

    # ------------------------------------------------------------------
    # 2) Plot group-level : une figure, 3 subplots (1 par canal)
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, n_channels, figsize=(5 * n_channels, 5), sharey=True)

    if n_channels == 1:
        axes = [axes]

    for ch_idx, ch_name in enumerate(ch_names):
        ax = axes[ch_idx]

        for ss_idx, stage in enumerate(ss_order):
            ss_color = ss_colors[stage]
            subj_arrays = group_data_ap_removed_log[stage][ch_idx]

            if len(subj_arrays) == 0:
                # Aucun sujet n'a de données pour ce stade et ce canal en linéaire
                continue

            # empilement: shape (n_subjects_with_stage, n_freqs)
            data_lin = np.vstack(subj_arrays)

            # # moyenne & SEM en espace linéaire (si espace linéaire utilisé)
            mean_lin = data_lin.mean(axis=0)
            sem_lin = data_lin.std(axis=0, ddof=1) / np.sqrt(data_lin.shape[0])

            lower_lin = np.maximum(mean_lin - sem_lin, np.finfo(float).eps) # avoid values <= 0 (log does not like them)
            upper_lin = mean_lin + sem_lin

            # plot in log
            ax.plot(
                freqs,
                mean_lin,
                label=f'{stage} (n={data_lin.shape[0]})',
                color=ss_color
            )

            ax.fill_between(
                freqs,
                lower_lin,
                upper_lin,
                alpha=0.3,
                color=ss_color
            )
            
            ax.set_title(f'{ch_name} - Group mean', fontsize=13)
            ax.set_xlabel('Frequency (Hz)', fontsize=12)
            if ch_idx == 0:
                ax.set_ylabel('PSD (dB/Hz)', fontsize=12)
            ax.set_xlim(freqs[0], freqs[-1])
            ax.tick_params(axis='both', labelsize=10)
            ax.legend(fontsize=9)
            

    fig.suptitle("Group-level PSD after removing aperiodic component per sleep stage (average in log space)", fontsize=14)
    fig.tight_layout()

    # Sauvegarde de la figure de groupe
    group_plot_path_aplog = os.path.join(reports_dir, "group_PSD_ap_removed_per_stage_fromlog.png")
    fig.savefig(group_plot_path_aplog, dpi=300)
    # plt.close(fig)
   
    print(f"✅ Figure groupe sauvegardée : {group_plot_path_aplog}")

#______________________________________________________________________________
# looking for difference between mean from linear and mean from log (post-check because group plot are weird)
for pi, psds in enumerate(psds_list):
    # PSD en linéaire (µV²/Hz)
    psds_lin = psds.get_data() * 1e12  # shape: (n_epochs, n_channels, n_freqs)
    psds_lin = np.where(psds_lin > 0, psds_lin, np.nan) # put nan where psds <= 0
    # psds_log = 10*np.log10(np.maximum(psds_lin.copy(), eps))
    psds_log = 10*np.log10(psds_lin)
    events = psds.events[:, 2]         # codes de stades (entiers)
    event_id = psds.event_id           # mapping nom -> code
    
    for stage in ss_order:
        if stage not in event_id:
            continue  # ce sujet n'a pas ce stade du tout dans ses event_id

        ss_code = event_id[stage]
        mask = (events == ss_code)

        if not np.any(mask):
            # ce sujet n'a pas d'epoch dans ce stade
            continue

        # moyenne par stade ET par canal (on garde la dimension freq)
        # résultat: (n_channels, n_freqs)
        subj_mean_stage_lin = np.nanmean(psds_lin[mask, :, :], axis=0)
        # subj_mean_stage_log = psds_log[mask, :, :].mean(axis=0)
        subj_mean_stage_log = np.nanmean(psds_log[mask, :, :], axis=0)
        
        # quick comparison
        # au niveau epoch (pour un sujet au canal frontal)
        P = psds_lin[mask, 0, :]  # epochs x freqs
        
        m1 = 10*np.log10(np.nanmean(P, axis=0))      # log(mean_lin)
        m2 = np.nanmean(10*np.log10(P), axis=0)      # mean_log
        
        diff = m1 - m2
        imax = np.nanargmax(diff)
        print(f"p{pi}, {stage}: max diff (dB) = {np.round(diff[imax], 3)} at freq {freqs[imax]}, std = {np.round(np.nanstd(diff), 3)}")
        
# quick check of flat signal
thr = 1e-20  # à ajuster
stage = "W"
for pi, psds in enumerate(psds_list):
    psds_lin = psds.get_data() * 1e12
    events = psds.events[:, 2]
    code = psds.event_id.get(stage, None)
    if code is None:
        continue
    mask = events == code

    P = psds_lin[mask, :, :]  # epochs x ch x freqs
    frac_tiny = np.mean(P < thr)
    print(f"p{pi}: frac(P<{thr}) in W = {frac_tiny:.4f}")
    
    # participant with an anormal fraction: 64, 60, 52, 47, 31, 21, 20
