import os
from os.path import join as ospj
import json
import glob
import numpy as np
import pandas as pd
import sys
import nibabel as nib
from statsmodels.tsa.ar_model import AutoReg
from scipy import signal
import pingouin as pg
import matplotlib.pyplot as plt
import math
import re

def main():

    #### Set dynamic connectivity parameters
    dyn_conn_len = 60 # s
    win_len = 10 # s
    win_overlap = 7.5 # s

    num_windows = num_wins(dyn_conn_len, win_len, win_overlap)
    win_displacement = win_len - win_overlap

    # Get current directory
    code_path = os.getcwd()
    
    # Get paths from config file
    with open(ospj(code_path, "config.json"), "rb") as f:
        config = json.load(f)
    repo_path = config["repoPath"]
    recon_path = config["reconPath"]
    data_path = ospj(repo_path, "outputs", "ieeg")
    
    sub_dirs = glob.glob(ospj(data_path, "sub-*"))
    for sub_dir in sub_dirs:
        
        sub = str(os.path.basename(sub_dir))

        print('-----------------------------------')
        print('Subject: ' + sub)
        print('-----------------------------------')

        # For debugging
        #if sub not in ['sub-RID0139', 'sub-RID0658', 'sub-RID0679']:
        #    continue

        # Get subdirectory paths
        ictal_dir = ospj(sub_dir,"raw","ictal")
        interictal_dir = ospj(sub_dir,"raw","interictal")
        preictal_dir = ospj(sub_dir,"raw","preictal")
    
        # Create output paths
        out_paths = {'ictal': ospj(data_path, sub, "dynamic_connectivity", "ictal"),
                     'preictal': ospj(data_path, sub, "dynamic_connectivity", "preictal"),
                     'interictal': ospj(data_path, sub, "dynamic_connectivity", "interictal"),
                     'coords': ospj(data_path, sub, "coords")}
        for out_path in out_paths.values():
            try:
                os.makedirs(out_path)
            except FileExistsError:
                pass

        # Get iEEG file lists
        ieeg_files = {'ictal': glob.glob(ospj(ictal_dir, "*.csv")),
                      'preictal': glob.glob(ospj(preictal_dir, "*preictal_time-matched*.csv")),
                      'interictal': glob.glob(ospj(interictal_dir, "*interictal_time-matched*.csv"))}
    
        # Extract ieeg_recon electrode labels
        recon_sub_path = ospj(recon_path, sub, 'derivatives', 'ieeg_recon', 'module3')
        recon_sub_file = glob.glob(ospj(recon_sub_path, "*.csv"))[0]
        recon_sub = pd.read_csv(recon_sub_file)
        recon_sub_elect_labels = set(recon_sub['name'])
    
        # Extract ieeg.org electrode labels from a test file
        sub_test_file = pd.read_csv(ieeg_files['ictal'][0])
        ieeg_elect_labels_raw = set(sub_test_file.columns)

        # Save a list of the beginnings of labels needing unique conversions (should be of the form GRIDxx)
        grid_label_prefixes = ["GT", "GI", "GM", "GS"]

        # Clean electrode labels (sub-RID0520 is a special case - see below)
        label_conv = {}
        ieeg_elect_labels = ieeg_elect_labels_raw.copy()
        for ieeg_label in ieeg_elect_labels_raw:
            
            # Save original electrode label
            label_conv[ieeg_label] = ieeg_label

            # Do unique conversion for older datasets with labels of the form ("EEG RPH 07-Ref")
            if ieeg_label.startswith('EEG'):

                new_label = ieeg_label.lstrip('EEG ')
                new_label = new_label.rstrip('-Ref')
                new_label = new_label.replace(" ","")

                lab = parse_electrode_labels(new_label)

                # Check for double digit number
                if len(lab[1]) == 2:

                    # Check for leading zero 
                    if lab[1][0] == str(0):

                        new_label = str(lab[0] + lab[1][1])

                        # Save label mapping
                        label_conv[ieeg_label] = new_label

                        # Save electrode label
                        ieeg_elect_labels.remove(ieeg_label)
                        ieeg_elect_labels.add(new_label)

            
            # Do unique conversion where electrodes should have "GRID" prefix according to ieeg_recon
            elif any(ieeg_label.startswith(substring) for substring in grid_label_prefixes):

                new_label = 'GRID' + ieeg_label[2:]

                # Save label mapping
                label_conv[ieeg_label] = new_label

                # Save electrode label
                ieeg_elect_labels.remove(ieeg_label)
                ieeg_elect_labels.add(new_label)


            else:
            
                # Parse label
                lab = parse_electrode_labels(ieeg_label)
    
                if sub != "sub-RID0520":

                    # Check for double digit number
                    if len(lab[1]) == 2:

                        # Check for leading zero 
                        if lab[1][0] == str(0):

                            new_label = str(lab[0] + lab[1][1])

                            # Save label mapping
                            label_conv[ieeg_label] = new_label

                            # Save electrode label
                            ieeg_elect_labels.remove(ieeg_label)
                            ieeg_elect_labels.add(new_label)

                # Do unique conversion for sub-RID0520's LG electrodes (labels have the form "LGr40" in ieeg_recon outputs)
                else:

                    # Check for LG label
                    if lab[0] == 'LG':

                        if len(lab[1]) == 2:

                            # Check for leading zero
                            if lab[1][0] == str(0):

                                new_label = str(lab[0] + 'r' + lab[1][1])

                            else:

                                new_label = str(lab[0] + 'r' + lab[1])

                            # Save label mapping
                            label_conv[ieeg_label] = new_label

                            # Replace electrode label
                            ieeg_elect_labels.remove(ieeg_label)
                            ieeg_elect_labels.add(new_label)

                        else:

                            new_label = str(lab[0] + 'r' + lab[1])

                            # Save label mapping
                            label_conv[ieeg_label] = new_label

                            # Replace electrode label
                            ieeg_elect_labels.remove(ieeg_label)
                            ieeg_elect_labels.add(new_label)

                    elif len(lab[1]) == 2:

                        # Check for leading zero
                        if lab[1][0] == str(0):

                            new_label = str(lab[0] + lab[1][1])

                            # Save label mapping
                            label_conv[ieeg_label] = new_label

                            # Replace electrode label
                            ieeg_elect_labels.remove(ieeg_label)
                            ieeg_elect_labels.add(new_label)

        ieeg_elect_labels_wm = []
        ieeg_elect_labels_gm = []
        elect_coords_wm = pd.DataFrame(columns = ['label', 'x', 'y', 'z'])
        elect_coords_gm = pd.DataFrame(columns = ['label', 'x', 'y', 'z'])
        for label in ieeg_elect_labels:

            # Find ieeg_recon outputs corresponding to label
            recon_elect_atlas_label = recon_sub.loc[recon_sub['name'] == label,'label'].to_string(index=False)

            if recon_elect_atlas_label  == 'Series([], )':
                continue

            # Filter by atlas labels with "White-Matter" tag
            if "White-Matter" in recon_elect_atlas_label:

                # Save electrode label to list
                ieeg_elect_labels_wm.append(recon_sub.loc[recon_sub['name'] == label,'name'].to_string(index=False))

                # Save electrode labels + coordinates
                x = recon_sub.loc[recon_sub['name'] == label,'x'].to_string(index=False)
                y = recon_sub.loc[recon_sub['name'] == label,'y'].to_string(index=False)
                z = recon_sub.loc[recon_sub['name'] == label,'z'].to_string(index=False)

                elect_coords_wm.loc[len(elect_coords_wm)] = [label, x, y, z]

            elif recon_elect_atlas_label != "Empty-Label":

                # Save electrode label to list
                ieeg_elect_labels_gm.append(recon_sub.loc[recon_sub['name'] == label,'name'].to_string(index=False))

                # Save electrode labels + coordinates
                x = recon_sub.loc[recon_sub['name'] == label,'x'].to_string(index=False)
                y = recon_sub.loc[recon_sub['name'] == label,'y'].to_string(index=False)
                z = recon_sub.loc[recon_sub['name'] == label,'z'].to_string(index=False)

                elect_coords_gm.loc[len(elect_coords_gm)] = [label, x, y, z]


        #print("White matter electrodes:")
        #print(elect_coords_wm)
        #print('============================================')
        #print("Grey matter electrodes:")
        #print(elect_coords_gm)
        #print('============================================')

        # Save coordinates
        elect_coords_wm.to_csv(out_paths['coords'] + '/' + sub + '_wm_electrode_coords.csv', index=False, header=True)
        elect_coords_gm.to_csv(out_paths['coords'] + '/' + sub + '_gm_electrode_coords.csv', index=False, header=True)

        # Concatenate
        elect_coords_wm['tissue'] = 'WM'
        elect_coords_gm['tissue'] = 'GM'
        elect_coords = pd.concat([elect_coords_wm, elect_coords_gm])

        # Compute connectivity matrices
        periods = ['ictal', 'preictal', 'interictal']
        for period in periods:

            print('--Computing ' + period  + ' connectivity matrices...')

            for file in ieeg_files[period]:

                if period != 'interictal':
                    sz_id_pattern = r'ieeg_id-\d+_sz-\d+'
                    sz_id_match = re.search(sz_id_pattern, file)
                    sz_id = sz_id_match.group()
                else:
                    sz_id_pattern = r'ieeg_id-\d+_window-\d+'
                    sz_id_match = re.search(sz_id_pattern, file)
                    sz_id = sz_id_match.group()
                    sz_id = sz_id.replace('window', 'sz')

                print('Seizure ID: ' + str(sz_id))

                

                # Make output directory
                out_path_win_param = out_paths[period] + '/' + str(sz_id) + '/dyn_conn_len-' + str(int(dyn_conn_len*1000)) + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000))
                if os.path.exists(out_path_win_param):
                    print('----Skipping: Connectivity matrices already made')
                    continue
                else:
                    os.makedirs(out_path_win_param)

                # Load file
                f = pd.read_csv(file)

                if np.shape(f)[0] < int(round(250*dyn_conn_len,0)):
                    print('Time series too short - only has ' + str(np.shape(f)[0]) + '/' + str(int(round(250*dyn_conn_len,0)))  + ' necessary samples')
                    continue

                # Drop inconsistent electrodes across files for sub-RID0490
                if sub == "sub-RID0490":
                    if 'ieeg_id-1' in file:
                        f = f.drop('RE06', axis=1)
                    elif 'ieeg_id-2' in file:
                        f = f.drop('RE07', axis=1)

                # Extract ieeg.org electrode labels
                ieeg_elect_labels = list(f.columns)

                # Find white matter ieeg.org electrode labels
                f_wm = f.copy()
                for label in ieeg_elect_labels:

                    # Extract white matter labels
                    if label_conv[label] in ieeg_elect_labels_wm:

                        # Rename electrode label to match ieeg_recon outputs
                        f_wm.rename(columns={label: label_conv[label]}, inplace=True) 

                    else:

                        # Drop electrodes
                        f_wm = f_wm.drop(label, axis=1)

                # Find grey matter ieeg.org electrode labels
                f_gm = f.copy()
                for label in ieeg_elect_labels:

                    # Extract grey matter labels
                    if label_conv[label] in ieeg_elect_labels_gm:

                        # Rename electrode label to match ieeg_recon outputs
                        f_gm.rename(columns={label: label_conv[label]}, inplace=True)

                    else:
               
                        # Drop electrodes
                        f_gm = f_gm.drop(label, axis=1)

                ### Preprocess iEEG signals

                # Pre-whiten
                wm_ieeg = f_wm.apply(prewhiten_ar, axis=0)
                gm_ieeg = f_gm.apply(prewhiten_ar, axis=0)

                # Apply notch filter
                wm_ieeg = wm_ieeg.apply(notch_line, f0=60, axis=0)
                gm_ieeg = gm_ieeg.apply(notch_line, f0=60, axis=0)

                # Extract mean signal
                ieeg_concat = pd.concat([wm_ieeg, gm_ieeg], axis=1)
                ieeg_avg_signal = ieeg_concat.mean(axis=1)

                # Use common average referencing (CAR) montage
                wm_ieeg = wm_ieeg.apply(car_montage, avg_signal=ieeg_avg_signal, axis=0)
                gm_ieeg = gm_ieeg.apply(car_montage, avg_signal=ieeg_avg_signal, axis=0)

                # Include only the last 'dyn_conn_len' (default: 10s) of the signal
                wm_ieeg = wm_ieeg.iloc[-dyn_conn_len*250:]
                gm_ieeg = gm_ieeg.iloc[-dyn_conn_len*250:]

                ########### Save broadband connectivity matrices

                # Initialize window
                window_start = 0
                window_end = win_len

                for window in range(1,num_windows+1):

                    # Define signal window
                    sample_start = int(window_start*250)
                    sample_end = int(window_end*250)

                    # Extract iEEG data from window
                    wm_ieeg_window = wm_ieeg.iloc[sample_start:sample_end,:]
                    gm_ieeg_window = gm_ieeg.iloc[sample_start:sample_end,:]

                    # Create connectivity matrices
                    wm_conn_pearson = wm_ieeg_window.corr(method='pearson')
                    wm_conn_spearman = wm_ieeg_window.corr(method='spearman')
                    gm_conn_pearson = gm_ieeg_window.corr(method='pearson')
                    gm_conn_spearman = gm_ieeg_window.corr(method='spearman')

                    # Only keep upper triangle
                    wm_conn_pearson_upper = wm_conn_pearson.mask(np.triu(np.ones(wm_conn_pearson.shape, dtype=np.bool_)))
                    wm_conn_spearman_upper = wm_conn_spearman.mask(np.triu(np.ones(wm_conn_spearman.shape, dtype=np.bool_)))
                    gm_conn_pearson_upper = gm_conn_pearson.mask(np.triu(np.ones(gm_conn_pearson.shape, dtype=np.bool_)))
                    gm_conn_spearman_upper = gm_conn_spearman.mask(np.triu(np.ones(gm_conn_spearman.shape, dtype=np.bool_)))

                    # Save broadband connectivity matrix
                    wm_conn_pearson_upper.to_csv(out_path_win_param + '/' + sub + '_wm_connectivity-pearson_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_broadband.csv', index=True, header=True)
                    wm_conn_spearman_upper.to_csv(out_path_win_param + '/' + sub + '_wm_connectivity-spearman_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_broadband.csv', index=True, header=True)
                    gm_conn_pearson_upper.to_csv(out_path_win_param + '/' + sub + '_gm_connectivity-pearson_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_broadband.csv', index=True, header=True)
                    gm_conn_spearman_upper.to_csv(out_path_win_param + '/' + sub + '_gm_connectivity-spearman_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_broadband.csv', index=True, header=True)
            

                    ########### Save band-specific connectivity matrices
                    freq_bands = {'delta': [2,4], 'theta': [4,8], 'alpha': [8,12], 'beta': [12,30], 'gamma': [30,60]}
                    for band in freq_bands.keys():

                        wm_ieeg_bp = wm_ieeg_window.copy()
                        gm_ieeg_bp = gm_ieeg_window.copy()

                        w0 = freq_bands[band][0]
                        w1 = freq_bands[band][1]

                        # Apply bandpass filter
                        wm_ieeg_bp = wm_ieeg_bp.apply(apply_bandpass, w0=w0, w1=w1, axis=0)
                        gm_ieeg_bp = gm_ieeg_bp.apply(apply_bandpass, w0=w0, w1=w1, axis=0)

                        # Create correlation matrix
                        wm_bp_conn_pearson = wm_ieeg_bp.corr(method='pearson')
                        wm_bp_conn_spearman = wm_ieeg_bp.corr(method='spearman')
                        gm_bp_conn_pearson = gm_ieeg_bp.corr(method='pearson')
                        gm_bp_conn_spearman = gm_ieeg_bp.corr(method='spearman')

                        # Only keep upper triangle
                        wm_bp_conn_pearson_upper = wm_bp_conn_pearson.mask(np.triu(np.ones(wm_bp_conn_pearson.shape, dtype=np.bool_)))
                        wm_bp_conn_spearman_upper = wm_bp_conn_spearman.mask(np.triu(np.ones(wm_bp_conn_spearman.shape, dtype=np.bool_)))
                        gm_bp_conn_pearson_upper = gm_bp_conn_pearson.mask(np.triu(np.ones(gm_bp_conn_pearson.shape, dtype=np.bool_)))
                        gm_bp_conn_spearman_upper = gm_bp_conn_spearman.mask(np.triu(np.ones(gm_bp_conn_spearman.shape, dtype=np.bool_)))

                        # Save bandpass filtered connectivity matrix
                        wm_bp_conn_pearson_upper.to_csv(out_path_win_param + '/' + sub + '_wm_connectivity-pearson_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_' + band + '.csv', index=True, header=True)
                        wm_bp_conn_spearman_upper.to_csv(out_path_win_param + '/' + sub + '_wm_connectivity-spearman_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_' + band + '.csv', index=True, header=True)
                        gm_bp_conn_pearson_upper.to_csv(out_path_win_param + '/' + sub + '_gm_connectivity-pearson_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_' + band + '.csv', index=True, header=True)
                        gm_bp_conn_spearman_upper.to_csv(out_path_win_param + '/' + sub + '_gm_connectivity-spearman_' + period + '_winlen-' + str(int(win_len*1000)) + '_winoverlap-' + str(int(win_overlap*1000)) + '_win-' + str(window) + '_' + band + '.csv', index=True, header=True)

                    # Slide window
                    window_start = window*win_displacement
                    window_end = window*win_displacement + win_len

def num_wins(dyn_conn_len, win_len, win_overlap):

    win_displacement = win_len - win_overlap

    # Initialize window
    window_start = 0
    window_end = win_len

    n_windows = 0
    while (window_start + win_len) <= dyn_conn_len:

        # Update number of windows
        n_windows = n_windows + 1

        # Slide window
        window_start = round(window_start + win_displacement,3)
        window_end = round(window_end + win_displacement,3)

    return n_windows

# Parse electrode labels for cleaning
def parse_electrode_labels(s):
    
    head = s.rstrip('0123456789')
    tail = s[len(head):]

    return head, tail

# Pre-whiten signal through an AR(1) model
def prewhiten_ar(data):

    prewhitened_signal = AutoReg(data,lags=1).fit().resid

    return prewhitened_signal

# Apply bandpass filter within specific frequency bands
def apply_bandpass(data,w0,w1):
    
    order = 3
    nyq = 0.5 * 250
    low = w0 / nyq
    high = w1 / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    data_bp = data.copy()
    data_bp = signal.filtfilt(b,a,data_bp)

    return data_bp

# Apply notch filters at f0 and its harmonics
def notch_line(data, f0):

    # Quality factor
    Q = 30

    # Sampling rate
    fs = 250

    # Normalize
    nyq = 0.5 * fs
    w0 = f0 / nyq

    data_notch = data.copy()
    while w0 <= 1:
        b,a = signal.iirnotch(w0, Q=Q)
        data_notch = signal.filtfilt(b,a,data_notch)
        w0 = w0 + w0

    return data_notch

def compute_coherence_matrix(df, fs, freq_range):

    coherence_matrix = pd.DataFrame(index=df.columns, columns=df.columns)
    for i in range(len(coherence_matrix)):
        for j in range(i, len(coherence_matrix)):
            f, Cxy = signal.coherence(df.iloc[:,i], df.iloc[:,j], fs, nperseg=fs/2)
            coherence_matrix.iloc[i,j] = max(Cxy[(f>=freq_range[0]) & (f<=freq_range[1])])
            coherence_matrix.iloc[j,i] = coherence_matrix.iloc[i,j]
            
    return coherence_matrix

def compute_partial_corr_matrix(df, corr_method):

    parcorr_matrix = pd.DataFrame(index=df.columns, columns=df.columns)
    for i in range(len(parcorr_matrix)):
        for j in range(i, len(parcorr_matrix)):

            if i != j:
                covars = list(set(df.columns) - set([df.columns[i], df.columns[j]]))
                parcorr_matrix.iloc[i,j] = pg.partial_corr(data=df, x=df.columns[i], y=df.columns[j], covar=covars, method=corr_method)['r'][0]
                parcorr_matrix.iloc[j,i] = parcorr_matrix.iloc[i,j]

    return parcorr_matrix

def car_montage(data, avg_signal):
    data_car = data - avg_signal
    return data - avg_signal

def save_montage(data):

    
    # Set the figure size and resolution
    fig = plt.figure(figsize=(50, 20), dpi=100)

    # Set the number of rows and columns in the plot
    n_rows = int(math.ceil(data.shape[1]/3))
    n_cols = 3

    # Loop through the columns and plot each signal
    for i, col in enumerate(data.columns):
        ax = fig.add_subplot(n_rows, n_cols, i+1)
        ax.plot(data[col])
        ax.set_title(col, fontsize=8)
        ax.set_yticks([])
        ax.set_xticks([])

    # Adjust the spacing between subplots
    fig.subplots_adjust(hspace=1)

    return fig

if __name__ == "__main__":
    main()
