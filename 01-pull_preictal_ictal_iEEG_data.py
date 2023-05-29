from ieeg.auth import Session
import os
from os.path import join as ospj
import sys
import numpy as np
import pandas as pd
import json
import random
import itertools
from ieeg.auth import Session
from scipy.io import loadmat
from scipy import signal
from tqdm import tqdm

def main():

    # Get current directory
    code_path = os.getcwd()
    
    # Get absolute path from config file
    with open(ospj(code_path, "config.json"), "rb") as f:
        config = json.load(f)
    repo_path = config["repoPath"]
    ieeg_username = config["iEEG_username"]
    data_path = ospj(repo_path, "source_data")
    metadata_path = ospj(data_path, "metadata")
    
    # Load subject list (contains both RIDs and HUP IDs) and convert to a dictionary
    subject_list = pd.read_csv(ospj(metadata_path, "subjects_w_fmri_ieeg.csv"))
    subject_dict = dict(zip(subject_list['hup'],subject_list['rid']))
    
    # Load metadata file
    metadata = pd.read_excel(
            ospj(metadata_path, "Manual validation.xlsx"), sheet_name="AllSeizureTimes"
    ).dropna(how="all")
    
    # Load sleep data and convert to data frame
    sleep_mat = loadmat(ospj(metadata_path, "sleeptimes.mat"))['sleepdata']
    sleep_df = pd.DataFrame({'Patient': np.squeeze(sleep_mat['name']), 
                             'IEEG_file_num': np.squeeze(sleep_mat['file']), 
                             'Times': np.squeeze(sleep_mat['t']), 
                             'Sleep': np.squeeze(sleep_mat['sleep'])})
    
    # Define length of preictal window
    preictal_window = 60  # 60 seconds
  
    # Start iEEG session
    with open(ospj(code_path, "ieeg_login", "ieeglogin.bin"), "r") as f:
        session = Session(ieeg_username, f.read())

    # Create time windows output directory
    metadata_out_path = ospj(repo_path, "outputs", "ieeg", "time_windows")
    try:
        os.makedirs(metadata_out_path)
    except FileExistsError:
        pass

    # Create file headers
    if not os.path.isfile(ospj(metadata_out_path, 'ictal_time_windows.csv')):
        with open(ospj(metadata_out_path, "ictal_time_windows.csv"),"a") as f:
            f.write("Patient,IEEGID,IEEGname,sz,start,end\n")
    if not os.path.isfile(ospj(metadata_out_path, 'preictal_time_windows.csv')):
        with open(ospj(metadata_out_path, "preictal_time_windows.csv"),"a") as f:
            f.write("Patient,IEEGID,IEEGname,sz,start,end\n")
    if not os.path.isfile(ospj(metadata_out_path, 'preictal_time_windows_matched.csv')):
        with open(ospj(metadata_out_path, "preictal_time_windows_matched.csv"),"a") as f:
            f.write("Patient,IEEGID,IEEGname,sz,start,end\n")

    ictal_times = pd.DataFrame(columns = ['Patient', 'IEEGID', 'IEEGname', 'sz', 'start', 'end'])
    preictal_times = pd.DataFrame(columns = ['Patient', 'IEEGID', 'IEEGname', 'sz', 'start', 'end'])
    for sub in subject_list['hup']:
    
        # Ignore NAs
        if pd.notnull(sub):

            print('-----------------------------------')
            print('Subject: ' + sub + ' =', subject_dict[sub])
            print('-----------------------------------')

            # Extract sleep data for subject
            sleep_sub = sleep_df[sleep_df["Patient"] == sub]
    
            # Check if sleep data is available for subject
            if len(list(sleep_sub["IEEG_file_num"])) == 0:
                print('Skipping - No sleep data available for ' + str(sub))
                continue

            # Skip subjects with missing iEEG recon module 3 outputs
            if sub in ["HUP191","HUP201"]:
                print('Skipping - No module 3 outputs')
                continue

            # Define raw signal output directories
            out_paths = {'ictal': ospj(repo_path, "outputs", "ieeg", subject_dict[sub], "raw", "ictal"),
                         'preictal': ospj(repo_path, "outputs", "ieeg", subject_dict[sub], "raw", "preictal")}

            # Create raw signal output directories
            for out_path in out_paths.values():
                try:
                    os.makedirs(out_path)
                except FileExistsError:
                    pass

            # Extract metadata for subject
            sub_metadata = metadata.loc[metadata['Patient'] == sub]
    
            ieeg_filename_list = pd.unique(sub_metadata['IEEGname'])
   
            # Iterate over iEEG files
            for filenames in range(len(ieeg_filename_list)):
    
                ieeg_filename = ieeg_filename_list[filenames]
    
                # For debugging
                #if 'HUP101' not in ieeg_filename:
                #    continue
    
                # Skip this dataset since its sampling rate was unusually low (256Hz)
                if 'HUP195_phaseII_D01' in ieeg_filename:
                    continue

                # Extract sleep data for corresponding iEEG file
                sleep_file = [item for sublist in list(sleep_sub["IEEG_file_num"])[0] for item in sublist]
                sleep_file_indices = [i for i, x in enumerate(list(sleep_file)) if x == (filenames+1)]
                sleep_status = [item for sublist in list(sleep_sub["Sleep"])[0] for item in sublist]
                sleep_status = [sleep_status[i] for i in sleep_file_indices]

                # Check if sleep data is available for iEEG file
                if (sum(sleep_status) == 0):
                    print('Skipping - No sleep data available for ' + str(ieeg_filename))
                    continue
    
                # Open iEEG dataset
                dataset = session.open_dataset(ieeg_filename)
    
                # Extract channel labels, indices, and time series details
                labels = dataset.get_channel_labels()
                indices = dataset.get_channel_indices(labels)
                ts_details = dataset.get_time_series_details(labels[0])
   
                # Extract metadata associated with iEEG file
                file_metadata = sub_metadata.loc[sub_metadata['IEEGname'] == ieeg_filename]

                # Fix incorrect IEEGIDs
                if 'HUP181_phaseII_D01' in ieeg_filename:
                    file_metadata.loc[file_metadata['IEEGID'] == 2.0, 'IEEGID'] = 1.0
                if 'HUP179_phaseII_D02' in ieeg_filename:
                    file_metadata.loc[file_metadata['IEEGID'] == 1.0, 'IEEGID'] = 2.0

                print('-----------------------------------')
                print('iEEG dataset name: ' + ieeg_filename)
                print('-----------------------------------')
    
                print('Saving ictal and preictal iEEG data')
    
                # Iterate over seizures
                num_sz = len(file_metadata)
                for sz in range(num_sz):
    
                    print('--Seizure: ' + str(sz+1) + '/' + str(num_sz))
   
                    # Check if iEEG data was already successfully pulled
                    if os.path.exists(out_paths['ictal'] + '/' + subject_dict[sub] + '_ictal_ieeg_id-' + str(filenames+1) + '_sz-' + str(sz+1)  + '.csv') and os.path.exists(out_paths['preictal'] + '/' + subject_dict[sub] + '_preictal_ieeg_id-' + str(filenames+1)  + '_sz-' + str(sz+1)  + '.csv'):
                        print('------Skipping: Data already pulled')
                        continue

                    # Skip HUP093_phaseII seizure #3 since preictal raw data is missing
                    if 'HUP093_phaseII' in ieeg_filename:
                        if sz == 2:
                            print('------Skipping: Missing raw preictal data')
                            continue

                    # Extract seizure metadata
                    ictal_metadata = file_metadata.iloc[sz,]

                    # Compute seizure duration
                    ictal_duration = float(ictal_metadata.loc['end']) - float(ictal_metadata.loc['start']) # seconds

                    # Skip seizures missing start/end times or with negative duration
                    if np.isnan(float(ictal_metadata.loc['end'])) or np.isnan(float(ictal_metadata.loc['start'])) or ictal_duration < 0:
                        print('----Skipping ' + ieeg_filename + ', Seizure #' + str(sz+1) + ' (Missing or nonsensible start/end time)')
                        continue

                    ### ICTAL ###

                    # Pull ictal iEEG data
                    ictal_ieeg = load_full_channels(dataset, ictal_metadata.loc['start'], ictal_duration, ts_details.sample_rate, indices)

                    # Downsample ictal iEEG data to 250Hz
                    ictal_ieeg_resamp = np.empty(shape=(int(250*ictal_duration), np.shape(ictal_ieeg)[1]))
                    for i in range(np.shape(ictal_ieeg)[1]):
                        ictal_ieeg_resamp[:,i] = signal.resample(ictal_ieeg[:,i], num=int(250*ictal_duration))

                    # Print warning if downsampling induced NAs
                    if np.sum(np.isnan(ictal_ieeg_resamp)) != 0:
                        print('----WARNING: NAs detected in downsampled ictal signal')
                        for i in range(0,np.shape(ictal_ieeg_resamp)[1]):
                            if np.sum(np.isnan(ictal_ieeg_resamp[:,i])) != 0:
                                print('------ ' + str(np.sum(np.isnan(ictal_ieeg_resamp))) + ' NAs detected after resampling ictal data for channel ' + str(np.array(labels)[i]))

                    # Add header
                    ictal_ieeg_resamp = np.concatenate([np.array(labels)[:,None].T, ictal_ieeg_resamp], axis=0)

                    # Save ictal iEEG data
                    np.savetxt(out_paths['ictal'] + '/' + subject_dict[sub] + '_ictal_ieeg_id-' + str(filenames+1) + '_sz-' + str(sz+1)  + '.csv', ictal_ieeg_resamp, delimiter=',', fmt='%s')

                    # Save ictal metadata
                    ictal_metadata_list = [subject_dict[ictal_metadata.loc['Patient']],
                                           str(ictal_metadata.loc['IEEGID']),
                                           ictal_metadata.loc['IEEGname'],
                                           str(sz+1),
                                           str(ictal_metadata.loc['start']),
                                           str(ictal_metadata.loc['end'])]
                    with open(ospj(metadata_out_path, "ictal_time_windows.csv"),"a") as f:
                        f.write(','.join(ictal_metadata_list) + "\n")

                    ### PREICTAL ###

                    # Pull preictal iEEG data
                    preictal_ieeg = load_full_channels(dataset, ictal_metadata.loc['start']-preictal_window, preictal_window, ts_details.sample_rate, indices)

                    # Downsample preictal iEEG data to 250Hz
                    preictal_ieeg_resamp = np.empty(shape=(int(250*preictal_window), np.shape(preictal_ieeg)[1]))
                    for i in range(np.shape(preictal_ieeg)[1]):
                        preictal_ieeg_resamp[:,i] = signal.resample(preictal_ieeg[:,i], num=int(250*preictal_window))

                    # Print warning if downsampling induced NAs
                    if np.sum(np.isnan(preictal_ieeg_resamp)) != 0:
                        print('----WARNING: NAs detected in downsampled preictal signal')
                        for i in range(np.shape(preictal_ieeg_resamp)[1]):
                            if np.sum(np.isnan(preictal_ieeg_resamp[:,i])) != 0:
                                print('------ ' + str(np.sum(np.isnan(preictal_ieeg_resamp))) + ' NAs detected after resampling preictal data for channel ' + str(np.array(labels)[i]))

                    # Add header
                    preictal_ieeg_resamp = np.concatenate([np.array(labels)[:,None].T, preictal_ieeg_resamp], axis=0)

                    # Save preictal iEEG data
                    np.savetxt(out_paths['preictal'] + '/' + subject_dict[sub] + '_preictal_ieeg_id-' + str(filenames+1) + '_sz-' + str(sz+1)  + '.csv', preictal_ieeg_resamp, delimiter=',', fmt='%s')

                    # Save preictal metadata
                    preictal_metadata_list = [subject_dict[ictal_metadata.loc['Patient']],
                                              str(ictal_metadata.loc['IEEGID']),
                                              ictal_metadata.loc['IEEGname'],
                                              str(sz+1),
                                              str(ictal_metadata.loc['start']-preictal_window),
                                              str(ictal_metadata.loc['start'])]
                    with open(ospj(metadata_out_path, "preictal_time_windows.csv"),"a") as f:
                        f.write(','.join(preictal_metadata_list) + "\n")

                    ### PREICTAL - DURATION MATCHED WITH ICTAL PERIOD ###

                    # Pull preictal iEEG data
                    preictal_ieeg_matched = load_full_channels(dataset, ictal_metadata.loc['start']-ictal_duration, ictal_duration, ts_details.sample_rate, indices)

                    # Downsample preictal iEEG data to 250Hz
                    preictal_ieeg_matched_resamp = np.empty(shape=(int(250*ictal_duration), np.shape(preictal_ieeg_matched)[1]))
                    for i in range(np.shape(preictal_ieeg_matched)[1]):
                        preictal_ieeg_matched_resamp[:,i] = signal.resample(preictal_ieeg_matched[:,i], num=int(250*ictal_duration))

                    # Print warning if downsampling induced NAs
                    if np.sum(np.isnan(preictal_ieeg_matched_resamp)) != 0:
                        print('----WARNING: NAs detected in downsampled preictal (time-matched) signal')
                        for i in range(np.shape(preictal_ieeg_matched_resamp)[1]):
                            if np.sum(np.isnan(preictal_ieeg_matched_resamp[:,i])) != 0:
                                print('------ ' + str(np.sum(np.isnan(preictal_ieeg_matched_resamp))) + ' NAs detected after resampling preictal (time-matched) data for channel ' + str(np.array(labels)[i]))

                    # Add header
                    preictal_ieeg_matched_resamp = np.concatenate([np.array(labels)[:,None].T, preictal_ieeg_matched_resamp], axis=0)

                    # Save preictal iEEG data
                    np.savetxt(out_paths['preictal'] + '/' + subject_dict[sub] + '_preictal_time-matched_ieeg_id-' + str(filenames+1) + '_sz-' + str(sz+1)  + '.csv', preictal_ieeg_matched_resamp, delimiter=',', fmt='%s')

                    # Save preictal metadata
                    preictal_metadata_matched_list = [subject_dict[ictal_metadata.loc['Patient']],
                                              str(ictal_metadata.loc['IEEGID']),
                                              ictal_metadata.loc['IEEGname'],
                                              str(sz+1),
                                              str(ictal_metadata.loc['start']-ictal_duration),
                                              str(ictal_metadata.loc['start'])]
                    with open(ospj(metadata_out_path, "preictal_time_windows_matched.csv"),"a") as f:
                        f.write(','.join(preictal_metadata_matched_list) + "\n")

def load_full_channels(dataset, start, duration_secs, sampling_rate, chn_idx):
    """
    Adapted from Brian's Spring 2023 BCI course, HW3

    Loads the entire channel from IEEG.org
    Input:
        dataset: the IEEG dataset object
        start: the start of the time window to load, in seconds
        duration_secs: the duration of the time window to load, in seconds
        sampling_rate: the sampling rate of the channel, in Hz
        chn_idx: the indicies of the m channels you want to load, as an array-like object
    Returns:
        [n, m] ndarry of the channels' values.
    """

    # Stores the segments of the channel's data
    chn_segments = []

    # How many segments do we expect?
    num_segments = int(np.ceil(duration_secs * sampling_rate / 6e4))

    # Segment start times and the step
    seg_starts, step = np.linspace(start*1e6, start*1e6 + duration_secs*1e6, num_segments, endpoint=False, retstep=True)

    # Get the segments
    for seg_start in seg_starts:
        chn_segments.append(dataset.get_data(seg_start, step, chn_idx))

    # Concatenate the segments vertically
    return np.vstack(chn_segments)

if __name__ == "__main__":
    main()
