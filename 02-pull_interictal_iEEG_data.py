from ieeg.auth import Session
import os
from os.path import join as ospj
import sys
import numpy as np
import pandas as pd
import json
import random
import itertools
import glob
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
    data_path = ospj(repo_path, "outputs", "ieeg")
    metadata_path = ospj(repo_path, "source_data", "metadata")

    # Load subject list (contains both RIDs and HUP IDs) and convert to a dictionary
    subject_list = pd.read_csv(ospj(metadata_path, "subjects_w_fmri_ieeg.csv"))
    subject_dict = dict(zip(subject_list['rid'],subject_list['hup']))

    # Load sleep data and convert to data frame
    sleep_mat = loadmat(ospj(metadata_path, "sleeptimes.mat"))['sleepdata']
    sleep_df = pd.DataFrame({'Patient': np.squeeze(sleep_mat['name']), 
                             'IEEG_file_num': np.squeeze(sleep_mat['file']), 
                             'Times': np.squeeze(sleep_mat['t']), 
                             'Sleep': np.squeeze(sleep_mat['sleep'])})

    # Load seizure times metadata and convert to data frame
    metadata = pd.read_csv(ospj(data_path, 'time_windows', 'ictal_time_windows.csv'))
    metadata_df = pd.DataFrame({'Patient': metadata['Patient'],
                                'IEEGID': metadata['IEEGID'],
                                'IEEGname': metadata['IEEGname'],
                                'sz': metadata['sz'],
                                'start': metadata['start'],
                                'end': metadata['end']})

    # Buffer between interictal period and perictal/ictal periods
    interictal_buffer = 6 * 60 * 60  # 6 hours

    # Initialize seed for reproducible interictal period selections
    random.seed(12345)

    # Start iEEG session
    with open(ospj(code_path, "ieeg_login", "ieeglogin.bin"), "r") as f:
        session = Session(ieeg_username, f.read())

    # Create time windows output directory
    metadata_out_path = ospj(data_path, "time_windows")
    try:
        os.makedirs(metadata_out_path)
    except FileExistsError:
        pass

    # Create file header
    if not os.path.isfile(ospj(metadata_out_path, "interictal_time_windows_matched.csv")):
        with open(ospj(metadata_out_path, "interictal_time_windows_matched.csv"),"a") as f:
            f.write("Patient,IEEGID,IEEGname,start,end\n")

    for sub_dir in glob.glob(ospj(data_path, 'sub-*')):
        
        sub = os.path.basename(sub_dir)

        # For debugging
        #if sub != 'sub-RID0529':
        #    continue

        print('-----------------------------------')
        print('Subject: ' + sub + ' =', subject_dict[sub])
        print('-----------------------------------')

        # Create raw signal output directory
        out_path = ospj(repo_path, "outputs", "ieeg", sub, "raw", "interictal")
        try:
            os.makedirs(out_path)
        except FileExistsError:
            pass

        # Extract sleep data for subject
        sleep_sub = sleep_df[sleep_df["Patient"] == subject_dict[sub]]

        # Extract metadata for subject's ictal iEEG data
        metadata_sub = metadata_df[metadata_df["Patient"] == sub]
        ieeg_file_names = np.unique(metadata_sub['IEEGname'])
        ieeg_files = np.unique(metadata_sub['IEEGID'])

        # Iterate over iEEG files
        for ieeg_file in ieeg_files:

            ieeg_fname = str(ieeg_file_names[np.where(ieeg_files == ieeg_file)][0])

            print('-----------------------------------')
            print('iEEG dataset name: ' + ieeg_fname)
            print('-----------------------------------')

            # Subset metadata_sub into metadata_sub_file
            metadata_sub_file = metadata_sub.loc[metadata_sub['IEEGID'] == ieeg_file]

            # Skip if data is already pulled
            if os.path.exists(out_path + '/' + sub + '_interictal_time-matched_ieeg_id-' + str(int(ieeg_file)) + '_window-' + str(len(metadata_sub_file))  + '.csv'):
                print('----Skipping: Data already pulled')
                continue

            # Open iEEG dataset
            dataset = session.open_dataset(ieeg_fname)

            # Extract channel labels, indices, and time series details
            labels = dataset.get_channel_labels()
            indices = dataset.get_channel_indices(labels)
            ts_details = dataset.get_time_series_details(labels[0])

            # Extract sleep data for corresponding iEEG file
            sleep_file = [item for sublist in list(sleep_sub["IEEG_file_num"])[0] for item in sublist]
            sleep_file_indices = [i for i, x in enumerate(list(sleep_file)) if x == ieeg_file]
            sleep_status = [item for sublist in list(sleep_sub["Sleep"])[0] for item in sublist]
            sleep_status = [sleep_status[i] for i in sleep_file_indices]
            sleep_time = [item for sublist in list(sleep_sub["Times"])[0] for item in sublist]
            sleep_time = [sleep_time[i] for i in sleep_file_indices]

            # Determine all awake time points
            awake_times = pd.DataFrame(columns = ['Patient', 'IEEGname', 'start', 'end', 'duration'])
            i = 0
            while i < (len(sleep_status)-1):

                # If participant is awake
                if sleep_status[i] == 0:

                    # Determine the extent of awake window
                    j = i + 1
                    while sleep_status[j] != 1 and j < (len(sleep_status)-1):
                        j = j + 1

                    # Save window
                    awake_times_df = pd.DataFrame({'Patient': [sub],
                                                    'IEEGname': [ieeg_file],
                                                    'start': [sleep_time[i]],
                                                    'end': [sleep_time[j-1]],
                                                    'duration': [(sleep_time[j-1] - sleep_time[i])]})
                    awake_times = pd.concat([awake_times, awake_times_df], axis=0)

                    i = j + 1

                else:

                    i = i + 1

            # Determine all time points within interictal buffers
            interictal_buffer_times = []
            for i in range(len(metadata_sub_file)):

                ictal_window = metadata_sub_file.iloc[i]
                interictal_buffer_window = np.arange(int(ictal_window.loc['start'])-interictal_buffer, int(ictal_window.loc['end'])+interictal_buffer + 1, 1)
                interictal_buffer_times = np.append(interictal_buffer_times, interictal_buffer_window)

            interictal_buffer_times = np.unique(interictal_buffer_times)

            # Determine possible interictal times
            possible_interictal_times = pd.DataFrame(columns = ['Patient', 'IEEGname', 'start', 'end', 'duration'])
            for i in range(len(awake_times)):

                awake_window = awake_times.iloc[i]
                awake_window_times = np.linspace(int(awake_window.loc['start']), int(awake_window.loc['end']), int(awake_window.loc['end'])-int(awake_window.loc['start']))

                # Check for overlap between awake window and interictal buffer times
                intersection = np.intersect1d(interictal_buffer_times, awake_window_times)
                if len(intersection) != 0:

                    # Case 1: Buffer overlaps with beginning (but not end) of awake window
                    if int(awake_window.loc['start']) == int(intersection[0]) and int(awake_window.loc['end']) != int(intersection[len(intersection)-1]):

                        possible_interictal_times_df = pd.DataFrame({'Patient': [sub],
                                                                    'IEEGname': [ieeg_file],
                                                                    'start': [intersection[len(intersection)-1]],
                                                                    'end': [awake_window.loc['end']],
                                                                    'duration': [awake_window.loc['end'] - intersection[len(intersection)-1]]})
                        possible_interictal_times = pd.concat([possible_interictal_times, possible_interictal_times_df], axis=0)

                    # Case 2: Buffer overlaps with end (but not beginning) of awake window
                    if int(awake_window.loc['start']) != int(intersection[0]) and int(awake_window.loc['end']) == int(intersection[len(intersection)-1]):

                        possible_interictal_times_df = pd.DataFrame({'Patient': [sub],
                                                                    'IEEGname': [ieeg_file],
                                                                    'start': [awake_window.loc['start']],
                                                                    'end': [intersection[0]],
                                                                    'duration': [intersection[0] - awake_window.loc['start']]})
                        possible_interictal_times = pd.concat([possible_interictal_times, possible_interictal_times_df], axis=0)

                    # Case 3: Buffer splits awake window
                    if int(awake_window.loc['start']) != int(intersection[0]) and int(awake_window.loc['end']) != int(intersection[len(intersection)-1]):

                        # First half
                        possible_interictal_times_df = pd.DataFrame({'Patient': [sub],
                                                                    'IEEGname': [ieeg_file],
                                                                    'start': [awake_window.loc['start']],
                                                                    'end': [intersection[0]],
                                                                    'duration': [intersection[0] - awake_window.loc['start']]})
                        possible_interictal_times = pd.concat([possible_interictal_times, possible_interictal_times_df], axis=0)

                        # Second half
                        possible_interictal_times_df = pd.DataFrame({'Patient': [sub],
                                                                    'IEEGname': [ieeg_file],
                                                                    'start': [intersection[len(intersection)-1]],
                                                                    'end': [awake_window.loc['end']],
                                                                    'duration': [awake_window.loc['end'] - intersection[len(intersection)-1]]})
                        possible_interictal_times = pd.concat([possible_interictal_times, possible_interictal_times_df], axis=0)

                    # Case 4: Entire awake window is inside of interictal buffer
                    if int(awake_window.loc['start']) == int(intersection[0]) and int(awake_window.loc['end']) == int(intersection[len(intersection)-1]):
                        continue

                # Entire awake window is outside of interictal buffer
                else:

                    possible_interictal_times_df = pd.DataFrame({'Patient': [sub],
                                                                'IEEGname': [ieeg_file],
                                                                'start': [awake_window.loc['start']],
                                                                'end': [awake_window.loc['end']],
                                                                'duration': [awake_window.loc['end'] - awake_window.loc['start']]})
                    possible_interictal_times = pd.concat([possible_interictal_times, possible_interictal_times_df], axis=0)

                i = i + 1

            # Sort possible interictal windows by duration
            possible_interictal_times_sorted = possible_interictal_times.sort_values('duration', ascending=False)

            # Remove short interictal windows (<120s)
            possible_interictal_times_sorted = possible_interictal_times_sorted.loc[possible_interictal_times_sorted['duration'] > 120,:]

            # For sub-RID0529, skip a possible interictal window with too many NAs
            if sub == 'sub-RID0529':
                possible_interictal_times_sorted = possible_interictal_times_sorted[possible_interictal_times_sorted['start'] != 516250]

            ### Determine interictal periods (time-matched with ictal duration) ###

            print('Determining interictal windows (time-matched)...')

            # Sort seizure metadata by duration
            metadata_sub_file['duration'] = metadata_sub_file['end'] - metadata_sub_file['start']
            metadata_sub_file = metadata_sub_file.sort_values(by='duration', ascending=False)

            interictal_times = pd.DataFrame(columns = ['Patient', 'IEEGID', 'IEEGname', 'start', 'end'])
            for i in tqdm(range(len(metadata_sub_file))):

                # Define ictal duration
                metadata_sub_file_sz = metadata_sub_file.iloc[i,:]
                ictal_duration = metadata_sub_file_sz.loc['end'] - metadata_sub_file_sz.loc['start']

                # If possible, pull each interictal period from a different possible interictal window
                if len(possible_interictal_times_sorted) >= len(metadata_sub_file):
                    num_interict = i

                # Otherwise, pull extra interictal periods from longest windows
                else:
                    num_interict = i % len(possible_interictal_times_sorted)

                # Define window start/end
                window = possible_interictal_times_sorted.iloc[num_interict,:]
                window_start = window.loc['start']
                window_end = window.loc['end']

                #print('----Search Window = Start: ' + str(window_start) + ' - End: ' + str(window_end))

                # Randomly select start of interictal period, ensuring no overlap with existing interictal periods
                searching = True
                while searching:

                    # Double check this line, since ictal_duration can be a float
                    interict_start = random.randint(int(window_start), int(window_end-ictal_duration))

                    num_overlaps = 0
                    for w in range(len(interictal_times)):

                        # Extract range of each existing interictal period
                        overlap_with = interictal_times.iloc[w,:]
                        overlap_with_start =  int(round(overlap_with['start'],0))
                        overlap_with_end = int(round(overlap_with['end'],0))
                        overlap_with_range = set(range(overlap_with_start, overlap_with_end))

                        # Check for overlap
                        if not overlap_with_range.intersection(range(interict_start, int(round(interict_start + ictal_duration,0)))) == set():
                            num_overlaps = num_overlaps + 1

                    # Load signal for all channels
                    available_signal = load_full_channels(dataset, interict_start, ictal_duration, ts_details.sample_rate, indices)

                    # Ensure no overlap with existing window or missing data
                    if num_overlaps == 0 and np.sum(np.isnan(available_signal)) == 0:
                        searching=False

                # Determine end of interictal period
                interict_end = interict_start + ictal_duration

                # Save interictal start and end times
                interictal_times_df = pd.DataFrame({'Patient': [sub],
                                                    'IEEGID': np.unique(metadata_sub_file['IEEGID']),
                                                    'IEEGname': np.unique(metadata_sub_file['IEEGname']),
                                                    'start': [interict_start],
                                                    'end': [interict_end]})
                interictal_times = pd.concat([interictal_times, interictal_times_df], axis = 0)

                #print('----Selected Window = Start: ' + str(interict_start) + ' - End: ' + str(interict_end))

            # Sort interictal times by start
            interictal_times = interictal_times.sort_values('start')

            # Pull interictal iEEG data
            print('Saving interictal windows (time-matched)...')
            for i in tqdm(range(len(interictal_times))):

                # Define ictal duration
                interictal_window = interictal_times.iloc[i,:]
                ictal_duration = interictal_window.loc['end'] - interictal_window.loc['start']

                window = interictal_times.iloc[i]
                window_start = window['start']

                # Pull
                interictal_ieeg = load_full_channels(dataset, window_start, ictal_duration, ts_details.sample_rate, indices)

                # Downsample iEEG data to 250Hz
                interictal_ieeg_resamp = np.empty(shape=(int(round(250*ictal_duration,0)), np.shape(interictal_ieeg)[1]))
                for j in range(np.shape(interictal_ieeg)[1]):
                    interictal_ieeg_resamp[:,j] = signal.resample(interictal_ieeg[:,j], num=int(round(250*ictal_duration,0)))

                # Print warning if downsampling induced NAs
                if np.sum(np.isnan(interictal_ieeg_resamp)) != 0:
                    print('----WARNING: NAs detected in downsampled interictal signal')
                    for j in range(np.shape(interictal_ieeg_resamp)[1]):
                        if np.sum(np.isnan(interictal_ieeg_resamp[:,j])) != 0:
                            print('------ ' + str(np.sum(np.isnan(interictal_ieeg_resamp))) + ' NAs detected after resampling interictal data for channel ' + str(np.array(labels)[j]))

                # Add header
                interictal_ieeg_resamp = np.concatenate([np.array(labels)[:,None].T, interictal_ieeg_resamp], axis=0)

                # Save interictal iEEG data
                np.savetxt(out_path + '/' + sub + '_interictal_time-matched_ieeg_id-' + str(int(np.unique(interictal_times['IEEGID'])[0])) + '_window-' + str(i+1)  + '.csv', interictal_ieeg_resamp, delimiter=',', fmt='%s')

            # Save interictal metadata
            interictal_times.to_csv(ospj(metadata_out_path, "interictal_time_windows_matched.csv"), mode='a', index=False, header=False)

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
    num_segments = int(np.ceil(duration_secs * sampling_rate / 6e2))

    # Segment start times and the step
    seg_starts, step = np.linspace(start*1e6, start*1e6 + duration_secs*1e6, num_segments, endpoint=False, retstep=True)

    # Get the segments
    for seg_start in seg_starts:
        chn_segments.append(dataset.get_data(seg_start, step, chn_idx))

    # Concatenate the segments vertically
    return np.vstack(chn_segments)

if __name__ == "__main__":
    main()
