import os
from os.path import join as ospj
import json
import glob
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

def main():

    # Get current directory
    code_path = os.getcwd()

    # Get absolute path from config file
    with open(ospj(code_path, "config.json"), "rb") as f:
        config = json.load(f)
    repo_path = config["repoPath"]
    data_path = ospj(repo_path, "outputs", "ieeg")

    sub_dirs = glob.glob(ospj(data_path, "sub-*"))
    for sub_dir in sub_dirs:

        sub = str(os.path.basename(sub_dir))

        print('-----------------------------------')
        print('Subject: ' + sub)
        print('-----------------------------------')

        if sub == "sub-RID0595":
            print('--Skipping - no module3 outputs')
            continue

        # Read in iEEG electrode coordinates
        coords_path = ospj(sub_dir, 'coords')
        coords_wm = pd.read_csv(ospj(coords_path, sub + '_wm_electrode_coords.csv'))
        coords_gm = pd.read_csv(ospj(coords_path, sub + '_gm_electrode_coords.csv'))

        # Read in iEEG connectivity matrix and extract label order
        ieeg_connectivity_wm = pd.read_csv(ospj(data_path, sub, 'connectivity', 'ictal', sub + '_wm_connectivity-pearson_ictal_alpha.csv'))
        ieeg_connectivity_gm = pd.read_csv(ospj(data_path, sub, 'connectivity', 'ictal', sub + '_gm_connectivity-pearson_ictal_alpha.csv'))
        ieeg_connectivity_wm_labels = ieeg_connectivity_wm.columns.tolist()[1:len(ieeg_connectivity_wm.columns.tolist())]
        ieeg_connectivity_gm_labels = ieeg_connectivity_gm.columns.tolist()[1:len(ieeg_connectivity_gm.columns.tolist())]

        # Calculate the pairwise distances between rows
        dist_array_wm = pdist(coords_wm[['x', 'y', 'z']])
        dist_array_gm = pdist(coords_gm[['x', 'y', 'z']])

        # Convert to square matrix
        dist_matrix_wm = squareform(dist_array_wm)
        dist_matrix_gm = squareform(dist_array_gm)

        # Convert to dataframe and add row and column labels
        distances_wm = pd.DataFrame(dist_matrix_wm, columns=coords_wm['label'], index=coords_wm['label'])
        distances_gm = pd.DataFrame(dist_matrix_gm, columns=coords_gm['label'], index=coords_gm['label'])

        # Reorder to match connectivity matrices
        distances_wm = distances_wm.loc[ieeg_connectivity_wm_labels, ieeg_connectivity_wm_labels]
        distances_gm = distances_gm.loc[ieeg_connectivity_gm_labels, ieeg_connectivity_gm_labels]

        # Only keep upper triangle
        distances_wm_upper = distances_wm.mask(np.triu(np.ones(distances_wm.shape, dtype=np.bool_)))
        distances_gm_upper = distances_gm.mask(np.triu(np.ones(distances_gm.shape, dtype=np.bool_)))

        # Save distance adjacency matrix
        distances_wm_upper.to_csv(ospj(coords_path, sub + '_wm_electrode_coords_euclidean.csv'))
        distances_gm_upper.to_csv(ospj(coords_path, sub + '_gm_electrode_coords_euclidean.csv'))

if __name__ == "__main__":
    main()
