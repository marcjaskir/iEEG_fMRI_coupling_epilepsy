## Characterizing iEEG-fMRI coupling in epilepsy patients

### Description
Using resting state fMRI and interictal, preictal, and ictal iEEG data from 32 epilepsy patients, we evaluated differences in average iEEG-fMRI correlations and the associations between their graph theoretic properties (global efficiency) across different tissue types, seizure states, and frequency bands. We also used dynamic functional connectivity to elucidate seizure state-specific changes in iEEG-fMRI coupling over time.

### Faculty Lead
Kathryn Davis

### Dependencies
- Python API for IEEG.org (https://github.com/ieeg-portal/ieegpy)
    - After an IEEG.org account is created, you may use the function in software/generate_ieegpy_password_file.py to create a binary password file. Add this file to a new directory in this repository called "ieeg_login"
- config.json should specify absolute paths on your file system along with your IEEG.org credentials:
    - "repoPath" should specify this repository
    - "reconPath" should specify the path to outputs from iEEG-Recon (https://github.com/penn-cnt/ieeg-recon)
    - "iEEG_username" should specify your iEEG.org username
- Python scripts should be run from a conda environment (Python 3.8) using the following dependencies: software/conda/environment.yml 
- ANTS (https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS)
