## Characterizing iEEG-fMRI coupling in epilepsy patients

### Description


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
