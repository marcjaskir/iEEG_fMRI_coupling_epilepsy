#!/bin/bash

# NOTE: This script was executed from bscsub1.pmacs.upenn.edu:/home/mjaskir/iEEG_fMRI_WM/code to transfer fmriprep and XCPEngine derivatives to borel, but is copied here for reference

fmriprepdir=/project/davis_group/allucas/projects/ieeg_fmri_penn/outputs/derivatives/fmriprep
xcpdir=/project/davis_group/allucas/projects/ieeg_fmri_penn/outputs/derivatives/fc-36p_no_wm_regress_despike
sublist="sub-RID0051 sub-RID0102 sub-RID0143 sub-RID0089 sub-RID0309 sub-RID0365 sub-RID0440 sub-RID0490 sub-RID0508 sub-RID0420 sub-RID0522 sub-RID0520 sub-RID0536 sub-RID0595 sub-RID0646 sub-RID0572 sub-RID0648 sub-RID0566 sub-RID0031 sub-RID0117 sub-RID0278 sub-RID0139 sub-RID0320 sub-RID0454 sub-RID0194 sub-RID0476"

for sub in ${sublist}; do

	subdir=${fmriprepdir}/${sub}
	if [ -d ${subdir} ]; then
		for sesdir in ${subdir}/ses-*; do
			
			ses=$(basename ${sesdir})

			# Move white matter segmentation
			if [ -f ${sesdir}/anat/${sub}_${ses}_label-WM_probseg.nii.gz ]; then

				rsync ${sesdir}/anat/${sub}_${ses}_label-WM_probseg.nii.gz mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			elif [ -f ${sesdir}/anat/${sub}_${ses}_acq-3D_label-WM_probseg.nii.gz ]; then

				rsync ${sesdir}/anat/${sub}_${ses}_acq-3D_label-WM_probseg.nii.gz mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			fi

			# Move grey matter segmentation
			if [ -f ${sesdir}/anat/${sub}_${ses}_label-GM_probseg.nii.gz ]; then

				rsync ${sesdir}/anat/${sub}_${ses}_label-GM_probseg.nii.gz mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			elif [ -f ${sesdir}/anat/${sub}_${ses}_acq-3D_label-GM_probseg.nii.gz ]; then

				rsync ${sesdir}/anat/${sub}_${ses}_acq-3D_label-GM_probseg.nii.gz mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			fi

			# Move T1w --> EPI transformation
			if [ -f ${sesdir}/func/${sub}_${ses}_task-rest_from-T1w_to-scanner_mode-image_xfm.txt ]; then

				rsync ${sesdir}/func/${sub}_${ses}_task-rest_from-T1w_to-scanner_mode-image_xfm.txt mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			fi

                        # Move MNI --> T1w transformation
                        if [ -f ${sesdir}/anat/${sub}_${ses}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 ]; then

                                rsync ${sesdir}/anat/${sub}_${ses}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			elif [ -f ${sesdir}/anat/${sub}_${ses}_acq-3D_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 ]; then

				rsync ${sesdir}/anat/${sub}_${ses}_acq-3D_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

                        fi

			# Move BOLD reference image
			if [ -f ${sesdir}/func/${sub}_${ses}_task-rest_space-T1w_boldref.nii.gz ]; then

				rsync ${sesdir}/func/${sub}_${ses}_task-rest_space-T1w_boldref.nii.gz mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			fi

			# Move BOLD time series
			if [ -f ${xcpdir}/${sub}/${ses}/regress/${sub}_${ses}_residualised.nii.gz ]; then

				rsync ${xcpdir}/${sub}/${ses}/regress/${sub}_${ses}_residualised.nii.gz mjaskir@borel.seas.upenn.edu:/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}

			fi

		done
	fi

done
