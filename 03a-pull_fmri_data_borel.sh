#!/bin/bash

ieeg_recon_dir=/mnt/leif/littlab/users/rsg20/ieeg_recon/BIDS
sublist="sub-RID0031
sub-RID0032
sub-RID0037
sub-RID0050
sub-RID0051
sub-RID0089
sub-RID0102
sub-RID0117
sub-RID0139
sub-RID0143
sub-RID0194
sub-RID0278
sub-RID0309
sub-RID0320
sub-RID0365
sub-RID0420
sub-RID0440
sub-RID0459
sub-RID0490
sub-RID0502
sub-RID0508
sub-RID0520
sub-RID0522
sub-RID0529
sub-RID0536
sub-RID0566
sub-RID0572
sub-RID0583
sub-RID0596
sub-RID0646
sub-RID0648
sub-RID0652
sub-RID0658
sub-RID0679"

for sub in ${sublist}; do

	echo ${sub}

	subdir=${ieeg_recon_dir}/${sub}
	if [ -d ${subdir} ]; then
	
		outdir=/mnt/leif/littlab/users/mjaskir/iEEG_fMRI_WM/source_data/fmri/${sub}
		if [ ! -d ${outdir} ]; then
			mkdir -p ${outdir}
		fi

		# Move T1w
		if [ -f ${subdir}/derivatives/ieeg_recon/module2/${sub}_ses-research3T_acq-3D_space-T00mri_T1w_ras.nii.gz ]; then
			cp ${subdir}/derivatives/ieeg_recon/module2/${sub}_ses-research3T_acq-3D_space-T00mri_T1w_ras.nii.gz ${outdir}
		elif [ -d ${subdir}/derivatives/ieeg_recon/module2/MRI_RAS  ]; then
			cp ${subdir}/derivatives/ieeg_recon/module2/MRI_RAS/${sub}_ses-research3T_acq-3D_space-T00mri_T1w.nii.gz ${outdir}
		else
			echo "${sub} has an unusual T1w file name or is missing module 2 outputs"
		fi

		# Move T1w with electrode spheres
		if [ -f ${subdir}/derivatives/ieeg_recon/module2/${sub}_ses-research3T_acq-3D_space-T00mri_T1w_ras_electrode_spheres.nii.gz ]; then
			cp ${subdir}/derivatives/ieeg_recon/module2/${sub}_ses-research3T_acq-3D_space-T00mri_T1w_ras_electrode_spheres.nii.gz ${outdir}
		elif [ -f ${subdir}/derivatives/ieeg_recon/module2/${sub}_ses-research3T_acq-3D_space-T00mri_T1w_electrode_spheres.nii.gz ]; then
			cp ${subdir}/derivatives/ieeg_recon/module2/${sub}_ses-research3T_acq-3D_space-T00mri_T1w_electrode_spheres.nii.gz ${outdir}
		else
			echo "${sub} has an unusual T1w + electrode spheres file name or is missing module 2 outputs"
		fi

	fi
done
