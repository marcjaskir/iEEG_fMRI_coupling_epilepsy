import os
from os.path import join as ospj
import sys
import json
import glob
import numpy as np
import pandas as pd
import scipy
import nibabel as nib
import subprocess

def main():

    # Get current directory
    code_path = os.getcwd()

    # Get absolute path from config file
    with open(ospj(code_path, "config.json"), "rb") as f:
        config = json.load(f)
    repo_path = config["repoPath"]
    ieeg_path = ospj(repo_path, "outputs", "ieeg")
    fmri_path = ospj(repo_path, "source_data", "fmri")
    atlas_path = ospj(repo_path, "source_data", "atlases")

    # Determine spherical ROI radiuses
    radii = [8] # [6, 8, 10]

    # Determine whether we want to constrain spherical ROI masks by a white/grey matter segmentation
    mask_constrain = 1 # No: 0, Yes: 1
    
    sub_dirs = glob.glob(ospj(fmri_path, "sub-*"))
    for sub_dir in sub_dirs:

        sub = str(os.path.basename(sub_dir))

        print('-----------------------------------')
        print('Subject: ' + sub)
        print('-----------------------------------')

        if sub == "sub-RID0595":
            print('--Skipping - no module3 outputs')
            continue

        if sub == "sub-RID0102":
            print('--Skipping for now - issue with resampling coordinates')
            continue

        # Check that iEEG data is present
        if not os.path.isdir(ospj(ieeg_path, sub)):
            print('--Skipping - no iEEG data available')
            continue

        # Create output directory
        if mask_constrain == 0:
            out_paths = {'raw': ospj(repo_path, "outputs", "fmri", sub, "raw"),
                        'connectivity': ospj(repo_path, "outputs", "fmri", sub, "connectivity"),
                        'masks': ospj(repo_path, "outputs", "fmri", sub, "masks")}
        elif mask_constrain == 1:
            out_paths = {'raw': ospj(repo_path, "outputs", "fmri_seg_masked", sub, "raw"),
                        'connectivity': ospj(repo_path, "outputs", "fmri_seg_masked", sub, "connectivity"),
                        'masks': ospj(repo_path, "outputs", "fmri_seg_masked", sub, "masks")}
        for out_path in out_paths.values():
            try:
                os.makedirs(out_path)
            except FileExistsError:
                pass

        # Check if already run
        if os.path.isfile(ospj(out_paths['connectivity'], sub + '_gm_fmri_connectivity-spearman_radius-8mm.csv')):
            print('--Skipping - outputs already created')
            continue

        # Load T1w image
        if len(glob.glob(ospj(sub_dir,"*T1w.nii.gz"))) > 0:
            for t1_file in glob.glob(ospj(sub_dir,"*T1w.nii.gz")):
                t1 = nib.load(t1_file)
        elif len(glob.glob(ospj(sub_dir,"*T1w_ras.nii.gz"))) > 0:
            for t1_file in glob.glob(ospj(sub_dir,"*T1w_ras.nii.gz")):
                t1 = nib.load(t1_file)
        t1_data = t1.get_fdata()

        # Load BOLD image
        for bold_file in glob.glob(ospj(sub_dir,"*residualised.nii.gz")):
            bold = nib.load(bold_file)
        bold_data = bold.get_fdata()
        bold_header = bold.header
        bold_nvols = bold_header['dim'][4]

        # Load BOLD reference image, extracting its resolution
        for bold_ref_file in glob.glob(ospj(sub_dir,"*T1w_boldref.nii.gz")):
            bold_ref = nib.load(bold_ref_file)
        bold_ref_data = bold_ref.get_fdata()
        bold_ref_header = bold_ref.header
        bold_ref_dims = bold_ref_header.get_zooms()

        # Load T1 to BOLD transformation
        xfm_file_t1_bold = glob.glob(ospj(sub_dir,"*from-T1w_to-scanner_mode-image_xfm.txt"))[0]

        # Load white matter segmentation
        for wm_seg_file in glob.glob(ospj(sub_dir,"*WM_probseg.nii.gz")):
            wm_seg = nib.load(wm_seg_file)
        wm_seg_data = wm_seg.get_fdata()

        # Load grey matter segmentation
        for gm_seg_file in glob.glob(ospj(sub_dir,"*GM_probseg.nii.gz")):
            gm_seg = nib.load(gm_seg_file)
        gm_seg_data = gm_seg.get_fdata()

        # Resample segmentations into BOLD space, thresholding the resampled image at 75% and binarizing
        wm_seg_resampled_fname = ospj(out_paths['masks'], sub + '_wm_seg_resamp.nii.gz')
        subprocess.call(['antsApplyTransforms',
                    '-i', wm_seg_file,
                    '-o', wm_seg_resampled_fname,
                    '-r', bold_ref_file,
                    '-n', 'Gaussian',
                    '-t', xfm_file_t1_bold])
        wm_seg_resampled = nib.load(wm_seg_resampled_fname)
        wm_seg_resampled_data = wm_seg_resampled.get_fdata()
        wm_seg_resampled_thresh_bin = np.where(wm_seg_resampled_data > 0.75, 1, 0)
        wm_seg_resampled_thresh_bin_img = nib.Nifti1Image(wm_seg_resampled_thresh_bin, affine=bold_ref.affine, dtype='int64')
        wm_seg_resampled_thresh_bin_fname = ospj(out_paths['masks'], sub + '_wm_seg_resamp_thresh-75pct_bin.nii.gz')
        nib.save(wm_seg_resampled_thresh_bin_img, wm_seg_resampled_thresh_bin_fname)

        # For grey matter mask, binarize resampled image and substract white matter mask
        gm_seg_resampled_fname = ospj(out_paths['masks'], sub + '_gm_seg_resamp.nii.gz')
        subprocess.call(['antsApplyTransforms',
                    '-i', gm_seg_file,
                    '-o', gm_seg_resampled_fname,
                    '-r', bold_ref_file,
                    '-n', 'Gaussian',
                    '-t', xfm_file_t1_bold])
        gm_seg_resampled = nib.load(gm_seg_resampled_fname)
        gm_seg_resampled_data = gm_seg_resampled.get_fdata()
        gm_seg_resampled_bin = np.where(gm_seg_resampled_data > 0, 1, 0)
        gm_wm_diff = gm_seg_resampled_bin - wm_seg_resampled_thresh_bin
        gm_seg_resampled_noWM_bin = np.where(gm_wm_diff > 0, 1, 0)
        gm_seg_resampled_noWM_bin_img = nib.Nifti1Image(gm_seg_resampled_noWM_bin, affine=bold_ref.affine, dtype='int64')
        gm_seg_resampled_noWM_bin_fname = ospj(out_paths['masks'], sub + '_gm_seg_resamp_noWM_bin.nii.gz')
        nib.save(gm_seg_resampled_noWM_bin_img, gm_seg_resampled_noWM_bin_fname)

        # Load in an example iEEG connectivity matrix for matching label order
        ieeg_connectivity_wm = pd.read_csv(ospj(ieeg_path, sub, 'connectivity', 'ictal', sub + '_wm_connectivity-pearson_ictal_alpha.csv'))
        ieeg_connectivity_gm = pd.read_csv(ospj(ieeg_path, sub, 'connectivity', 'ictal', sub + '_gm_connectivity-pearson_ictal_alpha.csv'))
        ieeg_connectivity_wm_labels = ieeg_connectivity_wm.columns.tolist()[1:len(ieeg_connectivity_wm.columns.tolist())]
        ieeg_connectivity_gm_labels = ieeg_connectivity_gm.columns.tolist()[1:len(ieeg_connectivity_gm.columns.tolist())]

        # Load coordinates
        coords_dir = ospj(repo_path, "outputs", "ieeg", sub, "coords")
        coords_wm = pd.read_csv(ospj(coords_dir, sub + "_wm_electrode_coords.csv"))
        coords_gm = pd.read_csv(ospj(coords_dir, sub + "_gm_electrode_coords.csv"))

        # Create dictionaries to hold BOLD timeseries within spheres and 
        bold_timeseries_wm = {}
        bold_timeseries_gm = {}
        spherical_rois_reference_wm = {}
        spherical_rois_reference_gm = {}
        for radius in radii:
            bold_timeseries_wm[radius] = np.empty(shape=(bold_nvols, len(coords_wm)))
            bold_timeseries_gm[radius] = np.empty(shape=(bold_nvols, len(coords_gm)))
            spherical_rois_reference_wm[radius] = np.zeros(shape=np.shape(bold_ref_data))
            spherical_rois_reference_gm[radius] = np.zeros(shape=np.shape(bold_ref_data))

        print('--White matter...')

        labels_wm = []
        for coord in range(len(coords_wm)):

            # Extract coordinates
            coord_line = coords_wm.iloc[coord,:]
            label,x,y,z = (coord_line['label'], coord_line['x'], coord_line['y'], coord_line['z'])
            labels_wm.append(label)

            # Generate small sherical mask around coordinate
            coord_df = pd.DataFrame({'x': x, 'y': y, 'z': z, 't': 0}, index=[0])
            coord_mask = generate_sphere(np.zeros(shape=np.shape(t1_data)),
                                         int(round(coord_df['x'],0)),
                                         int(round(coord_df['y'],0)),
                                         int(round(coord_df['z'],0)),
                                         1,
                                         1,
                                         1,
                                         0,
                                         wm_seg_resampled_thresh_bin,
                                         't1')
            coord_mask_img = nib.Nifti1Image(coord_mask, affine=t1.affine)
            coord_mask_fname = ospj(out_paths['masks'], sub + '_mask_wm_electrode-' + str(coord) + '.nii.gz')
            nib.save(coord_mask_img, coord_mask_fname)

            # Resample mask into BOLD space, extracting its COG
            coord_mask_resampled_fname = ospj(out_paths['masks'], sub + '_mask_wm_electrode-' + str(coord) + '_resamp.nii.gz')
            subprocess.call(['antsApplyTransforms',
                        '-i', coord_mask_fname,
                        '-o', coord_mask_resampled_fname,
                        '-r', bold_ref_file,
                        '-n', 'Gaussian',
                        '-t', xfm_file_t1_bold])
            coord_mask_resampled = nib.load(coord_mask_resampled_fname)
            coord_mask_resampled_data = coord_mask_resampled.get_fdata()
            coord_resampled = scipy.ndimage.center_of_mass(coord_mask_resampled_data)
            x_resamp,y_resamp,z_resamp = [int(round(coord_resampled[0],0)), int(round(coord_resampled[1],0)), int(round(coord_resampled[2],0))]

            # Iterate over radii
            for radius in radii:

                # Generate spherical ROIs in BOLD space
                if bold_ref_dims == (2.0, 2.0, 3.0):
                    sphere = generate_sphere_anisotropic(np.zeros(shape=np.shape(bold_ref_data)),
                                                         x_resamp,
                                                         y_resamp,
                                                         z_resamp,
                                                         radius,
                                                         1,
                                                         mask_constrain,
                                                         wm_seg_resampled_thresh_bin,
                                                         'bold')

                elif bold_ref_dims == (2.0, 2.0, 2.0):
                    sphere = generate_sphere(np.zeros(shape=np.shape(bold_ref_data)),
                                             x_resamp,
                                             y_resamp,
                                             z_resamp,
                                             radius,
                                             1,
                                             2,
                                             mask_constrain,
                                             wm_seg_resampled_thresh_bin,
                                             'bold')

                elif bold_ref_dims == (3.0, 3.0, 3.0):
                    sphere = generate_sphere(np.zeros(shape=np.shape(bold_ref_data)),
                                             x_resamp,
                                             y_resamp,
                                             z_resamp,
                                             radius,
                                             1,
                                             3,
                                             mask_constrain,
                                             wm_seg_resampled_thresh_bin,
                                             'bold')

                else:
                    print('Unusual BOLD resolution: ' + str(bold_ref_dims))
                    continue

                # Save a reference sphere for visualization
                if coord == 0:
                    sphere_img = nib.Nifti1Image(sphere, affine=bold_ref.affine)
                    sphere_fname = ospj(out_paths['masks'], sub + '_wm electrode_mask_reference_radius-' + str(radius) + 'mm.nii.gz')
                    nib.save(sphere_img, sphere_fname)

                # Compile ROIs for visualization
                spherical_rois_reference_wm[radius] = spherical_rois_reference_wm[radius] + sphere

                # Mask BOLD data by sphere
                bold_data_masked = np.zeros(np.shape(bold_data))
                for vols in range(bold_nvols):
                    vol = bold_data[:,:,:,vols]
                    bold_data_masked[:,:,:,vols] = vol*sphere
                
                # Average BOLD data
                num_voxels = np.sum(sphere)
                nonzero_voxels = np.where(sphere > 0)
                bold_timeseries_sphere = np.empty(shape=(int(num_voxels), bold_nvols))
                for vox in range(np.shape(nonzero_voxels)[1]):
                    x = nonzero_voxels[0][vox]
                    y = nonzero_voxels[1][vox]
                    z = nonzero_voxels[2][vox]
                    bold_timeseries_sphere[vox, :] = bold_data_masked[x,y,z,:]
                bold_timeseries_wm[radius][:,coord] = np.mean(bold_timeseries_sphere, axis=0)

            # Remove electrode masks
            os.remove(ospj(out_paths['masks'],  sub + '_mask_wm_electrode-' + str(coord) + '.nii.gz'))
            os.remove(ospj(out_paths['masks'],  sub + '_mask_wm_electrode-' + str(coord) + '_resamp.nii.gz'))

        for radius in radii:

            # Save compiled ROIs for visualization
            spherical_rois_reference_img = nib.Nifti1Image(spherical_rois_reference_wm[radius], affine=bold_ref.affine)
            nib.save(spherical_rois_reference_img, ospj(out_paths['masks'], sub + '_wm_electrode_masks_radius-' + str(radius) + 'mm.nii.gz'))
        
            # Save averaged BOLD timeseries
            np.savetxt(ospj(out_paths['raw'], sub + '_wm_bold_timeseries_radius-' + str(radius) + 'mm.csv'), bold_timeseries_wm[radius], delimiter=',')

            # Create correlation matrix from averaged BOLD timeseries
            bold_timeseries_df = pd.DataFrame(bold_timeseries_wm[radius], columns=labels_wm)
            bold_timeseries_df = bold_timeseries_df.filter(ieeg_connectivity_wm_labels)
            
            # Create correlation matrix
            bold_timeseries_conn_pearson = bold_timeseries_df.corr(method='pearson')
            bold_timeseries_conn_spearman = bold_timeseries_df.corr(method='spearman')

            # Only keep upper triangle
            bold_timeseries_conn_pearson_upper = bold_timeseries_conn_pearson.mask(np.triu(np.ones(bold_timeseries_conn_pearson.shape, dtype=np.bool_)))
            bold_timeseries_conn_spearman_upper = bold_timeseries_conn_spearman.mask(np.triu(np.ones(bold_timeseries_conn_spearman.shape, dtype=np.bool_)))

            # Save connectivity matrix (spherical ROIs)
            bold_timeseries_conn_pearson_upper.to_csv(out_paths['connectivity'] + '/' + sub + '_wm_fmri_connectivity-pearson_radius-' + str(radius) + 'mm.csv')
            bold_timeseries_conn_spearman_upper.to_csv(out_paths['connectivity'] + '/' + sub + '_wm_fmri_connectivity-spearman_radius-' + str(radius) + 'mm.csv')

        print('--Grey matter...')

        labels_gm = []
        for coord in range(len(coords_gm)):

            # Extract coordinates
            coord_line = coords_gm.iloc[coord,:]
            label,x,y,z = (coord_line['label'], coord_line['x'], coord_line['y'], coord_line['z'])
            labels_gm.append(label)

            # Generate small sherical mask around coordinate
            coord_df = pd.DataFrame({'x': x, 'y': y, 'z': z, 't': 0}, index=[0])
            coord_mask = generate_sphere(np.zeros(shape=np.shape(t1_data)),
                                         int(round(float(coord_df['x']),0)),
                                         int(round(float(coord_df['y']),0)),
                                         int(round(float(coord_df['z']),0)),
                                         1,
                                         1,
                                         1,
                                         0,
                                         gm_seg_resampled_noWM_bin,
                                         't1')
            coord_mask_img = nib.Nifti1Image(coord_mask, affine=t1.affine)
            coord_mask_fname = ospj(out_paths['masks'], sub + '_mask_gm_electrode-' + str(coord) + '.nii.gz')
            nib.save(coord_mask_img, coord_mask_fname)

            # Resample mask into BOLD space, extracting its COG
            coord_mask_resampled_fname = ospj(out_paths['masks'], sub + '_mask_gm_electrode-' + str(coord) + '_resamp.nii.gz')
            subprocess.call(['antsApplyTransforms',
                        '-i', coord_mask_fname,
                        '-o', coord_mask_resampled_fname,
                        '-r', bold_ref_file,
                        '-n', 'Gaussian',
                        '-t', xfm_file_t1_bold])
            coord_mask_resampled = nib.load(coord_mask_resampled_fname)
            coord_mask_resampled_data = coord_mask_resampled.get_fdata()
            coord_resampled = scipy.ndimage.center_of_mass(coord_mask_resampled_data)
            x_resamp,y_resamp,z_resamp = [int(round(coord_resampled[0],0)), int(round(coord_resampled[1],0)), int(round(coord_resampled[2],0))]

            # Iterate over radii
            for radius in radii:

                # Generate spherical ROIs in BOLD space
                if bold_ref_dims == (2.0, 2.0, 3.0):
                    sphere = generate_sphere_anisotropic(np.zeros(shape=np.shape(bold_ref_data)),
                                                         x_resamp,
                                                         y_resamp,
                                                         z_resamp,
                                                         radius,
                                                         1,
                                                         mask_constrain,
                                                         gm_seg_resampled_noWM_bin,
                                                         'bold')

                elif bold_ref_dims == (2.0, 2.0, 2.0):
                    sphere = generate_sphere(np.zeros(shape=np.shape(bold_ref_data)),
                                             x_resamp,
                                             y_resamp,
                                             z_resamp,
                                             radius,
                                             1,
                                             2,
                                             mask_constrain,
                                             gm_seg_resampled_noWM_bin,
                                             'bold')

                elif bold_ref_dims == (3.0, 3.0, 3.0):
                    sphere = generate_sphere(np.zeros(shape=np.shape(bold_ref_data)),
                                             x_resamp,
                                             y_resamp,
                                             z_resamp,
                                             radius,
                                             1,
                                             3,
                                             mask_constrain,
                                             gm_seg_resampled_noWM_bin,
                                             'bold')

                else:
                    print('Unusual BOLD resolution: ' + str(bold_ref_dims))
                    continue

                # Save a reference sphere for visualization
                if coord == 0:
                    sphere_img = nib.Nifti1Image(sphere, affine=bold_ref.affine)
                    sphere_fname = ospj(out_paths['masks'], sub + '_gm_electrode_mask_reference_radius-' + str(radius) + 'mm.nii.gz')
                    nib.save(sphere_img, sphere_fname)

                # Compile ROIs for visualization
                spherical_rois_reference_gm[radius] = spherical_rois_reference_gm[radius] + sphere

                # Mask BOLD data by sphere
                bold_data_masked = np.zeros(np.shape(bold_data))
                for vols in range(bold_nvols):
                    vol = bold_data[:,:,:,vols]
                    bold_data_masked[:,:,:,vols] = vol*sphere

                # Average BOLD data
                num_voxels = np.sum(sphere)
                nonzero_voxels = np.where(sphere > 0)
                bold_timeseries_sphere = np.empty(shape=(int(num_voxels), bold_nvols))
                for vox in range(np.shape(nonzero_voxels)[1]):
                    x = nonzero_voxels[0][vox]
                    y = nonzero_voxels[1][vox]
                    z = nonzero_voxels[2][vox]
                    bold_timeseries_sphere[vox, :] = bold_data_masked[x,y,z,:]
                bold_timeseries_gm[radius][:,coord] = np.mean(bold_timeseries_sphere, axis=0)

        for radius in radii:

            # Save compiled ROIs for visualization
            spherical_rois_reference_img = nib.Nifti1Image(spherical_rois_reference_gm[radius], affine=bold_ref.affine)
            nib.save(spherical_rois_reference_img, ospj(out_paths['masks'], sub + '_gm_electrode_masks_radius-' + str(radius) + 'mm.nii.gz'))

            # Save averaged BOLD timeseries
            np.savetxt(ospj(out_paths['raw'], sub + '_gm_bold_timeseries_radius-' + str(radius) + 'mm.csv'), bold_timeseries_gm[radius], delimiter=',')

            # Create correlation matrix from averaged BOLD timeseries
            bold_timeseries_df = pd.DataFrame(bold_timeseries_gm[radius], columns=labels_gm)
            bold_timeseries_df = bold_timeseries_df.filter(ieeg_connectivity_gm_labels)

            # Create correlation matrix
            bold_timeseries_conn_pearson = bold_timeseries_df.corr(method='pearson')
            bold_timeseries_conn_spearman = bold_timeseries_df.corr(method='spearman')

            # Only keep upper triangle
            bold_timeseries_conn_pearson_upper = bold_timeseries_conn_pearson.mask(np.triu(np.ones(bold_timeseries_conn_pearson.shape, dtype=np.bool_)))
            bold_timeseries_conn_spearman_upper = bold_timeseries_conn_spearman.mask(np.triu(np.ones(bold_timeseries_conn_spearman.shape, dtype=np.bool_)))

            # Save connectivity matrix (spherical ROIs)
            bold_timeseries_conn_pearson_upper.to_csv(out_paths['connectivity'] + '/' + sub + '_gm_fmri_connectivity-pearson_radius-' + str(radius) + 'mm.csv')
            bold_timeseries_conn_spearman_upper.to_csv(out_paths['connectivity'] + '/' + sub + '_gm_fmri_connectivity-spearman_radius-' + str(radius) + 'mm.csv')

def generate_sphere(A, x0,y0,z0, radius, value, dim, mask_constrain, seg, space):
    '''
        A: array where the sphere will be drawn
        radius : radius of circle inside A which will be filled with ones.
        x0,y0,z0: coordinates for the center of the sphere within A
        value: value to fill the sphere with
    '''
    AA = A

    for x in range(x0 - int(round(radius/dim,0)), x0 + int(round(radius/dim,0)) + 1):
        for y in range(y0 - int(round(radius/dim,0)), y0 + int(round(radius/dim,0)) + 1):
            for z in range(z0 - int(round(radius/dim,0)), z0 + int(round(radius/dim,0)) + 1):
                ''' deb: measures how far a coordinate in A is far from the center.
                        deb>=0: inside the sphere.
                        deb<0: outside the sphere.'''
                deb = int(round(radius/dim)) - ((x0-x)**2 + (y0-y)**2 + (z0-z)**2)**0.5
                if (deb) >= 0:

                    # Check that we're not looking outside of the image dimensions
                    if space == 'bold':
                        if x >= np.shape(seg)[0] or y >= np.shape(seg)[1] or z >= np.shape(seg)[2]:
                            continue

                    # Check if constraining by white matter segmentation is specified
                    if mask_constrain == 1:
                        if seg[x,y,z] == 1:
                            AA[x,y,z] = value
                    else:
                        AA[x,y,z] = value

    if np.sum(AA) < 5:
        print('WARNING - Small spherical ROI with ' + str(np.sum(AA)) + ' voxels')

    return AA

def generate_sphere_anisotropic(A, x0,y0,z0, radius, value, mask_constrain, seg, space):
    '''
        Analog of generate_spheres for anisotropic voxels, specifically 2x2x3mm
    '''
    AA = A

    for x in range(x0 - int(round(radius/2,0)), x0 + int(round(radius/2,0)) + 1):
        for y in range(y0 - int(round(radius/2,0)), y0 + int(round(radius/2,0)) + 1):
            for z in range(z0 - int(round(radius/3,0)), z0 + int(round(radius/3,0)) + 1):
                ''' deb: measures how far a coordinate in A is far from the center.
                        deb>=0: inside the sphere.
                        deb<0: outside the sphere.'''

                deb_xy = int(round(radius/2,0))  - ((x0-x)**2 + (y0-y)**2)**0.5
                deb_xz = int(round(radius/3,0)) + 1 - ((x0-x)**2 + (z0-z)**2)**0.5
                deb_yz = int(round(radius/3,0)) + 1 - ((y0-y)**2 + (z0-z)**2)**0.5

                if deb_xy >= 0 and deb_xz >= 0 and deb_yz >= 0:
                   
                    # Check that we're not looking outside of the image dimensions
                    if space == 'bold':
                        if x >= np.shape(seg)[0] or y >= np.shape(seg)[1] or z >= np.shape(seg)[2]:
                            continue

                    # Check if constraining by white matter segmentation is specified
                    if mask_constrain == 1:
                        if seg[x,y,z] == 1:
                            AA[x,y,z] = value
                    else:
                        AA[x,y,z] = value

    if np.sum(AA) < 5:
        print('WARNING - Small spherical ROI with ' + str(np.sum(AA)) + ' voxels')

    return AA

if __name__ == "__main__":
    main()
