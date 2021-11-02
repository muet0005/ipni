#!/bin/bash

idc=${1}


if [ "${idc}" == "" ] ; then
    echo "You must pass a valid procnum"
    echo "e.g., ./run_mk_fsf_luciana.sh 12345"
    exit
fi

#folder where the RSFMRI files live
DIR=/export/shared/monthly/MRIdata/nifti/155_ADOL_RT/ALL
#freesurfer subjects dir so that we can make a brain mask -> T1 brain for the feat registration
SUBJECTS_DIR=/export/shared/monthly/MRIdata/freesurfer_v4.5/155_ADOL_RT/ALL

#fsl setup
export FSLDIR=/export/shared/daily/data_analysis/rmuetzel/software/fsl/5.0.5
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=${PATH}:${FSLDIR}/bin

#python path
#we need this to use the modules in knicr/fmri/fmri_utils.py
pydir=/export/shared/daily/data_analysis/rmuetzel/scripts/genr_bitbucket/lisa
export PYTHONPATH=${PYTHONPATH}:${pydir}/py-site-packages

#set the name of the output folder for feat.
#this will be inisde of DIR/PROC/REST/
outdir_sfix='17Nov2013.feat'

#we can pre-specify a bunch of settings in this confix file.
#the python script will parse it. Take a look at the python script (knicr/fmri/fmri_utils.py) to see what can be parsed
config=/export/shared/daily/data_analysis/rmuetzel/scripts/genr_bitbucket/lisa/rsfmri/bootstrap/luciana.feat.cfg

#make a brain image quick
#we can simply use the brainmask from freesurfer
#first, move the brainmask.mgz image to native space
$FREESURFER_HOME/bin/mri_vol2vol --mov $SUBJECTS_DIR/$idc/mri/brainmask.mgz --targ $DIR/$idc/t1.nii.gz --regheader --o $DIR/$idc/REST/brainmask_2_native.nii.gz --no-save-reg
#next, make it binary
$FSLDIR/bin/fslmaths $DIR/$idc/REST/brainmask_2_native.nii.gz -bin $DIR/$idc/REST/brainmask_2_native_mask.nii.gz
#finally, multiply the brainmask_mask by the t1 to get a t1_brain
$FSLDIR/bin/fslmaths $DIR/$idc/t1.nii.gz -mul $DIR/$idc/REST/brainmask_2_native_mask.nii.gz $DIR/$idc/REST/t1_brain.nii.gz
#copy the original t1 into the REST dir....this is needed by flirt/BBR, as they need to the same basename, and live in the same folder
$FSLDIR/bin/fslmaths $DIR/$idc/t1.nii.gz $DIR/$idc/REST/t1.nii.gz
#remind yourself what you did here, since the names are all similar to what happens with using regular old BET:
touch $DIR/$idc/REST/README_T1.txt
echo "The t1 image was multiplied by the FreeSurfer brainmask image to get a brain-extracted image for FEAT preprocessin" >> $DIR/$idc/REST/README_T1.txt

t1=$DIR/$idc/REST/t1_brain.nii.gz

#make the FSF
python ${pydir}/rsfmri/bootstrap/mk_fsf_files.py -f $DIR -m ${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz -s ${idc} -o ${outdir_sfix} --t1 ${t1} --config ${config} --luciana

#run feat
${FSLDIR}/bin/feat $DIR/${idc}/REST/${idc}_${outdir_sfix}.fsf

#align the 4d output to  standard space
$FSLDIR/bin/flirt -in $DIR/$idc/REST/${idc}_${outdir_sfix}/filtered_func_data -ref $DIR/$idc/REST/${idc}_${outdir_sfix}/reg/standard -out $DIR/$idc/REST/${idc}_${outdir_sfix}/filtered_func_data_2_standard_2mm -dof 12 -applyxfm -init $DIR/$idc/REST/${idc}_${outdir_sfix}/reg/example_func2standard.mat

