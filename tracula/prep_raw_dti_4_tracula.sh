#!/bin/bash

shdir=/home/genr/software/bitbucket/lisa/tracula

DIR=$1
subj=$2

export FREESURFER_HOME=/home/genr/software/freesurfer/5.2.0
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

# First, fix the GE upsampling....
# The dicom header says the acquisition matrix is 110 x 110, but the nifti output are 256x 256
# The dicom header also says the FOV is 220mm.  This is double 110, which means the in-plane res should be 2x2mm
# The in-plane resolution of the nii data is .8594 x .8594.  to fit that into a 220mm fov, you need a 256 matrix
#the dicom header also says the slice thickness is 2mm.

###this needs to be made smarter...it should probe the vox dims to be sure....
if [ ! -e ${DIR}/${subj}/dti_idc_${subj}_2mm.nii.gz ] ; then
	#cp ${DIR}/${subj}/dti_idc_${subj}.nii.gz ${DIR}/${subj}/dti_idc_${subj}_raw.nii.gz
	$FREESURFER_HOME/bin/mri_convert ${DIR}/${subj}/dti_idc_${subj}.nii.gz -vs 2 2 2 ${DIR}/${subj}/dti_idc_${subj}_2mm.nii.gz
fi
#tracula wants the bvecs and bvals to be transposed in a 3 column format
#this is strange because bedpost still requires a 3-row format, which tracula must recreate anyway!
if [ ! -e ${DIR}/${subj}/dti_idc_${subj}.tracula.bvec ] ; then
	python $shdir/transpose_b_tbl.py ${DIR}/${subj}/dti_idc_${subj}.bvec ${DIR}/${subj}/dti_idc_${subj}.tracula.bvec
fi
if [ ! -e ${DIR}/${subj}/dti_idc_${subj}.tracula.bval ] ; then
	python $shdir/transpose_b_tbl.py ${DIR}/${subj}/dti_idc_${subj}.bval ${DIR}/${subj}/dti_idc_${subj}.tracula.bval
fi

