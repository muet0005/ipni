#!/bin/bash

idc=${1}

export FSLDIR=/home/genr/software/fsl/5.0.4
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=${PATH}:${FSLDIR}/bin

export FREESURFER_HOME=/home/genr/software/freesurfer/5.2.0
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

DIR=/home/genr/data/t1
scratch=$TMPDIR/rmuetzel/first


if [ ! -d $scratch/${idc} ] ; then
	mkdir -p $scratch/${idc}/first
fi


#rsync -rtvzL $DIR/${idc}/ ${scratch}/


$FSLDIR/bin/run_first_all -b -i ${DIR}/${idc}/t1_idc_${idc}_0.9mm_brain.nii.gz -o ${scratch}/${idc}/first/t1_idc_${idc}_0.9_brain_first


rsync -rtvzL  ${scratch}/${idc}/ ${DIR}/${idc}/
