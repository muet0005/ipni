#!/bin/bash

idc=${1}


export FSLDIR=/home/genr/software/fsl/5.0.4
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=${PATH}:${FSLDIR}/bin

export FREESURFER_HOME=/home/genr/software/freesurfer/5.2.0
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

if [ ! -d $TMPDIR/rmuetzel/t1/ ] ; then
	mkdir -p $TMPDIR/rmuetzel/t1
fi

rsync -rtvzL /home/genr/data/t1/${idc} $TMPDIR/rmuetzel/t1/

oDIR=$TMPDIR/rmuetzel/t1

if [ $# -ne 1 ] ; then
		echo "must supply at least 1 idc"
		exit
fi

if [ ! -d $oDIR/${idc}/average ] ; then
	echo "Cannot find $oDIR/${idc}/average....must exit..."
	exit
fi

refimg=knicr130lin2avg_t1_1mm

$FSLDIR/bin/fnirt --in=$oDIR/$idc/average/t1_idc_${idc}_1mm_restore --ref=$FSLDIR/data/custom/${refimg} --aff=$oDIR/$idc/average/t1_idc_${idc}_2_knicr130lin1avg_brain_12dof.mat --cout=$oDIR/$idc/average/t1_idc_${idc}_nonlin_2_${refimg}_cout --iout=$oDIR/$idc/average/t1_idc_${idc}_nonlin_2_${refimg} --logout=$oDIR/$idc/average/t1_idc_${idc}_nonlin_2_${refimg}.log --refmask=$FSLDIR/data/custom/knicr166_nonlin3_t1_1mm_brain_mask_rlm_edits.nii.gz

rsync -rtvzL  $TMPDIR/rmuetzel/t1/${idc}/ /home/genr/data/t1/${idc}/
