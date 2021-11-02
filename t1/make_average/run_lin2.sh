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

#run the 2nd linear reg.  This time use 12dof
if [ ! -e $oDIR/${idc}/average/t1_idc_${idc}_2_knicr130lin1avg_brain_12dof.nii.gz ] ; then
	echo "12dof reg"
	$FSLDIR/bin/flirt -in $oDIR/${idc}/average/t1_idc_${idc}_1mm_brain_restore -ref $FSLDIR/data/custom/knicr130lin1avg_brain -out $oDIR/${idc}/average/t1_idc_${idc}_2_knicr130lin1avg_brain_12dof -omat $oDIR/${idc}/average/t1_idc_${idc}_2_knicr130lin1avg_brain_12dof.mat -dof 12
fi
#apply the 6dof xfm to the whole head 
if [ ! -e $oDIR/${idc}/average/t1_idc_${idc}_2_knicr130lin1avg_12dof.nii.gz ] ; then
	$FSLDIR/bin/flirt -in $oDIR/${idc}/average/t1_idc_${idc}_1mm_restore -ref $FSLDIR/data/custom/knicr130lin1avg_brain -out $oDIR/${idc}/average/t1_idc_${idc}_2_knicr130lin1avg_12dof -applyxfm -init $oDIR/${idc}/average/t1_idc_${idc}_2_knicr130lin1avg_brain_12dof.mat -dof 12
fi

rsync -rtvzL  $TMPDIR/rmuetzel/t1/${idc}/ /home/genr/data/t1/${idc}/
