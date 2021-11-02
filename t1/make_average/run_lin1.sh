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

#idcs=(1 3 48 79 95 100 132 162 250 251 255 288 326 392 408 436 468 487 496 518 524 536 549 553 576 609 625 633 677 747 752 770 793 841 889 919 1019 1041 1072 1113 1265 1285 1326 1394 1399 1473 1505 1761 1776 1780 1844 2029 2030 2060 2127 2204 2465 2606 2774 2914 2956 2977 3049 3123 3150 3169 3188 3202 3226 3380 3414 3472 3521 3540 3573 3610 3628 3734 3781 3801 3829 3834 3846 3886 3916 3970 4049 4133 4137 4155 4237 4242 4304 4350 4447 4549 4577 4620 4652 4665 4944 4948 5052 5150 5164 5206 5221 5401 5447 5475 5530 5927 5932 6000 6039 6135 6203 6255 6321 6400 6404 6482 6734 6755 6818 6978 7005 7051 7062 7528 7839 7930 8274 8438)

if [ $# -ne 1 ] ; then
		echo "must supply at least 1 idc"
		exit
fi

if [ ! -d $oDIR/${idc}/average ] ; then
		mkdir $oDIR/${idc}/average
fi
	

echo $idc "linear 6dof"
if [ ! -e $oDIR/${idc}/average/t1_idc_${idc}_2_MNI152_T1_1mm_6dof.nii.gz -a ! -e $oDIR/${idc}/average/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
	#touch an isrunning so you can run on multiple cores..
	#resample from 512x512 -> 256x256, 0.9mm iso
	if [ ! -e $oDIR/${idc}/average/t1_idc_${idc}_1mm.nii.gz ] ; then
		echo "Resampling"
		$FSLDIR/bin/fslmaths ${oDIR}/${idc}/t1_idc_${idc}_0.9mm.nii.gz $oDIR/${idc}/average/t1_idc_${idc}_1mm.nii.gz
		$FSLDIR/bin/fslmaths ${oDIR}/${idc}/t1_idc_${idc}_0.9mm_brain.nii.gz $oDIR/${idc}/average/t1_idc_${idc}_1mm_brain.nii.gz
	fi
	#intensity normalize the brain image with FAST
	if [ ! -e $oDIR/${idc}/average/t1_idc_${idc}_1mm_brain_restore.nii.gz  ] ; then
		echo "intensity normalization"
		$FSLDIR/bin/fast -b -B -o $oDIR/${idc}/average/t1_idc_${idc}_1mm_brain $oDIR/${idc}/average/t1_idc_${idc}_1mm_brain
	fi
	#now run a 6dof reg with the instensity/bias corrected brain image
	if [ ! -e $oDIR/${idc}/average/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof.nii.gz ] ; then
		echo "Inital 6dof reg"
		$FSLDIR/bin/flirt -in $oDIR/${idc}/average/t1_idc_${idc}_1mm_brain_restore -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain -out $oDIR/${idc}/average/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof -omat $oDIR/${idc}/average/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof.mat -dof 6
	fi
	#correct the head image for intensity bias using the fast bias field output....NOTE...this will not correct the skull/non-brain!!!
	$FSLDIR/bin/fslmaths $oDIR/${idc}/average/t1_idc_${idc}_1mm -div $oDIR/${idc}/average/t1_idc_${idc}_1mm_brain_bias $oDIR/${idc}/average/t1_idc_${idc}_1mm_restore
	#apply the 6dof xfm to the whole head 
	$FSLDIR/bin/flirt -in $oDIR/${idc}/average/t1_idc_${idc}_1mm_restore -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain -out $oDIR/${idc}/average/t1_idc_${idc}_2_MNI152_T1_1mm_6dof -applyxfm -init $oDIR/${idc}/average/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof.mat -dof 6
fi

rsync -rtvzL  $TMPDIR/rmuetzel/t1/${idc}/ /home/genr/data/t1/${idc}/

