#!/bin/bash

idc=${1}

export FSLDIR=/home/genr/software/fsl/5.0.4
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=${PATH}:${FSLDIR}/bin

for i in t1 rsfmri ; do
	if [ ! -d $TMPDIR/rmuetzel/rsfmri ] ; then
		mkdir -p $TMPDIR/rmuetzel/${i}
	fi
done

outdir_sfix='27July2013.feat'
config=/home/genr/software/bitbucket/lisa/rsfmri/bootstrap/default.feat.cfg

rsync -rtvzL /home/genr/data/rsfmri/${idc} $TMPDIR/rmuetzel/rsfmri/
rsync -rtvzL /home/genr/data/t1/${idc}/t1_idc_${idc}_0.9mm.nii.gz $TMPDIR/rmuetzel/t1/
rsync -rtvzL /home/genr/data/t1/${idc}/t1_idc_${idc}_0.9mm_brain.nii.gz $TMPDIR/rmuetzel/t1/

python /home/genr/software/bitbucket/lisa/rsfmri/bootstrap/mk_fsf_files.py -f $TMPDIR/rmuetzel/rsfmri -m ${FSLDIR}/data/custom/knicr130_t1_2mm_brain.nii.gz -s ${idc} -o ${outdir_sfix} --t1 $TMPDIR/rmuetzel/t1/t1_idc_${idc}_0.9mm_brain.nii.gz --config ${config} --genr

${FSLDIR}/bin/feat $TMPDIR/rmuetzel/rsfmri/${idc}/idc_${idc}_${outdir_sfix}.fsf

$FSLDIR/bin/flirt -in $TMPDIR/rmuetzel/rsfmri/$idc/idc_${idc}_${outdir_sfix}/filtered_func_data -ref $TMPDIR/rmuetzel/rsfmri/$idc/idc_${idc}_${outdir_sfix}/reg/standard -out $TMPDIR/rmuetzel/rsfmri/$idc/idc_${idc}_${outdir_sfix}/filtered_func_data_2_standard -dof 12 -applyxfm -init $TMPDIR/rmuetzel/rsfmri/$idc/idc_${idc}_${outdir_sfix}/reg/example_func2standard.mat

rsync -rtvz $TMPDIR/rmuetzel/rsfmri/${idc}/ /home/genr/data/rsfmri/${idc}/
