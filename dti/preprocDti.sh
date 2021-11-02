#!/bin/bash

shDir=/home/genr/queues/dti

if [ $# != 3 ] ; then
	echo "Not enough arguments supplied..."
	echo "dataDir, subject and dwi are required!"
	exit
fi

dataDir=$1
subj=$2
dwi=$3

#set up fsl
export FSLDIR=/home/genr/software/fsl/5.0.9
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=${PATH}:${FSLDIR}/bin

$shDir/run_preproc.py --subject ${subj} --dir ${dataDir} --dwi ${dwi} --template $FSLDIR/data/standard/FMRIB58_FA_1mm.nii.gz

