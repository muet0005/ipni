#!/bin/bash

if [ $# != 1 ] ; then
    echo "Must specify a subject id!"
    exit
fi

#pass a subject ID
subj=${1}
#the folder where the DTI data live..
DIR=${HOME}/data/dti

####################-----------------######################
#tell this code where to find the other shells/python code it needs
shdir=/home/genr/software/bitbucket/lisa/dti

#which versions would you like to use?
#these are the same as FSLDIR and FREESURFER_HOME
fsl_version=/home/genr/software/fsl/5.0.5
fs_version=/home/genr/software/freesurfer/5.3.0

#point to the DTI scratch space (i.e., where you want it to go on the node)
dti_scratch=/scratch/rmuetzel/dti

#make the appropriate scratch folder
if [ ! -d $dti_scratch ] ; then
	mkdir -p $dti_scratch
fi

#sync the data to scratch
rsync -rtvzL ${DIR}/${idc} $dti_scratch/

#set up the FSL environment
export FSLDIR=$fsl_version
source $FSLDIR/etc/fslconf/fsl.sh
export SUBJECTS_DIR=$fs_scratch
export FREESURFER_HOME=$fs_version
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

#run bedpost
python /home/genr/software/bitbucket/lisa/dti/run_bedpost.py ${dti_scratch} ${subj}

#copy the data back to the home folder
cp -ruvL $dti_scratch/${idc}/ ${tracula_data}/
