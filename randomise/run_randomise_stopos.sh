#!/bin/bash

ic=${1}
statroot=${2}
seed=${3}
nperms=${4}

DIR=/home/genr/data/randomise_tests
oDIR=${TMPDIR}/rmuetzel/randomise/${ic}_${statroot}
gDIR=/global/rmuetzel

#make sure a seed is passed
#make sure number of permutations is passed
if [ $# != 4 ] ; then
    echo "Must supply the seed number and number of permutations for randomise..."
    exit
fi
#make they are integers
if [[ $seed != *[!0-9]* ]]; then
    echo "seed is: " ${seed}
else
    echo "Seed is not an integer...please enter a valid seed."
    exit
fi
if [[ $nperms != *[!0-9]* ]]; then
    echo "number of permutations is: " ${nperms}
else
    echo "number of permutations is not an integer...please enter a valid seed."
    exit
fi

#set up fsl
export FSLDIR=/home/genr/software/fsl/5.0.5
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=${PATH}:${FSLDIR}/bin

if [ ! -d $oDIR ] ; then
	mkdir -p $oDIR/parallel
fi
if [ ! -d $gDIR ] ; then
    mkdir -p $gDIR
fi

#echo "**************COPYING TO SCRATCH SPACE***************"
#echo "cp -r ${DIR}/${ic}.nii.gz ${gDIR}/"
#echo "cp -r ${DIR}/mask.nii.gz ${gDIR}/"
#echo "cp -r ${DIR}/${statroot}.* ${gDIR}/"
#if [! -e $gDIR/${ic}.nii.gz ] ; then
#    cp -r ${DIR}/${ic}.nii.gz ${gDIR}/
#fi
#if [ ! -e $gDIR/mask.nii.gz ] ; then
#    cp -r ${DIR}/mask.nii.gz ${gDIR}/
#fi
#if [ ! -e $gDIR/${statroot}.mat ] ; then
#    cp -r ${DIR}/${statroot}.* ${gDIR}/
#fi
#echo "**************FINISHED COPYING TO SCRATCH SPACE***************"

python $DIR/randomise_stopos.py -i ${DIR}/${ic}.nii.gz --dir ${oDIR} -d ${DIR}/${statroot}.mat -t ${DIR}/${statroot}.con -m ${DIR}/mask.nii.gz -o ${ic}_${statroot} --T --glm --nperm ${nperms} -s ${seed}

if [ ! -d $DIR/${ic}_${statroot}/parallel ] ; then
    mkdir -p $DIR/${ic}_${statroot}/parallel
fi
    

echo "**************COPYING TO HOME DIRECTORY***************"
echo "cp -r ${oDIR}/ ${DIR}/${ic}_${statroot}/"
cp -r "${oDIR}"/parallel/* ${DIR}/${ic}_${statroot}/parallel/