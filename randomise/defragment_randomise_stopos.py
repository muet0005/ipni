import os as os
import sys as sys
import csv as csv
import nibabel as nb
import numpy as np
from nibabel import nifti1 as nii
import nipype.interfaces.fsl as fsl
import socket as socket
import datetime as datetime
import argparse as argparse

DIR = '/home/genr/data/randomise_tests'
poolname = 'list.pool'
poolFile = open(os.path.join(DIR, poolname + '.data'))
stat_outs = ['vox_corrp_tstat', 'tfce_corrp_tstat']
contrasts = 4


__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "0.2"

def read_pool_data(pFile):
    data = []
    poolCSV = csv.reader(pFile, delimiter=' ')
    for line in poolCSV:
        data.append(line)
    return data


def parse_pool_data(pData):
    poolDict = {}
    for line in pData:
        inFile = line[0]
        statbase = line[1] ; seed = int(line[2]) ; nperm = int(line[3])
        if not poolDict.has_key(inFile):
            entry = [statbase, seed, nperm]
            poolDict[inFile] = entry
        else:
            if poolDict[inFile][0] == statbase and poolDict[inFile][1] < seed:
                poolDict[inFile][1] = seed
    return poolDict


poolData = read_pool_data(poolFile)
poolFile.close()
poolinfo = parse_pool_data(poolData)

for inFile in poolinfo.keys():
    design = poolinfo[inFile][0]
    ncores = poolinfo[inFile][1]
    nperms_per_core = poolinfo[inFile][2]
    print 'Input file: ', inFile
    print 'Design name: ', design
    print 'Number of seeds: ', ncores
    print 'Number of permutations per seed: ', nperms_per_core
    obn = inFile + '_' + design
    parallelDir = os.path.join(DIR, obn, 'parallel')
    defragDir = os.path.join(DIR, 'defragmented', obn, 'defragmented')
    tDir = os.path.join(DIR, 'defragmented', obn, 'tstat')
    if not os.path.exists(tDir):
        os.makedirs(tDir)
    if not os.path.exists(defragDir):
        os.makedirs(defragDir)
    for seed in range(1, ncores+1):
        for contrast in range(1, int(contrasts)+1):
            for stat_out in stat_outs:
                iFile = os.path.join(parallelDir, inFile + '_' + design +'_SEED_' + str(seed) + '_' + stat_out + str(contrast) + '.nii.gz')
                oFile = iFile.replace('.nii.gz', '_defrag.nii.gz')
                if not os.path.exists(oFile):
                    fslmaths = fsl.ImageMaths(in_file=iFile, op_string='-mul ' + str(nperms_per_core), out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
                    fslmaths.run()      
    for stat_out in stat_outs:
        for contrast in range(1, int(contrasts)+1):
            addlist = []
            add_str = ''
            for seed in range(1, int(ncores)+1):
                if not seed == 1:
                    add_str = add_str + ' -add ' + os.path.join(parallelDir, inFile + '_' + design + '_SEED_' + str(seed) + '_' + stat_out + str(contrast) + '_defrag.nii.gz')
                else:
                    addlist.append(os.path.join(parallelDir, inFile + '_' + design + '_SEED_' + str(seed) + '_' + stat_out + str(contrast) + '_defrag.nii.gz'))
            oFile = os.path.join(parallelDir, inFile +  '_' + design + '_' + stat_out + str(contrast) + '_defrag_sum.nii.gz')
            if not os.path.exists(oFile):
                fslmaths = fsl.ImageMaths(in_file=addlist[0], op_string=add_str, out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
                fslmaths.run()
            iFile = oFile
            oFile = os.path.join(defragDir, inFile +  '_' + design + '_' + stat_out + str(contrast) + '_defragmented.nii.gz')
            corrected_number_of_perms = ((nperms_per_core * ncores)-(ncores)) + 1
            if not os.path.exists(oFile):
                fslmaths = fsl.ImageMaths(in_file=iFile, op_string='-div ' + str(corrected_number_of_perms), out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
                fslmaths.run()
    for contrast in range(1, int(contrasts)+1):
        iFile = os.path.join(parallelDir, inFile + '_' + design + '_SEED_1_tstat' + str(contrast) + '.nii.gz')
        oFile = os.path.join(tDir, inFile + '_' + design +'_SEED_1_tstat' + str(contrast) + '.nii.gz')
        fslmaths = fsl.ImageMaths(in_file=iFile, out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
        fslmaths.run()