import os as os
import sys as sys
import csv as csv


dataDir = '/projects/0/qt06652/genr/data/xnat/nii/GenR_F9_MRI'

queueDataFile = '/home/genr/queues/dti/f9_preproc.data'

if os.path.exists(queueDataFile):
    os.remove(queueDataFile)
    

data = []
subjs = os.listdir(dataDir)
for subj in subjs:
    subjDir = os.path.join(dataDir, subj)
    dmri = os.path.join(subjDir, 'dmri')
    niis = os.listdir(subjDir)
    for nii in niis:
        if nii.startswith('dti_pa') and nii.endswith('_best.nii.gz') and not os.path.exists(dmri):
            dwi = os.path.join(subjDir, nii)
            data.append([dataDir, subj, dwi])
            break


queueData = open(queueDataFile, 'w')
queueDataWriter = csv.writer(queueData, delimiter=' ')
queueDataWriter.writerows(data)