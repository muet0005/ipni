import os, sys
import nibabel as nb
import argparse as argparse
import shutil as shutil
import tarfile as tarfile
import multiprocessing as mp
import logging
import ipni.dti.dti_utils as du




logging.getLogger().setLevel(logging.INFO)
#logging.getLogger().setLevel(logging.DEBUG)

dataDir = '/mnt/data/genr/mrdata/GenR_MRI/bids'
derDir = '/mnt/data/genr/mrdata/GenR_MRI/bids/derivatives/dwi/fsl'

subs = os.listdir(derDir)
#subs = ['sub-1', 'sub-1793', 'sub-2224', 'sub-2631', 'sub-3079', 'sub-3465', 'sub-3836', 'sub-4281', 'sub-4653',  'sub-5041', 'sub-5407', 'sub-5792', 'sub-6185', 'sub-6595']
sess = ['ses-F05', 'ses-F09', 'ses-F13']

if os.environ.has_key('SLURM_NTASKS_PER_NODE'):
    ncores = int(os.environ['SLURM_NTASKS_PER_NODE'])
    parallel = True
    logging.info('detected multiple cores...running in parallel', ncores)
else:
    parallel = False

fData = []
for sub in subs:
    for ses in sess:
        dwi = os.path.join(derDir, sub, ses, 'dmri', 'data.nii.gz')
        fa = os.path.join(derDir, sub, ses, 'dmri', 'dtifit_wls_FA.nii.gz')
        if os.path.exists(dwi) and not os.path.exists(fa):
            fInfo = [sub, ses, dwi, dataDir, derDir]
            logging.debug(fInfo)
            fData.append(fInfo)

def run_wls(fInfo):
    sub, ses, dwi, dataDir, derDir = fInfo
    d = du.knicrDTIprep(dataDir, sub, dwi = dwi, bidsSes = ses, bidsDerDir = derDir)
    d.fit(wls = True)




if not __name__ == '__main__':
    parallel = False

#run the parallelized version of randomise
if parallel:
    logging.info('setting up parallel pool...')
    pool = mp.Pool(processes=ncores)
    logging.info('launching processing across cores...')
    pool.map(run_wls, fData)