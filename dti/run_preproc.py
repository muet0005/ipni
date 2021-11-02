#!/usr/bin/env python

import knicr.dti.dti_utils as du
import os, sys
import nibabel as nb
import argparse as argparse


#First parse the command-line arguments that we need.
#set up the parser
parser = argparse.ArgumentParser(description='Run KNICR DTI processing pipeline')
#add the args you want parsed
parser.add_argument("-s", "--subject", help="subject ID number", required=True, nargs=1)
parser.add_argument("-d", "--dir", help="folder where input directories exist", required=True, nargs=1)
parser.add_argument("-t", "--template", help="Standard space image to be used for registration. Usually 1mm FMRIB58", required=True, nargs=1)
parser.add_argument("--dwi", help="the path to the dwi", required=True, nargs=1)
parser.add_argument("-g", "--nDwis", help="Number of Diffusion Weighted Images expected (default is set to 38)", required=False, nargs=1)
parser.add_argument("-v", "--voxelSize", help="desired size of voxels in mm (default is 2mm)", required=False, nargs=1)

#parse the args
args = parser.parse_args()
#required args
#what is the subject ID
subj = args.subject[0]
#where does this subject live?
dataDir = args.dir[0]
#which standard space brain to register to?
template_brain = args.template[0]
#the dwi
dwi = args.dwi[0]


if args.nDwis:
    nDwis = args.nDwis[0]
else:
    nDwis = 38

if args.voxelSize:
    voxelSize = float(args.voxelSize[0])
    trg_pix_dim = (voxelSize, voxelSize, voxelSize)
else:
    trg_pix_dim = (2.0,2.0,2.0)


#some quick error checking
if not os.path.exists(os.path.join(dataDir, subj)):
    print('cannot find subject', subj)
    print('dataDir = ', dataDir)
    sys.exit(0)
elif not os.path.exists(template_brain):
    print('cannot find template brain: ', template_brain)
elif not os.path.exists(dwi):
    print('cannot find dwi: ', dwi)

def run_dti_preproc(dtiinfo):
    dataDir, subj, dwi = dtiinfo
    dwi_nii = nb.load(dwi)
    bvecs = dwi.replace('.nii.gz', '.bvec')
    print(bvecs)
    print dataDir
    print subj
    print dwi
    try:
        nvols = dwi_nii.shape[3]
        if not nvols == nDwis:
            print 'ERROR:     Subject: ', subj, ' Only has: ', nvols, '     Expecting: ', nDwis
            return
    except:
        return
    #initialize the preprocess class
    d = du.knicrDTIprep(dataDir, subj, dwi = dwi)
    #check teh dimensions of the image
    #GE saves out super resolution, which is silly, and just makes the code run longer (i.e., bedpost especially)
    pix_dim = d.probe_pixdim()
    #if the pix dims don't match 2x2x2 (or whatever you specify in trg pix dim) this will resmaple/fix it
    if not pix_dim == trg_pix_dim:
        print 'Pixel dimensions do not match!...must reslice to target resolution:'		
        print 'current: ', pix_dim
        print 'target: ', trg_pix_dim
        d.resample_dti(trg_pix_dim)
    #run the BET brain extract tool on the first b=0 image
    d.bet()
    #do the standard FSL eddy current correction
    d.ecc()
    #rotate the gradient table based on the above ecc parameters
    d.rot_bvecs(bvec=bvecs)
    #fit the tensor
    d.fit()
    #fit the tensor again, with the RESTORE method
    #if not os.path.exists(os.path.join(dataDir, subj, 'dmri', 'camino_restore')):
    d.fit_restore(fmed=True)
    dtireg = du.knicrDTIreg(dataDir, subj)
    #register to standard space
    dtireg.stSpaceReg(template_brain)
    #run bedpost then autoPtx to get tracts
    #intialize the registrations
    d.bedpost()


run_dti_preproc([dataDir, subj, dwi])
