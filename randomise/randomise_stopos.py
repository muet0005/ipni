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


__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "0.2"


def get_randomise_info(i):
	inDir,iFile,mask,mat,con,nperm = i
	iFile = os.path.join(inDir, iFile)
	r = fsl.RandomiseGENR(in_file=iFile, mask=mask, tcon=con, design_mat=mat, num_perm=nperm, tfce=True, args='-Q', terminal_output='stream')
	o = r.run()
	return o.runtime.stdout.split(' ')	


def randomise_parallel(i):
	seed,inDir,iFile,obn,mask,mat,con,perms_per_core,glm,demean = i
	oDir = os.path.join(inDir, 'parallel')
	if not os.path.exists(oDir):
	    os.makedirs(oDir)
	obn = obn + '_SEED_' + str(seed)
	obn = os.path.join(oDir, obn)
	iFile = os.path.join(inDir, iFile)
	print iFile
	if glm:
		r = fsl.RandomiseGENR(in_file=iFile, mask=mask, tcon=con, design_mat=mat, seed=seed, num_perm=int(perms_per_core), demean=demean, tfce=True, vox_p_values=True, base_name=obn, terminal_output='stream', args='--glm_output')
	else:
		r = fsl.RandomiseGENR(in_file=iFile, mask=mask, tcon=con, design_mat=mat, seed=seed, num_perm=int(perms_per_core), demean=demean, tfce=True, vox_p_values=True, base_name=obn, terminal_output='stream')
	oInfo = open(obn + '.cmdline', 'w')
	oInfo.write(r.cmdline)
	oInfo.close()
	print r.cmdline
	r.run()


#First parse the command-line arguments that we need.
#set up the parser
parser = argparse.ArgumentParser(description='Run FSLs Randomise via python')
#add the args you want parsed
parser.add_argument("-i", "--input", help="4d input file for randomise.", required=True, nargs=1)
parser.add_argument("--dir", help="folder where input files exist", required=True, nargs=1)
parser.add_argument("-m", "--mask", help="Mask for analysis (mask.nii.gz in --dir is default)", required=False, nargs=1)
parser.add_argument("-d", "--design", help="design matrix (design.mat in --dir is default)",	 required=False, nargs=1)
parser.add_argument("-t", "--con", help="t contrast file (design.con in --dir is default)",	required=False, nargs=1)
parser.add_argument("-o", "--out", help="output basename", required=True, nargs=1)
parser.add_argument("-s", "--seed", help="seed number", required=True, nargs=1)
parser.add_argument("-p", "--nperm", help="number of permutations per seed", required=True, nargs=1)
parser.add_argument("--T", help="threshold free cluster enhancement", required=False, action="store_true")
parser.add_argument("--glm", help="output glm raw files", required=False, action="store_true")
parser.add_argument("--D", help="demean design matrix and MR data (v5+ ONLY!!)", required=False, action="store_true")
#parse the args
args = parser.parse_args()
#required args
iFile = args.input[0]
inDir = args.dir[0]
obn = args.out[0]
seed = int(args.seed[0])
nperm = int(args.nperm[0])
#optional args
if args.mask:
	mask = args.mask[0]
else:
	mask = os.path.join(inDir, 'mask.nii.gz')
if args.design:
	mat = args.design[0]
else:
	mat = os.path.join(inDir, 'design.mat')
if args.con:
	con = args.con[0]
else:
	con = os.path.join(inDir, 'design.con')
if args.glm:
	glm = True
else:
	glm = False
if args.D:
	demean = True
else:
	demean = False

if not os.path.exists(mat) or not os.path.exists(con) or not os.path.exists(mask):
    print 'Cannot locate design matrix/contrast or mask'
    print 'design matrix: ', mat
    print 'design contrast: ', con
    print 'Mask: ', mask
    print 'You must supply full path to the files, or ensure design.mat/.con / mask.nii.gz exist in the input directory!'
    sys.exit(0)

#n_spec_perms = the number of permutations specified by the user


seed_info = (seed, inDir, iFile, obn, mask, mat, con, nperm, glm, demean)
randomise_parallel(seed_info)

