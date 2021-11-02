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
parser.add_argument("-a", "--atpxLib", help="Where the autoPTX folder resides", required=True, nargs=1)

#parse the args
args = parser.parse_args()
#required args
#what is the subject ID
subj = args.subject[0]
#where does this subject live?
dataDir = args.dir[0]
#where are the autoptx files?
autoPtxLib = args.atpxLib[0]

#some quick error checking
if not os.path.exists(os.path.join(dataDir, subj)):
    print('cannot find subject', subj)
    print('dataDir = ', dataDir)
    sys.exit(0)
elif not os.path.exists(autoPtxLib):
    print('cannot find autoptx library: ', autoPtxLib)


def run_autoptx(dtiinfo):
	dataDir, subj = dtiinfo
	print dataDir
	print subj
	kdti_aptx = du.knicrAutoPtx(dataDir, subj, autoPtxLib)
	kdti_aptx.autoPtx_1_preproc()
	kdti_aptx.autoPtx_2_launchTractography()


run_autoptx([dataDir, subj])
