import os as os
import sys as sys
import knicr.fmri.fmri_utils as kfmriu
import argparse as argparse

#First parse the command-line arguments that we need.
#set up the parser
parser = argparse.ArgumentParser(description='Generate FSF files for FEAT analysis.')
#add the args you want parsed
parser.add_argument("-f", "--feat_dir", help="Location of all subjects' feat folders", required=True, nargs=1)
parser.add_argument("-m", "--standard_brain", help="Path to brain in standard space.", required=True, nargs=1)
parser.add_argument("-s", "--subject", help="subject id number within the feat_dir",  required=True, nargs=1)
parser.add_argument("-o", "--outdir", help="The name of the output folder that will go in feat_dir/subject/", required=True, nargs=1)
parser.add_argument("--t1", help="T1-weighted structural scan", required=False, nargs=1)
parser.add_argument("--niiPfix", help="prefix of the nifti image to be analyzed", required=False, nargs=1)
parser.add_argument("--niiSfix", help="suffix of the nifti image to be analyzed", required=False, nargs=1)
parser.add_argument("--genr", help="assume the outdir, nifti pfix and sfix are ordered in a certain way", required=False, action="store_true")
parser.add_argument("--luciana", help="assume the outdir, nifti pfix and sfix are ordered in a certain way", required=False, action="store_true")
parser.add_argument("--scan_id", help="if both a subject and scan session ID number is assigned, this is the scan-specific id", required=False, nargs=1)
parser.add_argument("--config", help="specify a configuration file with various feat parameters.", required=False, nargs=1)
#parse the args
args = parser.parse_args()

#grab the required ones, and give them reasonable names
feat_dir = args.feat_dir[0]
subject = args.subject[0]
standard_brain = args.standard_brain[0]
outdir_sfix = str('_') + args.outdir[0]

t1 = False

#check if this is the genr setup...if so, give it the default nomenclature
if args.genr:
	print "using GenR format for nifti and output files"
	fmrinii_pfix = 'rs-160_idc_'
	outdir = os.path.join(feat_dir, subject, 'idc_' + str(subject) + outdir_sfix)
	fmrinii = os.path.join(feat_dir, subject, fmrinii_pfix + str(subject) + '.nii.gz')
	fsfFile = os.path.join(feat_dir, subject, 'idc_' + str(subject) + outdir_sfix + '.fsf')
	if args.t1:
		t1 = args.t1[0]
		
#check if luciana setup is specified...give the default nomenclature
elif args.luciana:
	print 'Luciana lab flag set...using those presets...'
	fmrinii_pfix = 'rest_fmri'
	outdir = os.path.join(feat_dir, subject, 'REST', str(subject) + outdir_sfix)
	fmrinii = os.path.join(feat_dir, subject, 'REST', fmrinii_pfix + '.nii.gz')
	fsfFile = os.path.join(feat_dir, subject, 'REST', str(subject) + outdir_sfix + '.fsf')
	if args.t1:
		t1 = args.t1[0]
#otherwise give some kind of default (probably will never work, unelss the user pre-specifies this setup)
else:
	print "using default format for nifti and output files"
	fmrinii_pfix = 'resting_state_'
	outdir = os.path.join(feat_dir, subject, 'idc_' + str(subject) + outdir_sfix)
	fmrinii = os.path.join(feat_dir, subject, fmrinii_pfix + str(subject) + '.nii.gz')
	fsfFile = os.path.join(feat_dir, subject, 'idc_' + str(subject) + outdir_sfix + '.fsf')

config = False
if args.config:
	config = args.config[0]
		
print "output directory = ", outdir
print "fmri nifti input = ", fmrinii
print "fsf file = ", fsfFile

#first intialize the class
genFSF = kfmriu.genFSLfsf(fsfFile, fmrinii, outdir, standard_brain, highresniis=t1, config=config)
#run the function to make compute the number of TRs
nvols = kfmriu.calc_nvols(fmrinii)
#make the fsf file
genFSF.genFeatfsf(nvols)
