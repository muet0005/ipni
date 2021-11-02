import os as os
import sys as sys
import csv as csv
import nibabel as nb
import numpy as np
from multiprocessing import Pool
from nibabel import nifti1 as nii
import nipype.interfaces.fsl as fsl
import socket as socket
import datetime as datetime

__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "0.1"

#input dir
#must contain the fllowing:
#melodic_IC.nii.gz, subject_list.txt (full paths to filtered func data), 
INPUT = '/Volumes/rbraid/mr_data_idc/aug2013_final/rsfmri/dual_regression'
#output dir
OUTPUT = os.path.join(INPUT, 'output_lmeb_dysregulation_29Nov2013')

#feat folder prefix and suffix.	 This will look something like this in the end: $feat_pfix_$subjid_$feat_sfix
feat_pfix = 'idc_'
feat_sfix = '_27July2013.feat'
##name of the data file to feed into dual_reg.	FSL recommends the filtered_func_data in standard space.
ff_data_name = 'filtered_func_data_2_standard.nii.gz'

parallel = True
ncores = 5

def read_txt_file(txtFile, **kwargs):
	data = []
	delimiter = ' '
	iscsv = True
	for i in kwargs.keys():
		if i == 'csv':
			if kwargs[i]:
				iscsv = True
			else:
				iscsv = False
		if i == 'delimiter':
			delimiter = kwargs[i]
	f = open(txtFile, 'r')
	csv_reader = csv.reader(f, delimiter=delimiter)
	for line in csv_reader:
		rm_ws = True
		while rm_ws:
			try:
				line.remove('')
			except:
				rm_ws = False
		data.append(line)
	f.close()
	return data


def write_txt_file(data, txtFile, **kwargs):
	delimiter = '\t'
	iscsv = True
	for i in kwargs.keys():
		if i == 'csv':
			if kwargs[i]:
				iscsv = True
			else:
				iscsv = False
		if i == 'delimiter':
			delimiter = kwargs[i]
	f = open(txtFile, 'w')
	csv_writer = csv.writer(f, delimiter=delimiter)
	csv_writer.writerows(data)
	f.close()


def read_input_list(**kwargs):
	subj_list = os.path.join(OUTPUT, 'subject_list.txt')
	for i in kwargs.keys():
		if i == 'subj_list':
			subj_list =	 kwargs[i]
	if not subj_list or not os.path.exists(subj_list):
		print "Cannot locate input subject list.  Please ensure one of the following: "
		print "1.) The list is in input dir and called subject_list.txt"
		print "The list should have 1 subject ID per line."
		print "OR"
		print "2.) You add a full path to the read_input_list function as: subj_list=PATH_TO_SUBJECT_LIST"
		print "must exit..."
		sys.exit(0)
	subjects = []
	#open the file
	f = open(subj_list, 'r')
	#use the csv reader in case there are goofy characters present
	csv_reader = csv.reader(f, delimiter=' ')
	#read in each case
	for subj in csv_reader:
		subjects.append(subj[0])
	f.close()
	return subjects


def write_cmd_out(cmd, oFile):
	if os.path.exists(oFile):
		f = open(oFile, 'a')
	else:
		f = open(oFile, 'w')
	f.write('#######--------NEW INSTANCE CALLED--------#######' + '\n')
	f.write('timestamp=' + datetime.datetime.now().strftime('%d%b%Y_%H:%M:%S') + '\n')
	f.write('user=' + os.getlogin() + '\n')
	f.write('pwd=' + os.getcwd() + '\n')
	f.write('domainname=' + socket.getfqdn() + '\n')
	f.write('hostname=' + socket.gethostname() + '\n')
	f.write('osUname=' + str(os.uname()) + '\n')
	f.write('PythonVersion=' + sys.version + '\n')
	f.write('FSLversion=' + fsl.base.Info.version() + '\n')
	f.write('FSLCOMMAND=' + '\n')
	f.write(cmd + '\n')
	f.close()


def generate_mask(subj, **kwargs):
	iFile = os.path.join(os.path.dirname(INPUT), subj, feat_pfix + subj + feat_sfix, ff_data_name)
	if not os.path.exists(iFile):
		print 'Cannot locate input file....must exit'
		print 'input file: ', iFile
		print 'Please ensure that subject has the input file, or adjust the subject list to remove them.'
		sys.exit(0)
	oFile = os.path.join(OUTPUT, 'stage1', 'mask_idc_' + subj + '.nii.gz')
	if os.path.exists(oFile):
	    return
	cmd_out = os.path.join(OUTPUT, 'stage1', 'stage1_masking_idc_' + subj + '.out')
	fslmaths = fsl.ImageMaths(in_file=iFile, op_string= '-Tstd -bin', out_file=oFile, output_type='NIFTI_GZ', out_data_type='char')
	write_cmd_out(fslmaths.cmdline, cmd_out)
	fslmaths.run()	


def generate_common_mask(inputs):
	mask_list = []
	for subj in inputs:
		mask = os.path.join(OUTPUT, 'stage1', 'mask_idc_' + subj + '.nii.gz')
		mask_list.append(mask)
	allFile = os.path.join(OUTPUT, 'stage1', 'maskALL.nii.gz')	
	cmd_out = os.path.join(OUTPUT, 'stage1', 'stage1_maskall_idc_' + subj + '.out')
	fslmerge = fsl.Merge(dimension='t', terminal_output='stream',in_files=mask_list, merged_file=allFile, output_type='NIFTI_GZ')
	write_cmd_out(fslmerge.cmdline, cmd_out)
	fslmerge.run()
	oFile = os.path.join(OUTPUT, 'stage1', 'mask.nii.gz')
	fslmaths = fsl.ImageMaths(in_file=allFile, op_string='-Tmin', out_file=oFile)
	fslmaths.run()


def dr_stage1(i, **kwargs):
	subj, mask = i
	ics = os.path.join(OUTPUT, 'melodic_IC.nii.gz')
	desnorm = False
	demean=True
	for i in kwargs.keys():
		if i == 'desnorm':
			if kwargs[i]:
				desnorm = True
			else:
				desnorm = False
		if i == 'ics':
			ics = kwargs[i]
	if not os.path.exists(ics):
		print "cannot locate melodic component maps for dual regression...must exit.."
		print "looked in: ", ics
		sys.exit(0)
	iFile = os.path.join(os.path.dirname(INPUT), subj, feat_pfix + subj + feat_sfix, ff_data_name)
	oFile = os.path.join(OUTPUT, 'stage1', 'dr_stage1_idc_' + subj + '.txt')
	if os.path.exists(oFile):
	    return
	mask = os.path.join(OUTPUT, 'stage1', 'mask.nii.gz')
	fsl_glm = fsl.GLM(in_file=iFile, design=ics, terminal_output='stream', out_file=oFile, des_norm=desnorm, demean=demean, mask=mask, output_type='NIFTI_GZ')
	cmd_out = os.path.join(OUTPUT, 'stage1', 'stage1_fslglm_idc_' + subj + '.out')
	write_cmd_out(fsl_glm.cmdline, cmd_out)
	fsl_glm.run()


def dr_stage2(i, **kwargs):
	subj, mask = i
	regress_moco = True
	desnorm = True
	demean = True
	for i in kwargs.keys():
		if i == 'regress_moco':
			if kwargs[i]:
				regress_moco = True
			else:
				regress_moco = False
		if i == 'desnorm':
			if kwargs[i]:
				desnorm = True
			else:
				desnorm = False
	dr_s1_txt = os.path.join(OUTPUT, 'stage1', 'dr_stage1_idc_' + subj + '.txt')
	moco_txt = os.path.join(os.path.dirname(INPUT), subj, feat_pfix + subj + feat_sfix, 'mc', 'prefiltered_func_data_mcf.par')
	dr_s1 = read_txt_file(dr_s1_txt)
	if regress_moco:
		moco_pars = read_txt_file(moco_txt)
		for i in range(0, len(dr_s1)):
			for j in range(0,6):
				dr_s1[i].append(moco_pars[i][j])
		dr_s1_moco_txt = os.path.join(OUTPUT, 'stage1', 'dr_stage1_moco_idc_' + subj + '.txt')
		write_txt_file(dr_s1, dr_s1_moco_txt)
		designFile = dr_s1_moco_txt
	else:
		designFile = dr_s1_txt
	iFile = os.path.join(os.path.dirname(INPUT), subj, feat_pfix + subj + feat_sfix, ff_data_name)
	oFile = os.path.join(OUTPUT, 'stage2', 'dr_stage2_idc_' + subj + '.nii.gz')
	ozFile = os.path.join(OUTPUT, 'stage2', 'dr_stage2_Z_idc_' + subj + '.nii.gz')
	if os.path.exists(oFile) and os.path.exists(ozFile):
	    return
	mask = os.path.join(OUTPUT, 'stage1', 'mask.nii.gz')
	fsl_glm = fsl.GLM(in_file=iFile, design=designFile, terminal_output='stream', out_file=oFile, out_z_name=ozFile, mask=mask, des_norm=desnorm, demean=demean, output_type='NIFTI_GZ')
	cmd_out = os.path.join(OUTPUT, 'stage2', 'stage2_fslglm_idc_' + subj + '.out')
	write_cmd_out(fsl_glm.cmdline, cmd_out)
	fsl_glm.run()
	obname = os.path.join(OUTPUT, 'stage2', 'dr_stage2_idc_' + subj + '_ic')
	fslsplit = fsl.Split(dimension='t', in_file=oFile, out_base_name=obname, terminal_output='stream', output_type='NIFTI_GZ')
	fslsplit.run()
	obname = os.path.join(OUTPUT, 'stage2', 'dr_stage2_idc_Z' + subj + '_ic')
	fslsplit = fsl.Split(dimension='t', in_file=ozFile, out_base_name=obname, terminal_output='stream', output_type='NIFTI_GZ')
	fslsplit.run()


if not __name__ == '__main__' and parallel:
	parallel = False
	
#read the input list....this must be called subjects_list.txt in INPUT, or else specified as kwargs in this step...
inputs = read_input_list(subj_list=os.path.join(OUTPUT, 'subject_list.txt'))

if not os.path.exists(OUTPUT):
	os.makedirs(OUTPUT)
stages=['stage1', 'stage2']
for stage in stages:
    if not os.path.exists(os.path.join(OUTPUT, stage)):
	    os.makedirs(os.path.join(OUTPUT, stage))

print "Found: ", len(inputs), " data sets in subject list."
print "Making masks....mode parallel: ", parallel

#use multiprocessing function if possible (set True/False and ncores above...)
if not parallel:
	for subj in inputs:
		generate_mask(subj)
else:
	pool = Pool(processes=ncores)
	pool.map(generate_mask, inputs)

#merge each newly generated mask and make a mask suitable for all subejcts
generate_common_mask(inputs)

mask = os.path.join(OUTPUT, 'stage1', 'mask.nii.gz')
if not os.path.exists(mask):
	print "Cannot find mask...was in made yet?"
	print "looked here: ", mask
	print "must exit..."
	sys.exit(0)

inputs_masks = []
for subj in inputs:
	s = [subj, mask]
	inputs_masks.append(s)

print 'dual regression...'
if not parallel:
	for subj in inputs_masks:
		dr_stage1(subj)
		dr_stage2(subj)
else:
	pool = Pool(processes=ncores)
	pool.map(dr_stage1, inputs_masks)
	pool.map(dr_stage2, inputs_masks)


	

	



