import os as os
import sys as sys
import csv as csv
import nibabel as nb
import numpy as np
from multiprocessing import Pool
import knicr.fmri.fmri_utils as kfmriu


#input dir
#must contain the fllowing:
#melodic_IC.nii.gz, subject_list.txt (full paths to filtered func data), 
INPUT = '/Volumes/rbraid/mr_data_idc/aug2013_final/rsfmri/dual_reg_test'
#output dir
OUTPUT = os.path.join(INPUT, 'output_no_moco')

feat_sfix = '_27July2013.feat'
filtered_func_data = 'filtered_func_data_2_standard_3mm.nii.gz'

#initialize the class
dual_reg = kfmriu.dual_reg(INPUT, OUTPUT, featdir_sfix=feat_sfix, ff_data_name=filtered_func_data)
#read the input list....this must be called subjects_list.txt in INPUT, or else specified as kwargs in this step...
inputs = dual_reg.read_input_list()

if not os.path.exists(OUTPUT):
	os.makedirs(os.path.join(OUTPUT, 'stage1'))
	os.makedirs(os.path.join(OUTPUT, 'stage2'))


print "Found: ", len(inputs), " data sets in subject list."
print "Making mask..."
dual_reg.generate_mask(parallel=False, ncores=4)

print 'dual regression stage 1...'
dual_reg.dr_stage1()

print 'dual regression stage 2...'
dual_reg.dr_stage2()


#some code from the interwebs on how to fix the pickle/multiproce problem:
#http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
#
#
#def _pickle_method(method):
#	func_name = method.im_func.__name__
#	obj = method.im_self
#	cls = method.im_class
#	return _unpickle_method, (func_name, obj, cls)
#
#def _unpickle_method(func_name, obj, cls):
#	for cls in cls.mro():
#		try:
#			func = cls.__dict__[func_name]
#		except KeyError:
#			pass
#		else:
#			break
#	return func.__get__(obj, cls)
#
#import copy_reg
#import types
#copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


