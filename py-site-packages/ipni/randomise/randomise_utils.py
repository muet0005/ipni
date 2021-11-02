import os as os
import sys as sys
import csv as csv
from nibabel import nifti1 as nii
import nipype.interfaces.fsl as fsl
import random
import numpy as np
import pwd
import shutil as shutil
import nibabel as nb
from nipype.interfaces.base import Undefined
import subprocess as sp
import socket as socket
import datetime as datetime
import argparse as argparse
import MySQLdb as mysqldb


__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "0.2"


def calc_nvols(fmriniis):
    if type(fmriniis) is list:
        niiFile = fmriniis[0]
    else:
        niiFile = fmriniis
    if os.path.exists(niiFile):
        n = nii.load(niiFile)
        dims = n.get_shape()
        t = dims[3]
    else:
        t = False
    return t


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


def write_txt(data, txtFile):
    t = open(txtFile, 'w')
    if isinstance(data[0], list):
        csv_writer = csv.writer(t, delimiter=',')
        csv_writer.writerows(data)
    else:
        for i in data:
            t.write(i + '\n')
    t.close()


def write_log(cmd, oFile):
    """
    write_log
    Simple function to track data provenance via the nipype package
    Arguments:
    cmd: usually the nipype function.cmdline
    oFile: the file you want to write the logging to.
    """
    if os.path.exists(oFile):
        f = open(oFile, 'a')
    else:
        f = open(oFile, 'w')
    f.write('#######--------NEW INSTANCE CALLED--------#######' + '\n')
    f.write('timestamp=' + datetime.datetime.now().strftime('%d%b%Y_%H:%M:%S') + '\n')
    try:
        u = os.getlogin()
    except OSError:
        u = pwd.getpwuid(os.getuid())[0]
    f.write('user=' + str(u) + '\n')
    f.write('pwd=' + os.getcwd() + '\n')
    f.write('domainname=' + socket.getfqdn() + '\n')
    f.write('hostname=' + socket.gethostname() + '\n')
    f.write('osUname=' + str(os.uname()) + '\n')
    f.write('PythonVersion=' + sys.version + '\n')
    f.write(cmd + '\n')
    f.write('#######--------END OF LOG--------#######' + '\n')
    f.close()



class gen_mat():
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile
    
    def read_infile(self):
        data = []
        inFile = open(self.infile, 'rU')
        csv_reader = csv.reader(inFile, delimiter=',')
        print 'Reading in data...'
        for line in csv_reader:
            data.append(line)
        inFile.close()
        nrows= len(data)-1
        ncols = len(data[0])
        npdata = np.zeros([nrows, ncols])
        print 'nsubjects, nvars', nrows, ncols
        for i in range(0 ,nrows):
            listi = i + 1
            for j in range(ncols):
                npdata[i][j] = data[listi][j]
        return npdata
    
    def two_grpT(self, **kwargs):
        """Simple 2 group t-test mat maker."""
        grp_dict = {0: [1,0], 1: [0,1]}
        data = self.read_infile()
        mat = []
        mat_line = ['/NumWaves', 2]
        mat.append(mat_line)
        mat_line = ['/NumPoints', len(data)]
        mat.append(mat_line)
        mat_line = ['/PPheights', 1.0, 1.0]
        mat.append(mat_line)
        mat_line = ['/Matrix']
        mat.append(mat_line)        
        for i in range(0, len(data)):
            #assume subject id is column 1, group is column 2, and covariate is column 3
            #mat line: [grp_var, grp0, grp1, covar]
            mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1])]
            mat.append(mat_line)
        print 'Writing out mat file...' 
        oFile = open(self.outfile, 'w')
        csv_writer = csv.writer(oFile, delimiter = ' ')
        csv_writer.writerows(mat)
        oFile.close()
        return data
    
    def two_grpT_oneCovar_MainFx(self, **kwargs):
        grp_dict = {0: [1,0], 1: [0,1]}
        data = self.read_infile()
        demean = False
        for i in kwargs.keys():
            if i == 'demean':
                if kwargs[i]:
                    demean = True
        col2mean = np.mean(np.transpose(data)[2])
        print 'column 2 mean: ', col2mean
        if demean:
            print 'Demeaning cov1 with mean'
        mat = []
        mat_line = ['/NumWaves', 3]
        mat.append(mat_line)
        mat_line = ['/NumPoints', len(data)]
        mat.append(mat_line)
        mat_line = ['/PPheights', 1.0, 1.0, 1.0]
        mat.append(mat_line)
        mat_line = ['/Matrix']
        mat.append(mat_line)        
        for i in range(0, len(data)):
            #assume subject id is column 1, group is column 2, and covariate is column 3
            #mat line: [grp_var, grp0, grp1, covar]
            if demean:
                mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), data[i][2]-col2mean]
            else:
                mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), data[i][2]]
            mat.append(mat_line)
        print 'Writing out mat file...' 
        oFile = open(self.outfile, 'w')
        csv_writer = csv.writer(oFile, delimiter = ' ')
        csv_writer.writerows(mat)
        oFile.close()
        return data
    
    def two_grpT_oneCovar_Intx(self, **kwargs):
        grp_dict = {0: [1,0], 1: [0,1]}
        data = self.read_infile()
        col2mean = np.mean(np.transpose(data)[2])
        demean = False
        for i in kwargs.keys():
            if i == 'demean':
                if kwargs[i]:
                    demean = True
        col2mean = np.mean(np.transpose(data)[2])
        print 'column 2 mean: ', col2mean
        if demean:
            print 'Demeaning cov1 with mean'
        mat = []
        mat_line = ['/NumWaves', 4]
        mat.append(mat_line)
        mat_line = ['/NumPoints', len(data)]
        mat.append(mat_line)
        mat_line = ['/PPheights', 1.0, 1.0, 1.0, 1.0]
        mat.append(mat_line)
        mat_line = ['/Matrix']
        mat.append(mat_line)        
        for i in range(0, len(data)):
            #assume subject id is column 1, group is column 2, and covariate is column 3
            #mat line: [grp_var, grp0, grp1, covar]
            if demean:
                mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), (data[i][2]-col2mean)*int(grp_dict[int(data[i][1])][0])+0, (data[i][2]-col2mean)*int(grp_dict[int(data[i][1])][1])+0]
            else:
                mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), (data[i][2])*int(grp_dict[int(data[i][1])][0])+0, (data[i][2])*int(grp_dict[int(data[i][1])][1])+0]
            mat.append(mat_line)
        print 'Writing out mat file...' 
        oFile = open(self.outfile, 'w')
        csv_writer = csv.writer(oFile, delimiter = ' ')
        csv_writer.writerows(mat)
        oFile.close()
        return data
    
    def two_grpT_twoCovar_MainFx(self, **kwargs):
        grp_dict = {0: [1,0], 1: [0,1]}
        data = self.read_infile()
        col2mean = np.mean(np.transpose(data)[2])
        col3mean = np.mean(np.transpose(data)[3])
        demean = False
        for i in kwargs.keys():
            if i == 'demean':
                if kwargs[i]:
                    demean = True
        print 'column 2 mean: ', col2mean
        print 'column 3 mean: ', col3mean
        if demean:
            print 'Demeaning cov1 with mean'
            print 'Demeaning cov2 with mean'
        mat = []
        mat_line = ['/NumWaves', 4]
        mat.append(mat_line)
        mat_line = ['/NumPoints',  len(data)]
        mat.append(mat_line)
        mat_line = ['/PPheights', 1.0, 1.0, 1.0, 1.0]
        mat.append(mat_line)
        mat_line = ['/Matrix']
        mat.append(mat_line)        
        for i in range(0, len(data)):
            #assume subject id is column 1, group is column 2, and covariate is column 3
            #mat line: [grp_var, grp0, grp1, covar]
            if demean:
                mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), data[i][2]-col2mean, data[i][3]-col3mean]
            else:    
                mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), data[i][2], data[i][3]]
            mat.append(mat_line)
        print 'Writing out mat file...' 
        oFile = open(self.outfile, 'w')
        csv_writer = csv.writer(oFile, delimiter = ' ')
        csv_writer.writerows(mat)
        oFile.close()
        return data
    
    def two_grpT_twoCovar_Intx(self, **kwargs):
        covar_inter = 1
        demean = False
        for k in kwargs.keys():
            if k == 'covar':
                covar_inter = kwargs[k]
            elif k == 'demean':
                demean = kwargs[k]
        if not (covar_inter == 1 or covar_inter == 2):
            print 'This function only allows for you to select covariate 1 or 2 for the interaction with group'
            print 'You must set covar=1 or 2...'
            print 'Cannot continue...must exit'
            sys.exit(0)
        grp_dict = {0: [1,0], 1: [0,1]}
        data = self.read_infile()
        col2mean = np.mean(np.transpose(data)[2])
        col3mean = np.mean(np.transpose(data)[3])
        print 'column 2 mean: ', col2mean
        print 'column 3 mean: ', col3mean
        if demean:
            print 'Demeaning cov1 with mean'
            print 'Demeaning cov2 with mean'
        mat = []
        mat_line = ['/NumWaves', 5]
        mat.append(mat_line)
        mat_line = ['/NumPoints', len(data)]
        mat.append(mat_line)
        mat_line = ['/PPheights', 1.0, 1.0, 1.0, 1.0, 1.0]
        mat.append(mat_line)
        mat_line = ['/Matrix']
        mat.append(mat_line)        
        for i in range(0, len(data)):
            #assume subject id is column 1, group is column 2, and covariate is column 3
            #mat line: [grp_var, grp0, grp1, covar]
            if demean:
                if covar_inter == 1:
                    mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), (data[i][2]-col2mean)*int(grp_dict[int(data[i][1])][0])+0, (data[i][2]-col2mean)*int(grp_dict[int(data[i][1])][1])+0, data[i][3]-col3mean]
                else:
                    mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), data[i][2]-col2mean, (data[i][3]-col3mean)*int(grp_dict[int(data[i][1])][0])+0, (data[i][3]-col3mean)*int(grp_dict[int(data[i][1])][1])+0]
            else:
                if covar_inter == 1:
                    mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), (data[i][2])*int(grp_dict[int(data[i][1])][0])+0, (data[i][2])*int(grp_dict[int(data[i][1])][1])+0, data[i][3]]
                else:
                    mat_line = [int(grp_dict[int(data[i][1])][0]), int(grp_dict[int(data[i][1])][1]), data[i][2], (data[i][3])*int(grp_dict[int(data[i][1])][0])+0, (data[i][3])*int(grp_dict[int(data[i][1])][1])+0]
            mat.append(mat_line)
        print 'Writing out mat file...' 
        oFile = open(self.outfile, 'w')
        csv_writer = csv.writer(oFile, delimiter = ' ')
        csv_writer.writerows(mat)
        oFile.close()
        return data
    


class gen4d():
    def __init__(self, subj_list):
        self.subj_list = subj_list
    
    def dualreg4d(self, dr_dir, dr_pfix, ics, oDir):
        for ic in ics:
            print 'writing out 4d file for component: ', ic
            fourDlist = []
            for s in range(0, self.subj_list.shape[0]):
                subj = str(int(self.subj_list[s][0]))
                pe_file = os.path.join(dr_dir, 'stage2', dr_pfix + subj + '_ic' + str(ic) + '.nii.gz')
                if not os.path.exists(pe_file):
                    print 'CANNOT FIND filtered func data for :', subj
                    print 'Looked here: ', pe_file
                    print 'Cannot continue...must exit...'
                    sys.exit(0)
                fourDlist.append(pe_file)
            oFile = os.path.join(oDir, 'dr_stage2_merged_pe_ic' + str(ic) + '.nii.gz')
            fslmerge = fsl.Merge(dimension='t', terminal_output='stream',in_files=fourDlist, merged_file=oFile, output_type='NIFTI_GZ')
            fslmerge.run()
    
    def tbss4d(self, tbss_dir, oDir, template, **kwargs):
        """
        python function to generate 4d nii files for randomise from TBSS.
        
        Required arguments:
            tbss_dir --- where the 'FA' and 'stats' folder exists.  This by default on knicr servers is here:
            /Volumes/rbraid/mr_data_idc/aug2013_final/dti/tbss
        
            oDir --- this is where you want the 4d files to go. This should be specific to the analysis you are doing
            e.g., /home/tonya/my_data/tbss/dysregulation_analysis/
        
    
        Default is to assume scalars FA, MD, RD, AD are all processed.
        Also assumes the subject files are named as such:
        
            idc_#_SCALAR_to_FMRIB58_FA_2mm_nonlinear.nii.gz
            
        
        Optional Arguments:
            
            tbss_pfix --- default is "idc_"
            This is what goes before the subject id
            String type
            
            tbss_sfix --- default is "_to_FMRIB58_FA_2mm_nonlin.nii.gz"
            This is what goes after the subject id and scalar name (eg FA, MD)
            String type
            
            Scalars --- List of Scalars 
            default is: [FA, MD, RD, AD] 
            List type
            
            opfix --- prefix for the 4d output file. suffix will be the scalar type.
            default is tbss_all_
            string type.
            
            
        """
        tbss_pfix = 'idc_'
        tbss_sfix = '_to_FMRIB58_FA_2mm_nonlin.nii.gz'
        scalars = ['FA', 'MD', 'RD', 'AD']
        opfix = 'tbss_all_'
        threshold = 0.2
        for i in kwargs.keys():
            if i == 'tbss_pfix':
                tbss_pfix = kwargs[i]
            elif i == 'tbss_sfix':
                tbss_sfix = kwargs[i]
            elif i == 'opfix':
                opfix = kwargs[i]
            elif i == 'scalars':
                scalars = kwargs[i]
                if not type(scalars) is list:
                    print 'ERROR....scalar option must be a LIST!'
                    print 'Resetting scalars to be FA only...'
                    scalars = ['FA']
            elif i == 'threshold':
                threshold = kwarg[i]
        for scalar in scalars:
            print 'writing out 4d file for component: ', scalar
            fourDlist = []
            for s in range(0, self.subj_list.shape[0]):
                subj = str(int(self.subj_list[s][0]))
                tbss_file = os.path.join(tbss_dir, 'FA', tbss_pfix + subj + '_' + scalar + tbss_sfix)
                if not os.path.exists(tbss_file):
                    print 'CANNOT FIND tbss data for :', subj, scalar
                    print 'Looked here: ', tbss_file
                    print 'Cannot continue...must exit...'
                    sys.exit(0)
                fourDlist.append(tbss_file)
            #merge the file
            logFile = os.path.join(oDir, 'tbss_nonFA.log')
            merged4d = os.path.join(oDir, opfix + scalar + '.nii.gz')
            if not os.path.exists(merged4d):
                fslmerge = fsl.Merge(dimension='t', terminal_output='stream',in_files=fourDlist, merged_file=merged4d, output_type='NIFTI_GZ')
                fslmerge.run()
                write_log(fslmerge.cmdline, logFile)
                subjinfo =  merged4d.replace('.nii.gz', '_mergeInfo.csv')
                write_txt(fourDlist, subjinfo)
            if scalar.endswith('_FA'):
                #next we make the mean FA image
                op_str = '-max 0 -Tmin -bin'
                fslmaths = fsl.ImageMaths(in_file=merged4d, op_string=op_str, out_file=os.path.join(oDir, 'mean_FA_mask.nii.gz'), output_type='NIFTI_GZ', out_data_type='char')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
                op_str = '-mas ' + os.path.join(oDir, 'mean_FA_mask.nii.gz')
                fslmaths = fsl.ImageMaths(in_file=merged4d, op_string=op_str, out_file=merged4d.replace('.nii.gz', '_masked.nii.gz'), output_type='NIFTI_GZ')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
                op_str = '-Tmean'
                fslmaths = fsl.ImageMaths(in_file=merged4d.replace('.nii.gz', '_masked.nii.gz'), op_string=op_str, out_file=os.path.join(oDir, 'mean_FA.nii.gz'), output_type='NIFTI_GZ')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
                #now make the initial skeleton
                tbss_skeleton = fsl.TractSkeleton(in_file=os.path.join(oDir, 'mean_FA.nii.gz'), skeleton_file=os.path.join(oDir, 'mean_FA_skeleton.nii.gz'))
                tbss_skeleton.run()
                write_log(tbss_skeleton.cmdline, logFile)
                print 'Creating Skeleton Mask using threshold', str(threshold)
                op_str = '-thr ' + str(threshold) + ' -bin'
                fslmaths = fsl.ImageMaths(in_file=os.path.join(oDir, 'mean_FA_skeleton.nii.gz'), op_string=op_str, out_file=os.path.join(oDir, 'mean_FA_skeleton_mask.nii.gz'), output_type='NIFTI_GZ')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
                print 'Creating skeleton distancemap (for use in projection search)'
                op_str = '-mul 1 -add 1 -add ' + os.path.join(oDir, 'mean_FA_skeleton_mask.nii.gz')
                fslmaths = fsl.ImageMaths(in_file=os.path.join(oDir, 'mean_FA_mask.nii.gz'), op_string=op_str, out_file=os.path.join(oDir, 'mean_FA_skeleton_mask_dst.nii.gz'), output_type='NIFTI_GZ')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
                distancemap = fsl.DistanceMap(in_file=os.path.join(oDir, 'mean_FA_skeleton_mask_dst.nii.gz'), distance_map=os.path.join(oDir, 'mean_FA_skeleton_mask_dst.nii.gz'))
                distancemap.run()
                write_log(distancemap.cmdline, logFile)
                #to make things easier below, make a generic-named copy of the 4d fa map. needed for skeletonization.
                shutil.copyfile(merged4d.replace('.nii.gz', '_masked.nii.gz'), os.path.join(oDir, 'all_FA.nii.gz'))
            else:
                if not os.path.exists(os.path.join(oDir, 'mean_FA_mask.nii.gz')):
                    print oDir, 'mean_FA_mask.nii.gz does not exist'
                    print 'Make sure FA is the first scalar in your list. Must exit'
                    sys.exit(0)
                #mask other scalars
                op_str = '-mas ' + os.path.join(oDir, 'mean_FA_mask.nii.gz')
                fslmaths = fsl.ImageMaths(in_file=merged4d, op_string=op_str, out_file=merged4d.replace('.nii.gz', '_masked.nii.gz'), output_type='NIFTI_GZ')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
            #skeletenize
            nb_template = nb.load(template)
            x = abs(nb_template.get_affine()[0][0])
            y = abs(nb_template.get_affine()[1][1])
            z = abs(nb_template.get_affine()[2][2])
            cing_file1mm = os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'LowerCingulum_1mm.nii.gz')
            if x == 1 and y == 1 and z == 1:
                cing_file = cing_file1mm
            elif x == 2 and y == 2 and z == 2:
                cing_file = os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'LowerCingulum_2mm.nii.gz')
                if not os.path.exists(cing_file):
                    cing_file = os.path.join(oDir, 'LowerCingulum_2mm.nii.gz')
                    flirt = fsl.FLIRT(in_file=cing_file1mm, reference=os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'MNI152_T1_2mm.nii.gz'), out_file=cing_file, args='-applyisoxfm 2')
                    flirt.run()
                    write_log(flirt.cmdline, logFile)
            else:
                cing_file=False
            if not cing_file:
                print 'Cannot find cingulum file.'
                return
            print 'Projecting all FA data onto Skeleton'
            #run the tbss skeletonize thing
            if not os.path.exists(os.path.join(oDir, 'mean_FA_mask.nii.gz')):
                print oDir, 'mean_FA_mask.nii.gz does not exist'
                print 'Make sure FA is the first scalar in your list. Must exit'
                sys.exit(0)
            final4d = merged4d.replace('.nii.gz', '_skeletonized.nii.gz')    
            tbss_skeleton = fsl.TractSkeleton(threshold=threshold, distance_map=os.path.join(oDir, 'mean_FA_skeleton_mask_dst.nii.gz'), alt_data_file=merged4d.replace('.nii.gz', '_masked.nii.gz'), 
            data_file=os.path.join(oDir, 'all_FA.nii.gz'), projected_data=final4d, in_file=os.path.join(oDir, 'mean_FA.nii.gz'), project_data=True)
            tbss_skeleton.inputs.use_cingulum_mask = Undefined
            tbss_skeleton.inputs.search_mask_file=cing_file
            tbss_skeleton.run()
            write_log(tbss_skeleton.cmdline, logFile)


