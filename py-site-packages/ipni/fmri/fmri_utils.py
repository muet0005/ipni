import os as os
import sys as sys
import csv as csv
from nibabel import nifti1 as nii
import nipype.interfaces.fsl as fsl
import random
from multiprocessing import Pool
import numpy as np

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


class genFSLfsf():
    """
    Class for generating  FSL MELODIC FSF files
    
    The initialization takes 3 required arguments (i.e., the ones that don't have defaults):
    
    fsfFile = the name of the resulting fsfFile
    
    fmriniis = a single string fmri nifti file path, or a list of multiple paths to nifti files.    
    
    standard_brain = the path and file name of the standard space brain you're using.  
        
    **Option arguments:
    
    config='/path_to_your_configuration_file/config.txt'
    
    A file can be read in that will overwrite a number of defaults such as te, tr, standard_brain, etc.
    This file would have values that look like this:
    
    tr=2.0
    te=30
    regstandard=/home/fsl/data/standard/mni_152_t1_1mm_brain.nii.gz
    
    Most options have the same variable names as the ones you find in FSL's fsf files for feat.
    
    highresniis = a path to a single structural high res image, or a list with multiple.  THIS MUST MATCH SUBJECT ID ORDER OF THE FMRINII LIST!
    
    """
    def __init__(self, fsfFile, fmriniis, outdir, standard_brain, **kwargs):
        self.fsfFile = fsfFile; self.fmriniis = fmriniis; self.outdir = outdir; self.standard_brain = standard_brain
        #set a few default options.  These are overwirtten if specified via a config file
        self.cluster = True
        self.tr = 2.0; self.te = 30
        ###feat-specific options###
        if not self.cluster:
            self.featwatcher = True
        else:
            self.featwatcher = False
        self.sscleanup = False
        self.ndelete = 4
        self.smooth = 8
        self.reghighres = False
        self.regstandard = '\"' + standard_brain + '\"'
        self.dim = False
        self.ndimensions = 16
        self.reghighres_search = 90
        self.regstandard_search = 90
        self.regstandard_nonlinear_yn = False
        self.analysis = 7
        self.reghighres_dof = 'BBR'
        self.sliceTiming = 0
        #read in passed arguments...overwrite these defaults if needed.
        for i in kwargs.keys():
            if i == 'config':
                print 'Config file passed....reading options...'
                config_file = open(kwargs[i], 'r')
                for args in config_file:
                    arg,val = args.strip('\n').split('=')
                    if arg == 'tr':
                        print 'TR=', val
                        self.tr = val
                    elif arg == 'te':
                        print 'TE=', val
                        self.te = val
                    elif arg == 'cluster':
                        print 'Cluster=', val
                        self.cluster = bool(val)
                    elif arg == 'sscleanup':
                        print 'ssl cleanup=', val
                        self.sscleanup = bool(val)
                    elif arg == 'ndelete':
                        print 'delete vols=', val
                        self.ndelete = val
                    elif arg == 'smooth':
                        print 'smoothing kernel=', val
                        self.smooth = val
                    elif arg == 'reghighres':
                        print 'use high res=', val
                        self.reghighres = bool(val)
                    elif arg == 'regstandard':
                        print 'standard brain=', val
                        self.regstandard = '\"' + val + '\"'
                    elif arg == 'dim':
                        print 'dimensionality=', val
                        self.dim = bool(val)
                    elif arg == 'ndimensions':
                        print 'number of dimensions=', val
                        self.ndimensions = val
                    elif arg == 'reghighres_search':
                        print 'reg high res search=', val
                        self.reghighres_search = val
                    elif arg == 'regstandard_search':
                        print 'reg standard search=', val
                        self.regstandard_search = val
                    elif arg == 'regstandard_nonlinear_yn':
                        print 'reg standard nonlinear=', val
                        self.regstandard_nonlinear_yn = val
                    elif arg == 'analysis':
                        print 'analysis type=', val
                        self.analysis = val
                    elif arg == 'reghighres_dof':
                        print 'reghighres_dof=', val
                        self.reghighres_dof = val
                    elif arg == 'st':
                        print 'Slice timing= ', val
                        self.sliceTiming = val
                config_file.close()
            elif i == 'highresniis':
                print 'high res image=', kwargs[i]
                self.highresniis = kwargs[i]
                if not self.highresniis == False:
                    self.reghighres = True
    
    def genFeatfsf(self, nvols):
        if not nvols:
            print "cannot find nifti volume to calculate number of TRs...must exit..."
            sys.exit(0)
        fsf = []
        fsf.append('set fmri(version) 6.00' + '\n')
        fsf.append('set fmri(inmelodic) 0' + '\n')
        fsf.append('set fmri(level) 1' + '\n')
        fsf.append('set fmri(analysis) ' + str(self.analysis) + '\n')
        fsf.append('set fmri(relative_yn) 0' + '\n')
        fsf.append('set fmri(help_yn) 1' + '\n')
        fsf.append('set fmri(featwatcher_yn) ' + str(int(self.featwatcher)) + '\n')
        fsf.append('set fmri(sscleanup_yn) ' + str(int(self.sscleanup)) + '\n')
        fsf.append('set fmri(outputdir) \"' + self.outdir + '\"' + '\n')
        fsf.append('set fmri(tr) ' + str(self.tr) + '\n')
        fsf.append('set fmri(npts) ' + str(nvols) + '\n')
        fsf.append('set fmri(ndelete) ' + str(self.ndelete) + '\n')
        fsf.append('set fmri(tagfirst) 1' + '\n')
        if type(self.fmriniis) is list:
            fsf.append('set fmri(multiple) ' + str(len(self.fmriniis)) + '\n')
        else:
            fsf.append('set fmri(multiple) 1' + '\n')
        fsf.append('set fmri(inputtype) 1' + '\n')
        prestats = {1: True, 2: False, 3: True, 4: False, 6: False, 7: True}
        fsf.append('set fmri(filtering_yn) 1' + '\n')
        fsf.append('set fmri(brain_thresh) 10' + '\n')
        fsf.append('set fmri(critical_z) 5.3' + '\n')
        fsf.append('set fmri(noise) 0.66' + '\n')
        fsf.append('set fmri(noisear) 0.34' + '\n')
        fsf.append('set fmri(newdir_yn) 0' + '\n')
        fsf.append('set fmri(mc) 1' + '\n')
        fsf.append('set fmri(sh_yn) 0' + '\n')
        #TODO: need to work code in for this in the future; possibly consider using xnat dcm data base for params
        fsf.append('set fmri(regunwarp_yn) 0' + '\n')
        fsf.append('set fmri(dwell) 0.7' + '\n')
        fsf.append('set fmri(te) ' + str(self.te) + '\n')
        fsf.append('set fmri(signallossthresh) 10' + '\n')
        fsf.append('set fmri(unwarp_dir) y-' + '\n')
        fsf.append('set fmri(st) ' + str(self.sliceTiming) + '\n')
        fsf.append('set fmri(st_file) \"\"' + '\n')
        fsf.append('set fmri(bet_yn) 1' + '\n')
        fsf.append('set fmri(smooth) ' + str(self.smooth) + '\n')
        fsf.append('set fmri(norm_yn) 0' + '\n')
        fsf.append('set fmri(perfsub_yn) 0' + '\n')
        fsf.append('set fmri(temphp_yn) 1' + '\n')
        fsf.append('set fmri(templp_yn) 0' + '\n')
        fsf.append('set fmri(melodic_yn) 0' + '\n')
        fsf.append('set fmri(stats_yn) 1' + '\n')
        fsf.append('set fmri(prewhiten_yn) 1' + '\n')
        fsf.append('set fmri(motionevs) 0' + '\n')
        fsf.append('set fmri(robust_yn) 0' + '\n')
        fsf.append('set fmri(mixed_yn) 2' + '\n')
        fsf.append('set fmri(evs_orig) 1' + '\n')
        fsf.append('set fmri(evs_real) 2' + '\n')
        fsf.append('set fmri(evs_vox) 0' + '\n')
        fsf.append('set fmri(ncon_orig) 1' + '\n')
        fsf.append('set fmri(ncon_real) 1' + '\n')
        fsf.append('set fmri(nftests_orig) 0' + '\n')
        fsf.append('set fmri(nftests_real) 0' + '\n')
        fsf.append('set fmri(constcol) 0' + '\n')
        poststats = {1: False, 2: False, 3: False, 4: True, 6: True, 7: True}
        fsf.append('set fmri(poststats_yn) ' + str(int(poststats[self.analysis])) + '\n')
        fsf.append('set fmri(threshmask) \"\"' + '\n')
        fsf.append('set fmri(thresh) 3' + '\n')
        fsf.append('set fmri(prob_thresh) 0.05' + '\n')
        fsf.append('set fmri(z_thresh) 2.3' + '\n')
        fsf.append('set fmri(zdisplay) 0' + '\n')
        fsf.append('set fmri(zmin) 2' + '\n')
        fsf.append('set fmri(zmax) 8' + '\n')
        fsf.append('set fmri(rendertype) 1' + '\n')
        fsf.append('set fmri(bgimage) 1' + '\n')
        fsf.append('set fmri(tsplot_yn) 1' + '\n')
        fsf.append('set fmri(reg_yn) 1' + '\n')
        fsf.append('set fmri(reginitial_highres_yn) 0' + '\n')
        fsf.append('set fmri(reginitial_highres_search) 90' + '\n')
        fsf.append('set fmri(reginitial_highres_dof) 3' + '\n')
        fsf.append('set fmri(reghighres_yn) ' + str(int(self.reghighres)) + '\n')
        fsf.append('set fmri(reghighres_search) ' + str(self.reghighres_search) + '\n')
        fsf.append('set fmri(reghighres_dof) ' + str(self.reghighres_dof) + '\n')
        fsf.append('set fmri(regstandard_yn) 1' + '\n')
        fsf.append('set fmri(regstandard) ' + self.regstandard + '\n')
        fsf.append('set fmri(regstandard_search) ' + str(self.regstandard_search) + '\n')
        fsf.append('set fmri(regstandard_dof) 12' + '\n')
        fsf.append('set fmri(regstandard_nonlinear_yn) ' + str(int(self.regstandard_nonlinear_yn)) + '\n')
        fsf.append('set fmri(regstandard_nonlinear_warpres) 10' + '\n')
        fsf.append('set fmri(paradigm_hp) 100' + '\n')
        fsf.append('set fmri(ncopeinputs) 0' + '\n')            
        if not type(self.fmriniis) is list:
            fsf.append('set feat_files(1) \"' + self.fmriniis + '\"\n')
        else:
            for i in range(1, len(self.fmriniis)+1):
                 fsf.append('set feat_files(' + str(i) + ') \"' + self.fmriniis[i-1] + '\"\n')
        fsf.append('set fmri(confoundevs) 0' + '\n')
        if self.reghighres:
            if not type(self.highresniis) is list:
                fsf.append('set highres_files(1) \"' + self.highresniis + '\"\n' )
            else:
                for i in range(1, len(self.highresniis)+1):
                    fsf.append('set highres_files(' + str(i) + ') \"' + self.highresniis[i-1] + '\"\n')
        fsf.append('set fmri(evtitle1) \"\"' + '\n')
        fsf.append('set fmri(shape1) 0' + '\n')
        fsf.append('set fmri(convolve1) 2' + '\n')
        fsf.append('set fmri(convolve_phase1) 0' + '\n')
        fsf.append('set fmri(tempfilt_yn1) 1' + '\n')
        fsf.append('set fmri(deriv_yn1) 1' + '\n')
        fsf.append('set fmri(skip1) 0' + '\n')
        fsf.append('set fmri(off1) 30' + '\n')
        fsf.append('set fmri(on1) 30' + '\n')
        fsf.append('set fmri(phase1) 0' + '\n')
        fsf.append('set fmri(stop1) -1' + '\n')
        fsf.append('set fmri(gammasigma1) 3' + '\n')
        fsf.append('set fmri(gammadelay1) 6' + '\n')
        fsf.append('set fmri(ortho1.0) 0' + '\n')
        fsf.append('set fmri(ortho1.1) 0' + '\n')
        fsf.append('set fmri(con_mode_old) orig' + '\n')
        fsf.append('set fmri(con_mode) orig' + '\n')
        fsf.append('set fmri(conpic_real.1) 1' + '\n')
        fsf.append('set fmri(conname_real.1) \"\"' + '\n')
        fsf.append('set fmri(con_real1.1) 1' + '\n')
        fsf.append('set fmri(con_real1.2) 0' + '\n')
        fsf.append('set fmri(conpic_orig.1) 1' + '\n')
        fsf.append('set fmri(conname_orig.1) \"\"' + '\n')
        fsf.append('set fmri(con_orig1.1) 1' + '\n')
        fsf.append('set fmri(conmask_zerothresh_yn) 0' + '\n')
        fsf.append('set fmri(conmask1_1) 0' + '\n')
        fsf.append('set fmri(thresh_yn) 1' + '\n')
        fsf.append('set fmri(mmthresh) 0.5' + '\n')
        fsf.append('set fmri(ostats) 0' + '\n')
        fsf.append('set fmri(ts_model_mat) \"\"' + '\n')
        fsf.append('set fmri(ts_model_con) \"\"' + '\n')
        fsf.append('set fmri(subject_model_mat) \"\"' + '\n')
        fsf.append('set fmri(subject_model_con) \"\"' + '\n')
        fsf.append('set fmri(alternative_example_func) \"\"' + '\n')
        fsf.append('set fmri(alternative_mask) \"\"' + '\n')
        fsf.append('set fmri(init_initial_highres) \"\"' + '\n')
        fsf.append('set fmri(init_highres) \"\"' + '\n')
        fsf.append('set fmri(init_standard) \"\"' + '\n')
        fsf.append('set fmri(overwrite_yn) 0' + '\n')
        f = open(self.fsfFile, 'w')
        f.writelines(fsf)
        f.close()
    
    def genMelodicfsf(self, nvols):
        if not nvols:
            print "cannot find nifti volume to calculate number of TRs...must exit..."
            sys.exit(0)
        fsf = []
        fsf.append('set fmri(version) 6.00' + '\n')
        fsf.append('set fmri(inmelodic) 1' + '\n')
        fsf.append('set fmri(level) 1' + '\n')
        fsf.append('set fmri(analysis) ' + str(self.analysis) + '\n')
        fsf.append('set fmri(relative_yn) 0' + '\n')
        fsf.append('set fmri(help_yn) 1' + '\n')
        fsf.append('set fmri(featwatcher_yn) ' + str(int(self.featwatcher)) + '\n')
        fsf.append('set fmri(sscleanup_yn) ' + str(int(self.sscleanup)) + '\n')
        fsf.append('set fmri(outputdir) \"' + self.outdir + '\"' + '\n')
        fsf.append('set fmri(tr) ' + str(self.tr) + '\n')
        fsf.append('set fmri(npts) ' + str(nvols) + '\n')
        fsf.append('set fmri(ndelete) ' + str(self.ndelete) + '\n')
        fsf.append('set fmri(tagfirst) 1' + '\n')
        if type(self.fmriniis) is list:
            fsf.append('set fmri(multiple) ' + str(len(self.fmriniis)) + '\n')
        else:
            fsf.append('set fmri(multiple) 1' + '\n')
        fsf.append('set fmri(inputtype) 1' + '\n')
        prestats = {1: True, 2: False, 3: True, 4: False, 6: False, 7: True}
        fsf.append('set fmri(filtering_yn) 1' + '\n')
        fsf.append('set fmri(brain_thresh) 10' + '\n')
        fsf.append('set fmri(critical_z) 5.3' + '\n')
        fsf.append('set fmri(noise) 0.66' + '\n')
        fsf.append('set fmri(noisear) 0.34' + '\n')
        fsf.append('set fmri(newdir_yn) 0' + '\n')
        fsf.append('set fmri(mc) 1' + '\n')
        fsf.append('set fmri(sh_yn) 0' + '\n')
        #TODO: need to work code in for this in the future; possibly consider using xnat dcm data base for params
        fsf.append('set fmri(regunwarp_yn) 0' + '\n')
        fsf.append('set fmri(dwell) 0.7' + '\n')
        fsf.append('set fmri(te) ' + str(self.te) + '\n')
        fsf.append('set fmri(signallossthresh) 10' + '\n')
        fsf.append('set fmri(unwarp_dir) y-' + '\n')
        fsf.append('set fmri(st) ' + str(self.sliceTiming) + '\n')
        fsf.append('set fmri(st_file) \"\"' + '\n')
        fsf.append('set fmri(bet_yn) 1' + '\n')
        fsf.append('set fmri(smooth) ' + str(self.smooth) + '\n')
        fsf.append('set fmri(norm_yn) 0' + '\n')
        fsf.append('set fmri(perfsub_yn) 0' + '\n')
        fsf.append('set fmri(temphp_yn) 1' + '\n')
        fsf.append('set fmri(templp_yn) 0' + '\n')
        fsf.append('set fmri(melodic_yn) 0' + '\n')
        fsf.append('set fmri(stats_yn) 1' + '\n')
        fsf.append('set fmri(prewhiten_yn) 1' + '\n')
        fsf.append('set fmri(motionevs) 0' + '\n')
        fsf.append('set fmri(robust_yn) 0' + '\n')
        fsf.append('set fmri(mixed_yn) 2' + '\n')
        fsf.append('set fmri(evs_orig) 1' + '\n')
        fsf.append('set fmri(evs_real) 1' + '\n')
        fsf.append('set fmri(evs_vox) 0' + '\n')
        fsf.append('set fmri(ncon_orig) 1' + '\n')
        fsf.append('set fmri(ncon_real) 1' + '\n')
        fsf.append('set fmri(nftests_orig) 0' + '\n')
        fsf.append('set fmri(nftests_real) 0' + '\n')
        fsf.append('set fmri(constcol) 0' + '\n')
        poststats = {1: False, 2: False, 3: False, 4: True, 6: True, 7: True}
        fsf.append('set fmri(poststats_yn) ' + str(int(poststats[self.analysis])) + '\n')
        fsf.append('set fmri(threshmask) \"\"' + '\n')
        fsf.append('set fmri(thresh) 3' + '\n')
        fsf.append('set fmri(prob_thresh) 0.05' + '\n')
        fsf.append('set fmri(z_thresh) 2.3' + '\n')
        fsf.append('set fmri(zdisplay) 0' + '\n')
        fsf.append('set fmri(zmin) 2' + '\n')
        fsf.append('set fmri(zmax) 8' + '\n')
        fsf.append('set fmri(rendertype) 1' + '\n')
        fsf.append('set fmri(bgimage) 1' + '\n')
        fsf.append('set fmri(tsplot_yn) 1' + '\n')
        fsf.append('set fmri(reg_yn) 1' + '\n')
        fsf.append('set fmri(reginitial_highres_yn) 0' + '\n')
        fsf.append('set fmri(reginitial_highres_search) 90' + '\n')
        fsf.append('set fmri(reginitial_highres_dof) 3' + '\n')
        fsf.append('set fmri(reghighres_yn) ' + str(int(self.reghighres)) + '\n')
        fsf.append('set fmri(reghighres_search) ' + str(self.reghighres_search) + '\n')
        fsf.append('set fmri(reghighres_dof) ' + str(self.reghighres_dof) + '\n')
        fsf.append('set fmri(regstandard_yn) 1' + '\n')
        fsf.append('set fmri(regstandard) ' + self.regstandard + '\n')
        fsf.append('set fmri(regstandard_search) ' + str(self.regstandard_search) + '\n')
        fsf.append('set fmri(regstandard_dof) 12' + '\n')
        fsf.append('set fmri(regstandard_nonlinear_yn) ' + str(int(self.regstandard_nonlinear_yn)) + '\n')
        fsf.append('set fmri(regstandard_nonlinear_warpres) 10' + '\n')
        fsf.append('set fmri(paradigm_hp) 100' + '\n')
        fsf.append('set fmri(ncopeinputs) 0' + '\n')            
        if not type(self.fmriniis) is list:
            fsf.append('set feat_files(1) \"' + self.fmriniis + '\"\n')
        else:
            for i in range(1, len(self.fmriniis)+1):
                 fsf.append('set feat_files(' + str(i) + ') \"' + self.fmriniis[i-1] + '\"\n')
        fsf.append('set fmri(confoundevs) 0' + '\n')
        if self.reghighres:
            if not type(self.highresniis) is list:
                fsf.append('set highres_files(1) \"' + self.highresniis + '\"\n' )
            else:
                for i in range(1, len(self.highresniis)+1):
                    fsf.append('set highres_files(' + str(i) + ') \"' + self.highresniis[i-1] + '\"\n')
        fsf.append('set fmri(regstandard_res) 2.0' + '\n')
        fsf.append('set fmri(varnorm) 1' + '\n')
        fsf.append('set fmri(dim_yn) 1' + '\n')
        fsf.append('set fmri(dim) 1' + '\n')
        fsf.append('set fmri(icaopt) 1' + '\n')
        fsf.append('set fmri(thresh_yn) 1' + '\n')
        fsf.append('set fmri(mmthresh) 0.5' + '\n')
        fsf.append('set fmri(ostats) 0' + '\n')
        fsf.append('set fmri(ts_model_mat) \"\"' + '\n')
        fsf.append('set fmri(ts_model_con) \"\"' + '\n')
        fsf.append('set fmri(subject_model_mat) \"\"' + '\n')
        fsf.append('set fmri(subject_model_con) \"\"' + '\n')
        fsf.append('set fmri(alternative_mask) \"\"' + '\n')
        fsf.append('set fmri(init_initial_highres) \"\"' + '\n')
        fsf.append('set fmri(init_highres) \"\"' + '\n')
        fsf.append('set fmri(init_standard) \"\"' + '\n')
        fsf.append('set fmri(overwrite_yn) 0' + '\n')
        f = open(self.fsfFile, 'w')
        f.writelines(fsf)
        f.close()
    


class gen_bootstrap_samples():
    """
    Class to generate random list of subjects for melodic analyses
    
    Functions include:
    
    read_subj_list to read a list of Rnumbers/IDCs
    
    gen_samples to generate the list of random id numbers
    
    gen_melodic_lisa to generate a text file to be passed to FSLs melodic command line.
    """
    def __init__(self):
        self.subjs = []
        self.samples = {}
    
    def read_subj_list(self, subj_list):
        """
        Read in a subject list.  This assumes "subj_list" is a text file with each subject id number on a new line.
        
        The data read in here is passed to gen_samples and gen_melodic_list.
        """
        if len(self.subjs) > 0:
            print 'WARNING:     It appears read_subj_list has already been run.  To prevent problems, subject list has been reset...this is only a warning.'
            self.subjs = []
        sl = open(subj_list, 'r')
        for line in sl:
            if not str(line).replace('\n', '') == '':
                self.subjs.append(str(line).replace('\n', ''))
    
    def gen_samples(self, startSample, nSamples, N):
        """
        Script to generate a list of of random subject ids.
        
        Within this list, you have nSamples lists.  Each nSamples list holds N subjects.
        
        nSamples = the number of lists of random subjects you want
        N = the number of subjects per list.
        
        """
        nSubjs = len(self.subjs)
        if len(self.samples) > 0:
            print 'WARNING:     It appears gen_samples has already been run.  To prevent problems, samples list has been reset....this is only a warning.'
            self.samples = {}
        for n in range(startSample, startSample+nSamples):
            sample_subjs = []
            while len(sample_subjs) < N:
                rand_int = random.randint(0, nSubjs-1)
                try:
                    sample_subjs.index(self.subjs[rand_int])
                except:
                    sample_subjs.append(self.subjs[rand_int])
            self.samples[n] = sample_subjs
    
    def gen_melodic_lisa(self, startSample, fmri_dir, feat_pfix, feat_sfix, melodic_sample_dir, samples_sfix):
        """
        Will generate text files to feed the melodic command line utility with a list of filtered func data from feat.
        
        Args:
        
        fmri_dir = where the feat processed fmri datasets are located.  This assumes the following directory hierarchy:
            feat_dir/subject/subject.feat
        
        feat_sfix = The naming convention for the feat folder within the subject folder.  This assumes the folder starts with the subject id.
            feat_sfix = 1_Apr2013.feat, where "1" is the subject id number.
        
        melodic_sample_dir = where the subject lists for melodic will be saved out.
        
        samples_sfix = the suffix to be used for the above mentioned subject lists.
            
        """
        nSamples = len(self.samples)
        for s in range(startSample, startSample+nSamples):
            odir = os.path.join(melodic_sample_dir, 'sample.' + str(s) + '.' + samples_sfix)
            if not os.path.exists(odir):
                os.makedirs(odir)
            oFile = os.path.join(melodic_sample_dir, 'sample.' + str(s) + '.' + samples_sfix, 'sample.' + str(s))
            if not os.path.exists(oFile):
                f = open(oFile, 'w')
                for l in self.samples[s]:
                    f.write(os.path.join(fmri_dir, str(l),  feat_pfix + str(l) + feat_sfix, 'filtered_func_data_2_standard.nii.gz' + '\n'))
                f.close()
            else:
                print 'WARNING:     Output file: ', oFile, '  Already exists...will not overwrite! Please remove existing files, change your output suffix, or change the output folder.'
    
    def gen_FixedMelodic_lisa(self, startSample, fmri_dir, fourDsFix, melodic_sample_dir, samples_sfix):
        """
        Will generate text files to feed the melodic command line utility with a list of filtered func data from feat.
        
        This is an update to allow input of 4D files in a single folder (so we don't have to send the whole feat/melodic folder struc over)
        
        
        """
        nSamples = len(self.samples)
        for s in range(startSample, startSample+nSamples):
            odir = os.path.join(melodic_sample_dir, 'sample.' + str(s) + '.' + samples_sfix)
            if not os.path.exists(odir):
                os.makedirs(odir)
            oFile = os.path.join(melodic_sample_dir, 'sample.' + str(s) + '.' + samples_sfix, 'sample.' + str(s))
            if not os.path.exists(oFile):
                f = open(oFile, 'w')
                for l in self.samples[s]:
                    f.write(os.path.join(fmri_dir, 'idc_' + str(l) +  '_' + fourDsFix + '.nii.gz \n'))
                f.close()
            else:
                print 'WARNING:     Output file: ', oFile, '  Already exists...will not overwrite! Please remove existing files, change your output suffix, or change the output folder.'
    


class dual_reg():
    """
    Class with the fsl dual regression tools
    The initalize function takes 2 required arguments, an input folder and an output folder name.
    The assumption make is that the input folder is stored within the folder which holds all of your feat preprocessed data.
    E.g.,
    
    Feat dir will all of your subjects:
    
    /export/raid/my_data
    
    In dir would be:
    
    /export/raid/my_data/dual_reg
    
    The output dir is just a folder name that will go into indir.
    
    E.g.,
    
    output
    
    The directory structure is then made automatically:
    /export/raid/my_data/dual_reg --> where you can pre-load inputs, such as melodic_IC.nii.gz and subject_list.txt
    /export/raid/my_data/dual_reg/output -->the parent output folder
    /export/raid/my_data/dual_reg/output/stage1 --> stage 1 and masking output
    /export/raid/my_data/dual_reg/output/stage2 --> stage 2 output
    
    There are 4 main steps to this.
    1.) read_input_list will ready the file "subject_list.txt" in your indir.  If you want to pre-specify this, you can with the "subj_list" 
    flag when calling the function
    2.) generate_mask...this will generate a common mask based on all input images.
    3.) dr_stage1 will run the first fsl_glm, to obtain the temporal regressors based on the input ic maps.
    4.) dr_stage2 will run the second fsl_glm to obtain the subject specific spatial maps, based on the input temporal regressions.
        By default with dr_stage2, regress_moco=True.  Thus, the motion parameters from feat preprocessing "prefiltered_func_data_mcf.par"
        will be added to the model, and thus 6 components of (non)interest will be added to the output.  Setting regress_moco to false when
        calling the function will disable this to the default dual_regression command line behavior, but this is not recommended in most cases.
    """
    def __init__(self, indir, outdir, **kwargs):
        self.indir = indir
        self.outdir = outdir
        self.ics = os.path.join(self.indir, 'melodic_IC.nii.gz')
        if not os.path.exists(self.ics):
            print "Cannot locate IC inputs...These should be in input dir and called, melodic_IC.niigz"
            print "must exit..."
            sys.exit(0)
        self.subj_list = os.path.join(self.indir, 'subject_list.txt')
        if not os.path.exists(self.subj_list):
            print "WARNING....cannot find default subject list (subject_list.txt) in input dir...will check for input from read_input_list..."
            self.subj_list = False
        self.featdir_sfix = '.feat'
        self.ff_data_name = 'filtered_func_data_2_standard.nii.gz'
        #override some defaults set above
        for i in kwargs.keys():
            if i == 'featdir_sfix':
                self.featdir_sfix = kwargs[i]
            elif i == 'ff_data_name':
                self.ff_data_name = kwargs[i]
    
    def read_input_list(self, **kwargs):
        for i in kwargs.keys():
            if i == 'subj_list':
                self.subj_list = os.path.join(self.indir, kwargs[i])
        if not self.subj_list or not os.path.exists(self.subj_list):
            print "Cannot locate input subject list.  Please ensure one of the following: "
            print "1.) The list is in input dir and called subject_list.txt"
            print "The list should have 1 subject ID per line."
            print "OR"
            print "2.) You add a full path to the read_input_list function as: subj_list=PATH_TO_SUBJECT_LIST"
            print "must exit..."
            sys.exit(0)
        self.subjects = []
        #open the file
        f = open(self.subj_list, 'r')
        #use the csv reader in case there are goofy characters present
        csv_reader = csv.reader(f, delimiter=' ')
        #read in each case
        for subj in csv_reader:
            print subj
            self.subjects.append(subj[0])
        f.close()
        return self.subjects
    
    def generate_mask(self, **kwargs):
        self.mask_list = []
        parallel = False
        ncores=2
        for i in kwargs.keys():
            if i == 'parallel':
                if kwargs[i]:
                    parallel = True
            elif i == 'ncores':
                try:
                    ncores = int(kwargs[i])
                except:
                    print 'number of cores not specified properly..', ncores
        def dual_reg_mask(subj):
            iFile = os.path.join(os.path.dirname(self.indir), subj, 'idc_' + subj + self.featdir_sfix, self.ff_data_name)
            oFile = os.path.join(self.outdir, 'stage1', 'mask_idc_' + subj + '.nii.gz')
            self.mask_list.append(oFile)
            fslmaths = fsl.ImageMaths(in_file=iFile, op_string= '-Tstd -bin', out_file=oFile, output_type='NIFTI_GZ', out_data_type='char')
            fslmaths.run()
        
        if not parallel:
            for subj in self.subjects:
                dual_reg_mask(subj)
        #else:
        #   pool = Pool(processes=ncores)
        #   pool.map(dual_reg_mask, self.subjects)
        #once the masks are all made, then make the average:
        allFile = os.path.join(self.outdir, 'stage1', 'maskALL.nii.gz')
        fslmerge = fsl.Merge(dimension='t', terminal_output='stream',in_files=self.mask_list, merged_file=allFile, output_type='NIFTI_GZ')
        fslmerge.run()
        oFile = os.path.join(self.outdir, 'stage1', 'mask.nii.gz')
        fslmaths = fsl.ImageMaths(in_file=allFile, op_string='-Tmin', out_file=oFile)
        fslmaths.run()
    
    def dr_stage1(self, **kwargs):
        desnorm = True
        for i in kwargs.keys():
            if i == 'desnorm':
                if kwargs[i]:
                    desnorm = True
                else:
                    desnorm = False
        if desnorm:
            opts_str = '--demean --des_norm'
        else:
            opst_str = '--demean'
        for subj in self.subjects:
            iFile = os.path.join(os.path.dirname(self.indir), subj, 'idc_' + subj + self.featdir_sfix, self.ff_data_name)
            oFile = os.path.join(self.outdir, 'stage1', 'dr_stage1_idc_' + subj + '.txt')
            mask = os.path.join(self.outdir, 'stage1', 'mask.nii.gz')
            fsl_glm = fsl.GLM(in_file=iFile, design=self.ics, terminal_output='stream', out_file=oFile, mask=mask, options=opts_str, output_type='NIFTI_GZ')
            fsl_glm.run()
    
    def dr_stage2(self, **kwargs):
        regress_moco = True
        desnorm = True
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
        for subj in self.subjects:
            dr_s1_txt = os.path.join(self.outdir, 'stage1', 'dr_stage1_idc_' + subj + '.txt')
            moco_txt = os.path.join(os.path.dirname(self.indir), subj, 'idc_' + subj + self.featdir_sfix, 'mc', 'prefiltered_func_data_mcf.par')
            dr_s1 = read_txt_file(dr_s1_txt)
            if regress_moco:
                moco_pars = read_txt_file(moco_txt)
                for i in range(0, len(dr_s1)):
                    for j in range(0,6):
                        dr_s1[i].append(moco_pars[i][j])
                dr_s1_moco_txt = os.path.join(self.outdir, 'stage1', 'dr_stage1_moco_idc_' + subj + '.txt')
                write_txt_file(dr_s1, dr_s1_moco_txt)
                designFile = dr_s1_moco_txt
            else:
                designFile = dr_s1_txt
            iFile = os.path.join(os.path.dirname(self.indir), subj, 'idc_' + subj + self.featdir_sfix, self.ff_data_name)
            oFile = os.path.join(self.outdir, 'stage2', 'dr_stage2_idc_' + subj + '.nii.gz')
            ozFile = os.path.join(self.outdir, 'stage2', 'dr_stage2_Z_idc_' + subj + '.nii.gz')
            mask = os.path.join(self.outdir, 'stage1', 'mask.nii.gz')
            if desnorm:
                opts_str = '--demean --des_norm'
            else:
                opst_str = '--demean'
            fsl_glm = fsl.GLM(in_file=iFile, design=designFile, terminal_output='stream', out_file=oFile, mask=mask, options=opts_str, output_type='NIFTI_GZ')
            fsl_glm.run()
            obname = os.path.join(self.outdir, 'stage2', 'dr_stage2_idc_' + subj + '_ic')
            fslsplit = fsl.Split(dimension='t', in_file=oFile, out_base_name=obname, terminal_output='stream', output_type='NIFTI_GZ')
            fslsplit.run()
    



