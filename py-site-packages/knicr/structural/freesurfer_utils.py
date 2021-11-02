import os as os
import sys as sys
import csv as csv
import pwd
import numpy as np
import shutil as shutil
import nibabel as nb
from nibabel import nifti1 as nii
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as freesurfer
from nipype.interfaces.base import Undefined
import subprocess as sp
import socket as socket
import datetime as datetime
import argparse as argparse
import MySQLdb as mysqldb
import scipy.io as sio


__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "0.2"


class qdec_clusters():
    def __init__(self, subjects_dir):
        """
        Required Arguments:
            subjects_dir = location of freesurfer data foldres.
        """
        self.subjects_dir = subjects_dir
        self.exclude = ('knicr', 'fsaverage', 'lh', 'rh')
        self.white_surf_meas = ('thickness', 'volume', 'wsurarea', 'wirmcurv', 'wirgcurv', 'wfoldind', 'wintcurv')
        self.pial_surf_meas = ('psurarea', 'pirmcurv', 'pirgcurv', 'pfoldind', 'pintcurv')
        self.lgi_surf_meas = ('lgi')
        self.surf_meas = {'white': self.white_surf_meas, 'pial': self.pial_surf_meas, 'pial_lgi': self.lgi_surf_meas}
        os.environ['SUBJECTS_DIR'] = subjects_dir
    
    def mapQdecLabel2Subjects(self, label, hemi, **kwargs):
        """
        Use mri_label2label to map label file from qdec to each subject.
        
        Required arguments:
            
            label = full path to the + name of the label file.
            
            hemi = hemisphere. lh or rh.
        
        Optional arguments:
            
            fsaverage = default is fsaverage. can specify custom template. Wherever the label was drawn on.
            
            regmethod. default is surface. volume is option as well.
             
            subjects. List type. If you want to run on a subset.  Otherwise, the default to to use the full subjects_dir listing.
            
            overwrite. boolean. if specified, label files existing in subject folder will be overwritten. otherwise skipped. default is False.
            
        """
        subjects = os.listdir(self.subjects_dir)
        measure = 'thickness'
        fsaverage = 'fsaverage'
        regmethod = 'surface'
        overwrite = False
        for i in kwargs.keys():
            if i == 'measure':
                measure = kwargs[i]
            elif i == 'fsaverage':
                fsaverage = kwargs[i]
            elif i == 'regmethod':
                regmethod = kwargs[i]
            elif i == 'overwrite':
                overwrite = True
            elif i == 'subjects':
                if type(kwargs[i]) is list:
                    subjects = kwargs[i]
        print 'Mapping labels for subjects:'
        for s in subjects:
            if s.startswith(self.exclude) or not os.path.exists(os.path.join(self.subjects_dir, s, 'label')):
                #get rid of the riff raff
                continue
            if os.path.exists(os.path.join(self.subjects_dir, s, 'label', os.path.basename(label) + '.label')):
                if not overwrite:
                    continue
            if not os.path.exists(os.path.join(self.subjects_dir, s, 'surf', hemi + '.white')):
                continue
            print s
            label2label = freesurfer.MRILabel2Label(srclabel=label + '.label', srcsubject=fsaverage, trgsubject=s, trglabel=os.path.basename(label), hemi=hemi, regmethod=regmethod)
            label2label.run()
    
    def computeQdecLabelStats(self, label, hemi, **kwargs):
        """
        Required Arguments:
            label. Just give the full path to the original qdec label. I'll figure out the rest.
            
            hemi = lh or rh
        
        Optional Arguments:
        Measure: This is the kind of stats to extract.
                    Avaialble options are:
                    self.white_surf_meas = ('thickness', 'volume', 'wsurarea', 'wirmcurv', 'wirgcurv', 'wfoldind', 'wintcurv')
                    self.pial_surf_meas = ('psurarea', 'pirmcurv', 'pirgcurv', 'pfoldind', 'pintcurv')
                    self.lgi_surf_meas = ('lgi')
                    If you choose a measure from self.white_surf_meas, it will use the *.white.stats output file to provide the stats.
                    The same goes for the pial and LGI stats.
                
        """
        subjects = os.listdir(self.subjects_dir)
        measure = 'thickness'
        overwrite = False
        for i in kwargs.keys():
            if i == 'measure':
                measure = kwargs[i]
            elif i == 'overwrite':
                overwrite = True
            elif i == 'subjects':
                if type(kwargs[i]) is list:
                    subjects = kwargs[i]
                    print 'Using custom subjects list: ', subjects
        print 'Generating stats for subjects: '
        for s in subjects:
            labelfile = os.path.join(self.subjects_dir, s, 'label', os.path.basename(label) + '.label')
            if not os.path.exists(labelfile) or s.startswith(self.exclude):
                continue
            stats = freesurfer.MRISAnatomicalStats(labelfile=labelfile, tabular=True, hemi=hemi, subject=s)
            if measure.startswith(self.white_surf_meas):
                stats.inputs.thicknessfile = hemi + '.thickness'
                surf = 'white'
            elif measure.startswith(self.pial_surf_meas):
                stats.inputs.surface = 'pial'
                surf = 'pial'
            elif measure.startswith(self.lgi_surf_meas):
                stats.inputs.thicknessfile = hemi + '.pial_lgi'
                surf = 'pial_lgi'
            oFile = os.path.join(self.subjects_dir, s, 'stats', os.path.basename(label) + '.' + surf + '.stats')
            if os.path.exists(oFile):
                if not overwrite:
                    continue
            print s
            stats.inputs.tablefile = oFile
            stats.run()
    
    def dumpQdecLabelStatsMysql(self, cursor, label, hemi, **kwargs):
        """
        Dump data from label file (stats output) into mysql db.
        
        Required:
            cursor -> Set up a mysql DB cursos. (see example below)
            label -> full path to label file, just as with other funcs in this class
            hemi -> lh or rh
        Optional:
            tblroot -> this will appar at hte front of your mysql table name (prepends the label name). default is qdec_label_
            id_col -> this is the name of the id column. default is idc.
            measure -> Available measures are:
                    self.white_surf_meas = ('thickness', 'volume', 'wsurarea', 'wirmcurv', 'wirgcurv', 'wfoldind', 'wintcurv')
                    self.pial_surf_meas = ('psurarea', 'pirmcurv', 'pirgcurv', 'pfoldind', 'pintcurv')
                    self.lgi_surf_meas = ('lgi')
                    see: http://freesurfer.net/fswiki/mris_anatomical_stats#Outputs
            subjects -> default is listing of subjects_dir. must be list type.
        
        Example:
        #set up your mysql connection
        >>>import MySQLdb as mysqldb
        
        #select a database to use:
        >>>mydb = 'freesurfer_v51_tests'
        
        #establish the con:
        >>>con = mysqldb.connect(db = mydb, read_default_file=os.path.join(os.getenv("HOME"), '.my.cnf'))
        
        #set the cursor:
        >>>cursor = con.cursor()
        
        """
        tblroot = 'qdec_label_'
        id_col = 'idc'
        measure = 'thickness'
        tbl = tblroot + os.path.basename(label).replace('.', '_').replace('-', '_') + '_' + measure
        subjects = os.listdir(self.subjects_dir)
        for i in kwargs.keys():
            if i == 'id_col':
                id_col = kwargs[i]
                print 'Resetting table ID column to: ', id_col
            elif i == 'measure':
                measure = kwargs[i]
                print 'Resetting measure to be: ', measure
            elif i == 'subjects':
                 if type(kwargs[i]) is list:
                     subjects = kwargs[i]
                     print 'Using custom subjects list: ', subjects  
            elif i == 'table':
                tbl = kwargs[i]
        try:
            cursor.execute("""create table %s (%s int)""" % (tbl, id_col))
        except:
            pass
        print 'Inserting stats into MysQL db:'
        for s in subjects:
            if s.startswith(self.exclude):
                continue
            #specify the file to read in
            if measure.startswith(self.white_surf_meas):
                stats_file = os.path.join(self.subjects_dir, s, 'stats', os.path.basename(label) + '.white.stats')
            elif measure.startswith(self.pial_surf_meas):
                stats_file = os.path.join(self.subjects_dir, s, 'stats', os.path.basename(label) + '.pial.stats')
            elif measure.startswith(self.lgi_surf_meas):
                stats_file = os.path.join(self.subjects_dir, s, 'stats', os.path.basename(label) + '.pial_lgi.stats')
            if not os.path.exists(stats_file):
                continue
            print s
            #check the table to see if the subject exists
            cursor.execute("""select %s from %s where %s=\'%s\'""" % (id_col, tbl, id_col, s))
            pn_exist = cursor.fetchone()
            if not pn_exist:
                cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, id_col, s))
            statsFile = open(stats_file, 'r')
            stats_csv = csv.reader(statsFile, delimiter=' ')
            for row in stats_csv:
                if row[0] != '#':
                    for i in range(0, len(row)):
                        try:
                            row.index('')
                            row.remove('')
                        except:
                            continue
                        roi = row[0]
                        val_std = False
                        area_meas = ('wsurarea', 'psurarea')
                        irmcurv_meas = ('wirmcurv', 'pirmcurv')
                        irgcurv_meas = ('wirgcurv', 'pirgcurv')
                        foldind_meas = ('wfoldind', 'pfoldind')
                        intcurv = ('wintcurv', 'pintcurv')
                        col = roi.replace('.', '_') + '_' + measure   
                        if measure == 'thickness':
                            val = row[4]
                            val_std = row[5]
                            col = roi.replace('.', '_') + '_thickavg'
                            col_std = roi.replace('.', '_') + '_thickavg_std'
                        elif measure == 'lgi':
                            val = row[4]
                            val_std = row[5]
                            col = roi.replace('.', '_') + '_lgi'
                            col_std = roi.replace('.', '_') + '_lgi_std'
                        elif measure == 'volume':
                            val = row[3]
                        elif measure.startswith(area_meas):
                            val = row[2]
                        elif measure.startswith(irmcurv_meas):
                            val = row[6]
                        elif measure.startswith(irgcurv_meas):
                            val = row[7]
                        elif measure.startswith(foldind_meas):
                            val = row[8]
                        elif measure.startswith(intcurv):
                            val = row[9]
                        if val == '':
                            continue
                        val = float(val)
                        try:
                            cursor.execute("""alter table %s add column (%s float)""" % (tbl, col))
                        except:
                            pass
                        cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col, val, id_col, s))
                        if val_std:
                            try:
                                cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_std))
                            except:
                                pass
                            cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_std, val_std, id_col, s))
            statsFile.close()

class surfaceAnalysis():
    def __init__(self, **kwargs):
        #defaults
        self.hemi = 'lh'
        self.fwhm = 10
        self.template = 'fsaverage'
        self.measure = 'thickness'
        self.subjectsDir = False ; self.oDir = False ; self.outBaseName = False ; self.overwrite = False ; self.fsgdFile = False
        for i in kwargs.keys():
            if i == 'hemi':
                self.hemi = kwargs[i]
            elif i == 'fwhm':
                self.fwhm = kwargs[i]
            elif i == 'template':
                self.template = kwargs[i]
            elif i == 'measure':
                self.measure = kwargs[i]
            elif i == 'subjectsDir':
                self.subjectsDir = kwargs[i]
            elif i == 'oDir':
                self.oDir = kwargs[i]
            elif i == 'outBaseName':
                self.outBaseName = kwargs[i]
            elif i == 'fsgdFile':
                self.fsgdFile = kwargs[i]
            elif i == 'overwrite':
                self.overwrite = True
            elif i == 'measures':
                self.measure = kwargs[i]
        if not self.subjectsDir or not self.oDir or not self.outBaseName or not self.fsgdFile:
            print 'You must specify the subjectsDir, oDir AND outBaseName!'
            print 'must exit'
            sys.exit(0)
        os.environ['SUBJECTS_DIR'] = self.subjectsDir 
        self.data = os.path.join(self.oDir, self.hemi + '.' + self.outBaseName + '.' + self.measure + '.' + str(self.fwhm) + '.mgh')             
    
    def fix_newline(self):
        fsgdFileOut = self.fsgdFile.replace('.txt', '.fixed.txt')
        stdout = open(fsgdFileOut, 'w')
        stdin = open(self.fsgdFile, 'r')
        opts = ['tr', '\'\\r\'', '\'\\n\'']
        p = sp.Popen(opts, stdout=stdout, stdin=stdin)
        p.wait()
        stdout.close()
        stdin.close()
    
    def fsgd_2_mat(self):
        #open the fsgd file
        fsgd = open(self.fsgdFile, 'rU')
        #read it in as a tab delimited file
        fsgd_csv = csv.reader(fsgd, delimiter='\t')
        #make a blank list to temporarily store the data
        data = []
        #make a list of the subjects in the fsgd file, in the proper order
        #vars to QA
        len_rows = False
        test_row = False
        #loop over each row
        for row in fsgd_csv:
            #wait until you get to the "Inputs", which are real data points
            if row[0].startswith('Input'):
                #QA...find the length of the first row, if you haven't already
                if not len_rows:
                    len_row = len(row)
                #make an empty list to store data for this one row
                val_insert = []
                #loop over each item in the row
                for i in range(0, len_row):
                    #skipping the first two..."input" and the subject id number
                    if i == 0 or i == 1:
                        continue
                    #QA make sure every value in the row can be made into a float! 
                    if not test_row:
                        try:
                            float(row[i])
                        except:
                            #if not, stop here....the design matrix is useless...
                            print 'error....', row[i], ' can not be converted to int or float'
                            print row
                            sys.exit(0)
                    #otherwise, append that data point into the subject-specific insert list
                    val_insert.append(row[i])
                #once the subject-specific insert-list is full, append it to the full dataset with the other subjects
                data.append(val_insert)
                #this assumes a row has now been tested for float-compatibility...turn it off for future lines
                test_row = True
        #now, to make a matlab MAT format file, we need to start with numpy friendly structure
        nrows = len(data)
        ncols = len_row-2
        print 'subjects: ', nrows, ' Variables: ', ncols
        #make an empty matrix
        data_matrix = np.zeros([nrows, ncols])
        #now populate that matrix with the data from above
        for r in range(0, nrows):
            for c in range(0, ncols):
                data_matrix[r, c] = data[r][c]
        #make an out file name for the matlab MAT file...just rename .txt from the input fsgd to .mat
        self.fsgdMat = self.fsgdFile.replace('.txt', '.mat')
        print 'saving MAT file: ', self.fsgdMat
        sio.savemat(self.fsgdMat, {'X':data_matrix})
    
    def mrisPreproc(self, **kwargs):
        cacheIn = self.measure + '.' + 'fwhm' + str(self.fwhm) + '.' + self.template
        argStr = '--cache-in ' + cacheIn
        mrispp = freesurfer.MRISPreproc(fsgd_file = self.fsgdFile.replace('.txt', '.fixed.txt'), target = self.template, hemi = self.hemi, out_file = self.data, args = argStr)
        if os.path.exists(self.data) and not self.overwrite:
            return
        else:
            mrispp.run()
    
    def mriGlmFit(self, **kwargs):
        self.glmdir = False ; contrasts = False 
        for i in kwargs.keys():
            if i == 'glmdir':
                self.glmdir = kwargs[i]
            elif i == 'contrasts':
                contrasts = kwargs[i]
                if not isinstance(contrasts, list):
                    print 'Contrasts MUST BE LIST TYPE....must exit'
                    sys.exit(0) 
        if not self.glmdir or not contrasts:
            print 'must specify the glmdir for output and CONTRASTS!!!......must exit'
            sys.exit(0)
        glm = freesurfer.GLMFit(in_file = self.data, contrast = contrasts, design = self.fsgdMat, surf = True, subject_id = self.template, hemi = self.hemi, cortex = True, glm_dir = self.glmdir)
        glm.run()
    
    def mriGlmFitSim(self, **kwargs):
        sim = False
        bothHemis = True
        for i in kwargs.keys():
            if i == 'sim':
                sim = kwargs[i]
                if not isinstance(sim, list):
                    print 'SIM must be set as LIST type! must exit'
                    sys.exit(0)  
            if i == 'oneHemi':
                bothHemis = False              
        if not sim:
            print 'must specify the sim type...must exit'
            sys.exit(0)
        simType = sim[0]
        if simType == 'mc-z':
            thr = str(sim[1]) ; simSign = sim[2] ; cwp = str(sim[3])
            opts= ['mri_glmfit-sim', '--glmdir', self.glmdir, '--cache', thr, simSign, '--cwpvalthresh', cwp]
            if bothHemis:
                opts.append('--2spaces')
        else:
            print 'sim type not supported...must exit.', simType
            sys.exit(0)
            #this appears to do the full correction
            #opts= ['mri_glmfit-sim', '--glmdir', self.glmdir, '--sim', simType, nIter, thr, simBase, '--sim-sign', simSign]
        if self.overwrite and os.path.exists(os.path.join(self.glmdir, os.path.basename(self.data))):
            os.remove(os.path.join(self.glmdir, os.path.basename(self.data)))
        if not os.path.exists(os.path.join(self.glmdir, os.path.basename(self.data))):
            os.symlink(self.data, os.path.join(self.glmdir, os.path.basename(self.data)))
        p = sp.Popen(opts)
        p.wait()

class extractFreeSurferStats():
    def __init__(self, subjects_dir, **kwargs):
        self.subjects_dir = subjects_dir
        #some defaults
        #the name of the database
        self.mydb = 'freesurfer'
        #the prefix to the table Name
        self.tblroot = 'f5_freesurfer_' + datetime.date.today().strftime("%d_%m_%Y")
        #the name of the ID column (i.e., "subject")
        self.id_col = 'idc'
        #A suffix to apply to all variable names (default is none)
        self.varSfix = ''
        ##connect to the mysql server
        if kwargs.has_key('mydb'):
            self.mydb = kwargs['mydb']
        self.con = mysqldb.connect(db = self.mydb, read_default_file=os.path.join(os.getenv("HOME"), '.my.cnf'))
        ##set up a cursor to use
        self.cursor = self.con.cursor()
        for i in kwargs.keys():
            if i == 'con':
                con = kwargs[i]
                self.cursor = con.cursor()
            elif i == 'tblroot':
                self.tblroot = kwargs[i]
            elif i == 'id_col':
                self.id_col = kwargs[i]
            elif i == 'varSfix':
                self.varSfix = kwargs[i]        
    
    def dumpAparcStats(self, subject):
        stats_csv = ''
        stats_file = ''
        #setup the table info
        tbl = self.tblroot + '_aparc_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        hemis = ['lh', 'rh']
        for hemi in hemis:
            stats_dir = os.path.join(self.subjects_dir, subject, 'stats')
            stats_file = open(os.path.join(stats_dir, hemi + '.aparc.stats'), 'r')
            stats_csv = csv.reader(stats_file, delimiter=' ')
            for row in stats_csv:
                if row[0] != '#':
                    for i in range(0, len(row)):
                        try:
                            row.index('')
                            row.remove('')
                        except:
                            continue
                    roi = row[0]
                    vol = float(row[3])
                    thickavg = row[4]
                    surfarea = float(row[2])
                    col_vol = hemi + '_' + roi + '_vol' + self.varSfix
                    col_thk = hemi + '_' + roi + '_thickavg' + self.varSfix
                    col_sua = hemi + '_' + roi + '_surfarea' + self.varSfix
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                    except:
                        pass
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_thk))
                    except:
                        pass
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_sua))
                    except:
                        pass
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_thk, thickavg, self.id_col, subject))
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_sua, surfarea, self.id_col, subject))
            stats_file.close()
            stats_file = open(os.path.join(stats_dir, hemi + '.aparc.stats'), 'r')
            stats_csv = csv.reader(stats_file, delimiter=',')
            for row in stats_csv:
                if len(row) == 5:
                    if row[1].strip() == 'MeanThickness':
                        roi = hemi + '_'  + row[1].strip()
                        vol = float(row[3])
                        col_vol = roi
                        try:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                        except:
                            pass
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
            stats_file.close()            

    def dumpAparcPialStats(self, subject):
        stats_csv = ''
        stats_file = ''
        #setup the table info
        tbl = self.tblroot + '_aparc_pial_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        hemis = ['lh', 'rh']
        for hemi in hemis:
            stats_dir = os.path.join(self.subjects_dir, subject, 'stats')
            stats_file = open(os.path.join(stats_dir, hemi + '.aparc.pial.stats'), 'r')
            stats_csv = csv.reader(stats_file, delimiter=' ')
            for row in stats_csv:
                if row[0] != '#':
                    for i in range(0, len(row)):
                        try:
                            row.index('')
                            row.remove('')
                        except:
                            continue
                    roi = row[0]
                    vol = float(row[3])
                    thickavg = row[4]
                    surfarea = float(row[2])
                    col_vol = hemi + '_' + roi + '_vol' + self.varSfix
                    col_thk = hemi + '_' + roi + '_thickavg' + self.varSfix
                    col_sua = hemi + '_' + roi + '_psurfarea' + self.varSfix
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                    except:
                        pass
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_thk))
                    except:
                        pass
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_sua))
                    except:
                        pass
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_thk, thickavg, self.id_col, subject))
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_sua, surfarea, self.id_col, subject))
            stats_file.close()
            stats_file = open(os.path.join(stats_dir, hemi + '.aparc.pial.stats'), 'r')
            stats_csv = csv.reader(stats_file, delimiter=',')
            for row in stats_csv:
                if len(row) == 5:
                    if row[1].strip() == 'MeanThickness':
                        roi = hemi + '_'  + row[1].strip()
                        vol = float(row[3])
                        col_vol = roi
                    elif row[1].strip() == 'PialSurfArea':
                        roi = hemi + '_'  + row[1].strip()
                        vol = float(row[3])
                        col_vol = roi
                    else:
                        continue
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                    except:
                        pass
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
            stats_file.close()            

    def dumpLobeStats(self, subject):
        stats_csv = ''
        stats_file = ''
        #setup the table info
        tbl = self.tblroot + '_lobes_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        hemis = ['lh', 'rh']
        for hemi in hemis:
            stats_dir = os.path.join(self.subjects_dir, subject, 'stats')
            stats_file = open(os.path.join(stats_dir, hemi + '.lobes.stats'), 'r')
            stats_csv = csv.reader(stats_file, delimiter=' ')
            for row in stats_csv:
                if row[0] != '#':
                    for i in range(0, len(row)):
                        try:
                            row.index('')
                            row.remove('')
                        except:
                            continue
                    roi = row[0]
                    vol = float(row[3])
                    thickavg = row[4]
                    surfarea = float(row[2])
                    col_vol = hemi + '_' + roi + '_vol' + self.varSfix
                    col_thk = hemi + '_' + roi + '_thickavg' + self.varSfix
                    col_sua = hemi + '_' + roi + '_surfarea' + self.varSfix
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                    except:
                        pass
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_thk))
                    except:
                        pass
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_sua))
                    except:
                        pass
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_thk, thickavg, self.id_col, subject))
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_sua, surfarea, self.id_col, subject))
            stats_file.close()

    def dumpBAStats(self, subject):
        #setup the table info
        stats_csv = ''
        stats_file = ''
        tbl = self.tblroot + '_ba_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
            #specify the file to read in
            hemis = ['lh', 'rh']
            for hemi in hemis:
                stats_dir =os.path.join(self.subjects_dir, subject, 'stats')
                stats_file = os.path.join(stats_dir, hemi + '.BA.stats')
                if not os.path.exists(stats_file):
                    stats_file = os.path.join(stats_dir, hemi + '.BA_exvivo.stats')
                    if not os.path.exists(stats_file):
                        #was break
                        continue
                stats_file = open(stats_file, 'r')
                stats_csv = csv.reader(stats_file, delimiter=' ')
                for row in stats_csv:
                    if row[0] != '#':
                        for i in range(0, len(row)):
                            try:
                                row.index('')
                                row.remove('')
                            except:
                                continue
                        roi = row[0]
                        vol = float(row[3])
                        thickavg = row[4]
                        surfarea = float(row[2])
                        col_vol = hemi + '_' + roi + '_vol' + self.varSfix
                        col_thk = hemi + '_' + roi + '_thickavg' + self.varSfix
                        col_sua = hemi + '_' + roi + '_surfarea' + self.varSfix
                        try:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                        except:
                            pass
                        try:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_thk))
                        except:
                            pass
                        try:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_sua))
                        except:
                            pass
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_thk, thickavg, self.id_col, subject))
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_sua, surfarea, self.id_col, subject))
            stats_file.close()
    
    def dumpAsegStats(self, subject):
        stats_csv = ''
        stats_file = ''
        tbl = self.tblroot + '_aseg_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        stats_dir =os.path.join(self.subjects_dir, subject, 'stats')
        stats_file = open(os.path.join(stats_dir, 'aseg.stats'), 'r')
        stats_csv = csv.reader(stats_file, delimiter=' ')
        for row in stats_csv:
            if row[0] != '#':
                for i in range(0, len(row)):
                    try:
                        row.index('')
                        row.remove('')
                    except:
                        continue
                roi = row[4].replace('-', '_')
                if roi.split('_')[0] == '3rd' or roi.split('_')[0] == '4th' or roi.split('_')[0] == '5th':
                    roi = roi.replace('3rd', 'Third')
                    roi = roi.replace('4th', 'Fourth')
                    roi = roi.replace('5th', 'Fifth')
                vol = float(row[3])
                col_vol = roi + '_vol' + self.varSfix
                try:
                    self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                except:
                    pass
                self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
        stats_file.close()
    
    def dumpWmparcStats(self, subject):
        stats_csv = ''
        stats_file = ''
        tbl = self.tblroot + '_wmparc_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s=\'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        stats_dir =os.path.join(self.subjects_dir, subject, 'stats')
        stats_file = open(os.path.join(stats_dir, 'wmparc.stats'), 'r')
        stats_csv = csv.reader(stats_file, delimiter=' ')
        for row in stats_csv:
            if row[0] != '#':
                for i in range(0, len(row)):
                    try:
                        row.index('')
                        row.remove('')
                    except:
                        continue
                roi = row[4].replace('-', '_')
                vol = float(row[3])
                col_vol = roi + '_vol' + self.varSfix
                if roi.startswith('wm_lh') or roi.startswith('wm_rh'):
                    if roi.split('_')[2] != 'S' and roi.split('_')[2] != 'G':
                        try:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                        except:
                            pass
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
        stats_file.close()
    
    def dumpTbvStats(self, subject):
        stats_csv = ''
        stats_file = ''
        tbl = self.tblroot + '_tbv_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s=\'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        stats_dir =os.path.join(self.subjects_dir, subject, 'stats')
        stats_file = open(os.path.join(stats_dir, 'aseg.stats'), 'r')
        stats_csv = csv.reader(stats_file, delimiter=',')   
        for row in stats_csv:
            if len(row) == 5:
                col_vol = False
                if row[1].strip() == 'lhCortexVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'rhCortexVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'CortexVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'lhCerebralWhiteMatterVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'rhCerebralWhiteMatterVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'CerebralWhiteMatterVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'SubCortGrayVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix                
                elif row[1].strip() == 'TotalGrayVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'SupraTentorialVol':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'SupraTentorialVolNotVent':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'SupraTentorialVolNotVentVox':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                elif row[1].strip() == 'eTIV':
                    roi = row[1].strip()    
                    vol = float(row[3])
                    col_vol = roi + self.varSfix
                if not col_vol:
                    continue
                try:
                    self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                except:
                    pass
                self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
        stats_file.close()
        hemis = ['lh', 'rh']
        for hemi in hemis:
            stats_file = open(os.path.join(stats_dir, hemi + '.aparc.stats'), 'r')
            stats_csv = csv.reader(stats_file, delimiter=',')   
            for row in stats_csv:
                if len(row) == 5:
                    col_vol = False
                    if row[1].strip() == 'NumVert':
                        roi = row[1].strip()    
                        vol = float(row[3])
                        col_vol = hemi + '_' + roi + self.varSfix
                    elif row[1].strip() == 'WhiteSurfArea':
                        roi = row[1].strip()    
                        vol = float(row[3])
                        col_vol = hemi + '_' + roi + self.varSfix
                    elif row[1].strip() == 'MeanThickness':
                        roi = row[1].strip()    
                        vol = float(row[3])
                        col_vol = hemi + '_' + roi + self.varSfix
                    if not col_vol:
                        continue
                    try:
                        self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                    except:
                        pass
                    self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
        stats_file.close()

    def dumpBrainStemStats(self, subject):
        stats_csv = ''
        stats_file = ''
        tbl = self.tblroot + '_brainstem_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        stats_dir =os.path.join(self.subjects_dir, subject, 'mri')
        stats_file = open(os.path.join(stats_dir, 'brainstemSsVolumes.v10.txt'), 'r')
        stats_csv = csv.reader(stats_file, delimiter=' ')
        for row in stats_csv:
            for i in range(0, len(row)):
                try:
                    row.index('')
                    row.remove('')
                except:
                    continue
            roi = row[0]
            vol = float(row[1])
            col_vol = roi + '_vol' + self.varSfix
            try:
                self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
            except:
                pass
            self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
        stats_file.close()

    def dumpHippSubfieldStats(self, subject):
        stats_csv = ''
        stats_file = ''
        #setup the table info
        tbl = self.tblroot + '_hippfields_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        hemis = ['lh', 'rh']
        for hemi in hemis:
            stats_dir = os.path.join(self.subjects_dir, subject, 'mri')
            stats_file = open(os.path.join(stats_dir, hemi + '.hippoSfVolumes-T1.v10.txt'), 'r')
            stats_csv = csv.reader(stats_file, delimiter=' ')
            for row in stats_csv:
                for i in range(0, len(row)):
                    try:
                        row.index('')
                        row.remove('')
                    except:
                        continue
                roi = row[0].replace('-', '_')
                vol = float(row[1])
                col_vol = hemi + '_' + roi + '_vol' + self.varSfix
                try:
                    self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                except:
                    pass
                self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
            stats_file.close()

    def dumpThalNucStats(self, subject):
        stats_csv = ''
        stats_file = ''
        #setup the table info
        tbl = self.tblroot + '_thalnuc_stats'
        try:
            self.cursor.execute("""create table %s (%s varchar(50))""" % (tbl, self.id_col))
        except:
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s = \'%s\'""" % (self.id_col, tbl, self.id_col, subject))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, self.id_col, subject))
        #specify the file to read in
        stats_dir = os.path.join(self.subjects_dir, subject, 'mri')
        stats_file = open(os.path.join(stats_dir, 'ThalamicNuclei.v10.T1.volumes.txt'), 'r')
        stats_csv = csv.reader(stats_file, delimiter=' ')
        for row in stats_csv:
            for i in range(0, len(row)):
                try:
                    row.index('')
                    row.remove('')
                except:
                    continue
            roi = row[0].replace('-', '_').replace(')', '').replace('(', '_')
            vol = float(row[1])
            col_vol = roi + '_vol' + self.varSfix
            #print(row)
            #print(col_vol, vol)
            try:
                self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
            except:
                pass
            self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, vol, self.id_col, subject))
        self.con.commit()
        stats_file.close()
