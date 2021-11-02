"""KNICR / GenR DTI Processing Utilities"""
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
import nipype.interfaces.camino as camino
from nipype.interfaces.base import Undefined
import subprocess as sp
import socket as socket
import datetime as datetime
import argparse as argparse
try:
    import MySQLdb as mysqldb
except:
    pass


__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "1.0"



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


def read_txt(txtFile):
    t = open(txtFile, 'r')
    raw_data = [] ; data = []
    csv_reader = csv.reader(t, delimiter=',')
    for line in csv_reader:
        raw_data.append(line[0])
    num_rows = len(raw_data)
    for i in range(0, num_rows):
        if not i == 0:
            data.append(raw_data[i])
    t.close()
    return data


def write_txt(data, txtFile):
    t = open(txtFile, 'w')
    if isinstance(data[0], list):
        csv_writer = csv.writer(t, delimiter=',')
        csv_writer.writerows(data)
    else:
        for i in data:
            t.write(i + '\n')
    t.close()


def read_fs_clut(**kwargs):
    """
    lut = color look up table.
        options are: lobes, wmparc, aparc, and aseg.
        
    returns dictionary of freesurfer color lookup table
    
    key = index
    value = ROI
    """
    fslut = False
    #standard parcelleation options. DO NOT CHANGE. these default to the FreeSruer color lookup table
    #for some reason, when an annotation is made (e.g., for the lobes), the index options overlap with the FS color lut.
    #thus, you WILL have misclassified ROIs
    standardParcOpts = ('wmparc', 'aparc', 'aseg')
    for i in kwargs.keys():
        if i == 'lut':
            if kwargs[i] == 'lobes':
                fslut = open('/Users/rmuetzel/lobes.ctab', 'r')
            elif kwargs[i].startswith(standardParcOpts):
                fslut = open(os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt'), 'r')
    if not fslut:   
        print 'Unknown parctype...cannot find lookup table.  Must exit'
        sys.exit(0)
    fslut_csv = csv.reader(fslut, delimiter=' ')
    fs_clut = {}
    print 'Reading FS color lookup table...',
    for line in fslut_csv:
        if len(line) == 0:
            continue
        elif line[0].startswith('#'):
            continue
        replace_ws = True
        while replace_ws:
            try:
                line.remove('')
            except:
                replace_ws = False
        roi = line[1] ; index = line[0]
        fs_clut[index] = roi
    print 'Loaded: ', len(fs_clut), ' FreeSurfer Labels.'
    fslut.close()
    return fs_clut


def compute_stat(iData, mask, **kwargs):
    """
    Compute a weighted average of a mask
    
    Two inputs:
    iData = input file that is to be averaged (i.e., FA map, DR map, etc...). This should be a nibabel loaded file to save on reads.
        e.g.:
        iDdata = nb.load(iDataFile).get_data()
    
    mask = area of input file included in average (i.e., a freesurfer ROI). String, path to file.
    
    optional args:
    
    stat = string   Available opts:
                    wavg    ---> weighted average. Particularly useful in non-binary masks (i.e., probtrack normalized distirbutions)
                    mean    ---> mean within the mask (same result as wavg if mask is binary)
                    median  ---> median within mask
                    stdev   ---> standard deviation within the mask
                    vol     ---> volume of the mask (will return volume and number of voxels)
    threshold = float
                the threshold for the mask
                default = 1
                    
    rescale = bool. True means the mask will be rescaled to be 0->1. I don't think this makes any difference? But maybe needed for some distributions?
    """
    #set up some defaults
    stat = 'wavg'
    rescale = False
    threshold = 1
    val = False
    #check for optional args
    for i in kwargs.keys():
        if i == 'stat':
            stat = kwargs[i]
        elif i == 'rescale':
            rescale = kwargs[i]
        elif i == 'threshold':
            threshold = kwargs[i]
    #load in the mask data
    mask_data = nb.load(mask).get_data()
    #rescale the mask to be 0-> if requested
    if rescale:
        mask_data_rescaled = False
        if np.max(mask_data) > 0:
            mask_data_rescaled = mask_data / np.max(mask_data)
        else:
            print 'WARNING: Cannot rescale mask image.....probably because the mask is empty...'
    #compute the stat
    if stat == 'wavg':
        try:
            if rescale and isinstance(mask_data_rescaled, np.ndarray):
                val = np.ma.average(np.ma.masked_where(mask_data<threshold, iData), weights=np.ma.masked_where(mask_data<threshold, mask_data_rescaled))
            else:
                val = np.ma.average(np.ma.masked_where(mask_data<threshold, iData), weights=np.ma.masked_where(mask_data<threshold, mask_data))
        except ZeroDivisionError:
            val = False
    elif stat == 'mean':
        val = np.ma.mean(np.ma.masked_where(mask_data<threshold, iData))
    elif stat == 'median':
        val = np.ma.median(np.ma.masked_where(mask_data<threshold, iData))
    elif stat == 'stdev':
        val = np.ma.std(np.ma.masked_where(mask_data<threshold, iData))
    elif stat == 'vol':
        if threshold == 1:
            threshold = 0.99
        nvox = np.count_nonzero(mask_data>threshold)
        mask_hd = nb.load(mask).get_header()
        x = float(mask_hd['pixdim'][1]) ;y = float(mask_hd['pixdim'][2]) ; z = float(mask_hd['pixdim'][3])
        vol = nvox * (x * y * z)
        val = [nvox, vol]
    if  not stat == 'vol':
        if np.isnan(val):
            print 'Warning....nan detected'
            val = False
    return val


class knicrDTIprep():
    def __init__(self, dtiDir, subject, **kwargs):
        """
        Class to preprocess KNICR MRI data using FSL tools
        
        bet:  Brain Extraction Tool.
        
        ecc: Eddy Current Correction (eddy_correct)
        
        rot_bvecs:  xfmrot from FreeSurfer to rotate the gradient table based on the ecc.log
        
        fit: dtifit to fit the diffusion tensor
        
        Initialization arguments:
            dtiDir: Where the dti images reside
                    Assumes the dwi is called: dti_idc_${subject}_2mm.nii.gz
            
            subject: the Subject ID 
        """
        #updated to help with bids
        self.dwi = False
        self.bidsSes = False
        self.bidsDerDir = False
        for i in kwargs.keys():
            if i == 'bidsSes':
                self.bidsSes = kwargs[i]
            if i == 'dwi':
                self.dwi = kwargs[i]
            if i == 'bidsDerDir':
                self.bidsDerDir = kwargs[i]
        self.subject = subject
        if self.bidsSes:
            self.subjDir = os.path.join(dtiDir, self.subject, self.bidsSes)
        else:
            self.subjDir = os.path.join(dtiDir, self.subject)
        if self.bidsDerDir:
            self.dtifitDir = os.path.join(self.bidsDerDir, self.subject, self.bidsSes, 'dmri')
        else:
            self.dtifitDir = os.path.join(self.subjDir, 'dmri')
        if not os.path.exists(self.dtifitDir):
            os.makedirs(self.dtifitDir)
        if not self.dwi:
            self.dwi = os.path.join(self.subjDir, 'dti_idc_' + str(subject) + '_2mm.nii.gz')
        self.data = os.path.join(self.dtifitDir, 'data.nii.gz')
        self.bvecs = os.path.join(self.dtifitDir, 'bvecs.ecc')
        self.bvals = os.path.join(self.dtifitDir, 'bvals')
        self.nodif_brain_mask = os.path.join(self.dtifitDir, 'nodif_brain_mask.nii.gz')
    
    def probe_pixdim(self):
        loaded_nii = nb.load(self.dwi)
        pixdims = loaded_nii.get_header()['pixdim']
        return (pixdims[1], pixdims[2], pixdims[3])
    
    def resample_dti(self, voxel_size):
        def resampleMriConvert(img, oFile, vs):
            mri_convert = freesurfer.MRIConvert(in_file=img, out_file=oFile, vox_size=vs)
            mri_convert.run()
            logFile = os.path.join(self.dtifitDir, 'resample.log')
            write_log(mri_convert.cmdline, logFile)

        if not isinstance(self.dwi, list):
            oFile = os.path.join(self.dtifitDir, os.path.basename(self.dwi).replace('.nii.gz', '_' + str(voxel_size[0]) + 'mm.nii.gz'))
            if os.path.exists(oFile):
                self.dwi = oFile
                return
            #resample the DTI data
            resampleMriConvert(self.dwi, oFile, voxel_size)
            #copy the bvals and bvecs....a lot of the other code relies on the name of the self.dwi object for these
            shutil.copyfile(self.dwi.replace('.nii.gz', '.bvec'), oFile.replace('.nii.gz', '.bvec'))
            shutil.copyfile(self.dwi.replace('.nii.gz', '.bval'), oFile.replace('.nii.gz', '.bval'))
            #point to the new files
            self.dwi = oFile
            #self.bvecs = oFile.replace('.nii.gz', '.bvec')
            #self.bvals = oFile.replace('.nii.gz', '.bval')
            

        else:
            for i in self.dwi:
                img = os.path.join(self.subjDir, i)
                oFile = img.replace('.nii.gz', '_' + str(voxel_size[0]) + 'mm.nii.gz')
                if not os.path.exists(oFile):
                    resampleMriConvert(img, oFile, voxel_size)
            dwiRs = []
            for i in self.dwi:
                i = i.replace('.nii.gz', '_' + str(voxel_size[0]) + 'mm.nii.gz')
                dwiRs.append(i)
            self.dwi = dwiRs    
    
    def mergeDwis(self, **kwargs):
        """
        Merge multiple 4D DWI volumes into a single 4D
        Also Merge the bvecs and bvals
        
        Two options that can be passed are 'bTable' and 'transpose'
            bTable = merge ; merge the outputs of dcm2nii
                   = filename ; copy this file to a standard naming convention for future use
            transpose = {'x': -1, 'y': 1, 'z': 1} ; transpose the final btable. This is a dictionary with 1/-1 depending on what you want to do for xyz
        """
        if not type(self.dwi) is list:
            print 'Cannot merge DWIs...not a list'
            return
        bTable = 'merge' ; bFlipInfo = False
        oBvec = os.path.join(self.subjDir,  'dti_all_merged.bvec')
        oBval = os.path.join(self.subjDir,  'dti_all_merged.bval')
        oDwi = os.path.join(self.subjDir,  'dti_all_merged.nii.gz')
        bValue = 1000
        #self.dwi = sorted(self.dwi)
        for i in kwargs.keys():
            if i == 'bTable':
                bTable = kwargs[i]
            elif i == 'transpose':
                bFlipInfo = kwargs[i]
            elif i == 'bValue':
                bValue = kwargs[i]
        def mergeBtables():
            bvecCode = {'x': 0, 'y': 1, 'z': 2}
            bvecData = {}
            bvalData = []
            for v in self.dwi:
                rootName = v.strip('.nii.gz')
                vBvec = open(os.path.join(self.subjDir, rootName + '.bvec'), 'r')
                vBval = open(os.path.join(self.subjDir, rootName + '.bval'), 'r')
                vBvecCsv = csv.reader(vBvec, delimiter=' ')
                vBvalCsv = csv.reader(vBval, delimiter=' ')
                tmpBvecData = []
                for bvec in vBvecCsv:
                    tmpBvecData.append(bvec)
                for b in bvecCode.keys():
                    if not bvecData.has_key(b):
                        bvecData[b] = tmpBvecData[bvecCode[b]]
                    else:
                        for i in tmpBvecData[bvecCode[b]]:
                            if not i == '':
                                bvecData[b].append(i)
                for bval in vBvalCsv:
                    for i in bval:
                        if not i == '':
                            bvalData.append(i)
                vBvec.close()
                vBval.close()
            print bValue
            if type(bValue) is int:
                oBvalFile =   open(oBval, 'w')
                bvalWriter = csv.writer(oBvalFile, delimiter=' ')
                bvalWriter.writerow(bvalData)
                oBvalFile.close()
            oBvecFile =   open(oBvec, 'w')
            bvecWriter = csv.writer(oBvecFile, delimiter=' ')
            for d in bvecCode.keys():
                bvecWriter.writerow(bvecData[d])
            oBvecFile.close()

        def transposeBtable(bTable):
            b = open(bTable, 'r')
            csvB = csv.reader(b, delimiter=' ')
            bVecData = {}
            bValData = []
            for i in csvB:
                if len(i) == 1:
                    continue
                removeWs = True
                while removeWs:
                    try:
                        i.remove('')
                    except:
                        removeWs = False
                x=float(i[0])*bFlipInfo['x']
                y=float(i[1])*bFlipInfo['y']
                z=float(i[2])*bFlipInfo['z']
                if x == 0:
                    x = 0
                if y == 0:
                    y = 0
                if z == 0:
                    z = 0
                if x == 0 and y == 0 and z == 0:
                    bValData.append(0)
                else:
                    bValData.append(bValue)
                if not bVecData.has_key('x'):
                    bVecData['x'] = [x]
                else:
                    bVecData['x'].append(x)	
                if not bVecData.has_key('y'):
                    bVecData['y'] = [y]
                else:
                    bVecData['y'].append(y)		
                if not bVecData.has_key('z'):
                    bVecData['z'] = [z]
                else:
                    bVecData['z'].append(z)
            oBvecTable = open(oBvec, 'w')
            bVecWrite = csv.writer(oBvecTable, delimiter=' ', lineterminator='\n')
            bVecWrite.writerow(bVecData['x'])
            bVecWrite.writerow(bVecData['y'])
            bVecWrite.writerow(bVecData['z'])
            oBvecTable.close()
            b.close()
            if type(bValue) is int:
                oBvalTable = open(oBval, 'w')
                bValWrite = csv.writer(oBvalTable, delimiter=' ', lineterminator='\n')
                bValWrite.writerow(bValData)
                oBvalTable.close()
        if bTable == 'merge':
            mergeBtables()
        else:
            try:
                if not os.path.exists(bTable):
                    print "cannot find input Btable: ", bTable
            except TypeError:
                print 'you must specify a file name...'
        if type(bValue) is str:
            try:
                if os.path.exists(bValue):
                    print 'B-value:  ',  bValue
                    print 'oBval file: ', oBval
                    print  'final bval file: ', self.bvals
                    print 'Copying..'
                    shutil.copyfile(bValue, oBval)
                    shutil.copyfile(bValue, self.bvals)
                else:
                    print 'You must specify an existing bval table or a single integer bvalue'
                    return
            except:
                print 'You must specify an existing bval table or a single integer bvalue'
                return
        if bFlipInfo:
            transposeBtable(bTable)
        else:
            print 'add in a copy of the original B file here.....'
        mergeList = []
        for v in self.dwi:
            print v
            mergeList.append(os.path.join(self.subjDir, v))
        print mergeList
        fslmerge = fsl.Merge(dimension='t', in_files=mergeList, merged_file=oDwi)
        fslmerge.run()
        logFile = oDwi.replace('.nii.gz', '.mergeLog.log')
        write_log(fslmerge.cmdline, logFile)
        self.dwi = oDwi
        shutil.copyfile(oBvec, os.path.join(self.dtifitDir, 'bvecs'))
        
        
    def bet(self, **kwargs):
        #a default of 0.25 seems to work well
        f = 0.25
        #ability to overwrite with this
        for i in kwargs.keys():
            if i == 'f':
                f = kwargs[i]
        #check to see if the mask needs to be created.
        if not os.path.exists(self.nodif_brain_mask):
            print 'Running BET on: ', self.dwi
            nodif_brain = self.nodif_brain_mask.replace('_mask.nii.gz', '.nii.gz')
            bet = fsl.BET(frac=f, in_file=self.dwi, mask=True, out_file=nodif_brain)
            bet.run()
            logFile = os.path.join(self.dtifitDir, 'bet.log')
            write_log(bet.cmdline, logFile)
    
    def ecc(self):
        if not os.path.exists(self.data):
            print 'Running eddy_correct on:', self.dwi
            eddy_correct = fsl.EddyCorrect(in_file=self.dwi, out_file=self.data, ref_num=0, terminal_output='stream', ignore_exception=True)
            eddy_correct.run()
            logFile = os.path.join(self.dtifitDir, 'eddy_correct.log')
            write_log(eddy_correct.cmdline, logFile)
    
    def rot_bvecs(self, **kwargs):
        """
        Rotate the gradient bvecs table using the ECC output from FSL.
        
        Default is to assume bvecs file is simply called, "bvecs". and that ECC params come from the "ecc" function in this class
        
        
        Optional :
        bvec = your_bvec_file to overide the assumed naming convention.
        This assumes that the bvals file has the same naming convention!
        """
        #rotate the bvecs based on the ECC output.    
        #xfmrot ecclog in_bvecs out_bvecs
        eccmat = os.path.join(self.dtifitDir, 'data.ecclog')
        bvec_raw = os.path.join(self.subjDir, 'dti_idc_' + str(self.subject) + '.bvec')
        for i in kwargs.keys():
            if i == 'bvec':
                bvec_raw = kwargs[i]
        if not os.path.exists(self.bvecs):
            print 'Rotating gradient table for eddy current correction: '
            print 'Input: ', bvec_raw
            print 'output: ', self.bvecs
            if not os.path.exists(self.bvals) and bvec_raw.endswith('.bvec'):
                shutil.copyfile(bvec_raw.replace('.bvec', '.bval'), self.bvals)
            optStr = ["xfmrot", eccmat, bvec_raw, self.bvecs]
            print(optStr)
            sp.call(optStr)

    def eddy(self, **kwargs):
        pe = [0, 1, 0] #x, y, z
        dwellTime = 0.087
        for i in kwargs.keys():
            if i == 'pe':
                if not type(kwargs[i]) == list or not len(kwargs[i]) == 3:
                    print('ERROR:  phase encode pe argument not set correctly')
                    return False
                pe = kwargs[i]
            elif i == 'dwellTime':
                dwellTime = kwargs[i]
        ndwis = nb.load(self.dwi).shape[3]
        eddyIndices = self.data.replace('.nii.gz', '.indices.txt')
        eddyAcqParm = self.data.replace('.nii.gz', '.acqparm.txt')
        def mkEddyFile(nvols, oFile, eddyFileType, pe, dwellTime):
            d = False
            if eddyFileType == 'Indices':
                d = [1]*ndwis
            elif eddyFileType == 'AcqParm':
                pe.append(dwellTime)
                d = pe
            f = open(oFile, 'w')
            csvWriter = csv.writer(f, delimiter = ' ')
            csvWriter.writerow(d)
            f.close()
        #check if eddy has been run
        if os.path.exists(eddyIndices):
            return
        mkEddyFile(ndwis, eddyIndices, 'Indices', pe, dwellTime)
        mkEddyFile(ndwis, eddyAcqParm, 'AcqParm', pe, dwellTime)
        bvec_raw = self.dwi.replace('.nii.gz', '.bvec')
        bval_raw = self.dwi.replace('.nii.gz', '.bval')
        optStr = ['eddy', '--imain=' + self.dwi, '--mask=' + self.nodif_brain_mask, '--bvecs=' + bvec_raw, '--bvals=' + bval_raw, '--out=' + self.data, '--index=' + eddyIndices, '--acqp=' + eddyAcqParm, '--repol', '-v']
        sp.call(optStr)
        if not os.path.exists(self.bvals):
            shutil.copyfile(bval_raw, self.bvals)
        shutil.copyfile(self.data+'.eddy_rotated_bvecs', self.bvecs)
        logFile = os.path.join(self.dtifitDir, 'eddy.log')
        write_log(str(optStr), logFile)        
    
    
    def fit(self):
        obn = os.path.join(self.dtifitDir, 'dtifit')
        fa = obn + '_FA.nii.gz'
        if not os.path.exists(fa):
            print 'Fitting the tensor: ', fa
            #in the current version of nipype...it doesn't handle the out_files properly, in terms of finding them
            #thus, the ignore_exception flag is activated..
            dtifit = fsl.DTIFit(base_name=obn, bvals=self.bvals, bvecs=self.bvecs, dwi=self.data, mask=self.nodif_brain_mask, sse=True, terminal_output='stream', ignore_exception=True)
            dtifit.run()
            logFile = obn + '.log'
            write_log(dtifit.cmdline, logFile)
            #Next make the axial diffusion scan...a simple copy of the primary eigen value
            iFile = obn + '_L1.nii.gz'
            oFile = obn + '_AD.nii.gz'
            fslmaths = fsl.ImageMaths(in_file=iFile, out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
            fslmaths.run()
            logFile = obn + '.log'
            write_log(fslmaths.cmdline, logFile)
            #next make the radial diffusivity scan...an average of the secondary and tertiary eigenvalues
            iFile = obn + '_L2.nii.gz'
            rd_op_str = '-add ' + obn + '_L3.nii.gz' + ' -div 2'
            oFile = obn + '_RD.nii.gz'
            fslmaths = fsl.ImageMaths(in_file=iFile, op_string=rd_op_str, out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
            fslmaths.run()
            logFile = obn + '.log'
            write_log(fslmaths.cmdline, logFile)
    
    def fit_restore(self, **kwargs):
        fmed = False
        for i in kwargs.keys():
            if i == 'fmed':
                fmed = kwargs[i]
        oDir = os.path.join(self.dtifitDir, 'camino_restore')
        if fmed:
            oDir = oDir + '_fmed'
        logFile = os.path.join(oDir, 'camino_restore.log')
        if not os.path.exists(oDir):
            os.makedirs(oDir)
        #for some reason, the masking within modelfit is failing. so, we mask the dwi first.
        dwi_masked = os.path.join(oDir, os.path.basename(self.data.replace('.nii.gz', '_masked.nii.gz')))
        if not os.path.exists(dwi_masked):
            fslmaths = fsl.maths.ApplyMask(in_file=self.data, mask_file=self.nodif_brain_mask, out_file=dwi_masked,  output_type='NIFTI_GZ')
            fslmaths.run()
            write_log(fslmaths.cmdline, logFile)
        if fmed:
            print 'Applying median filter to data...'
            dwi_masked_fmed = dwi_masked.replace('.nii.gz', '_fmed.nii.gz')
            fslmaths = fsl.SpatialFilter(in_file=dwi_masked, out_file=dwi_masked_fmed, operation='median', terminal_output='stream', kernel_shape='gauss', kernel_size=2.0, )
            fslmaths.run()
            dwi_masked = dwi_masked_fmed
            write_log(fslmaths.cmdline, logFile)
        bvec_camino = os.path.join(oDir, 'bvecs.ecc.camino')
        if not os.path.exists(bvec_camino):
            print 'Converting bvecs file: ', bvec_camino
            stout = open(bvec_camino, 'w')
            opts = ['fsl2scheme', '-bvecfile', self.bvecs, '-bvalfile', self.bvals]
            p = sp.Popen(opts, stdout=stout)
            p.wait()
            stout.close()
            write_log(str(opts), logFile)
        dwi_camino = os.path.join(oDir, 'data.Bfloat')
        if not os.path.exists(dwi_camino):
            print 'Converting dwi to camino format: ', dwi_camino
            img2vox = camino.Image2Voxel(in_file=dwi_masked, out_file=dwi_camino, out_type='float', terminal_output='stream', ignore_exception=True)
            img2vox.run()
            write_log(img2vox.cmdline, logFile)
        snr = os.path.join(oDir, 'estimatesnr.out')
        if not os.path.exists(snr):
            print 'computing sigma for RESTORE:'
            opts = ['estimatesnr', '-inputfile', dwi_camino, '-bgmask', self.nodif_brain_mask, '-schemefile', bvec_camino]
            #opts = ['estimatesnr', '-inputfile', dwi_camino, '-schemefile', bvec_camino]
            stout = open(snr, 'w')
            p = sp.Popen(opts, stdout=stout)
            p.wait()
            stout.close()
            write_log(str(opts), logFile)
        sigmaFile = open(snr, 'r')
        sigmaFilecsv = csv.reader(sigmaFile, delimiter='\t')
        for line in sigmaFilecsv:
            if len(line) > 1:
                if line[0] == 'sigma mult:':
                    sigma = float(line[2])
                    print 'sigma: ', sigma
        sigmaFile.close()
        restore_tensor = os.path.join(oDir, 'tensor.restore.Bdouble')
        #outlier_map = os.path.join(oDir, 'tensor.restore.outlier.Bdouble')
        #residmap = os.path.join(oDir, 'tensor.restore.residual.Bdouble')
        #t = open(outlier_map, 'w') ; t.close()
        #t = open(residmap, 'w') ; t.close()
        restore_tensor = os.path.join(oDir, 'tensor.restore.Bdouble')
        if not os.path.exists(restore_tensor):
            print 'Fitting restore model'
            #modelFit = camino.ModelFit(in_file=dwi_camino, model='restore', scheme_file=bvec_camino, bgmask=ero_mask, sigma=sigma, out_file=restore_tensor, ignore_exception=True)
            modelFit = camino.ModelFit(in_file=dwi_camino, model='restore', scheme_file=bvec_camino, bgthresh=1,sigma=sigma, out_file=restore_tensor, ignore_exception=True)
            modelFit.run()
            write_log(modelFit.cmdline, logFile)
        #restore_tensor_nii_pfix = os.path.join(oDir, 'tensor.restore_')
        #restore_tensor_nii = restore_tensor_nii_pfix + 'dt.nii.gz'
        #if not os.path.exists(restore_tensor_nii):
        #    print 'saving out tensor as nifti', restore_tensor_nii
        #    dt2nii = camino.DT2NIfTI(in_file = restore_tensor, output_root=restore_tensor_nii_pfix, header_file=dwi, terminal_output='stream', ignore_exception=True)
        #    dt2nii.run()
        restore_eigen = os.path.join(oDir, 'eigen.restore.Bdouble')
        if not os.path.exists(restore_eigen):
            print 'saving out eigen system'
            compute_eigen = camino.ComputeEigensystem(in_file=restore_tensor, out_file=restore_eigen, terminal_output='stream', ignore_exception=True)
            compute_eigen.run()
            write_log(compute_eigen.cmdline, logFile)
        scalars = {'FA':'fa', 'MD':'md', 'AD':'l1','RD':'rd', 'RA':'ra', 'L1':'l1', 'L2':'l2', 'L3':'l3'}
        for scalar in scalars.keys():
            met = scalars[scalar]
            if not fmed:
                oFile = os.path.join(oDir, 'camino_restore_' + scalar + '.nii.gz')
            else:
                oFile = os.path.join(oDir, 'camino_restore_fmed_' + scalar + '.nii.gz')
            if not os.path.exists(oFile):
                dtmetric = camino.DTMetric(eigen_data=restore_eigen, metric=met, data_header=self.data, outputfile=oFile,terminal_output='stream', ignore_exception=True)
                dtmetric.run()
                write_log(dtmetric.cmdline, logFile)
        restoreTensorNii = os.path.join(oDir, 'camino_restore_fmed_tensor.nii.gz')
        if not os.path.exists(restoreTensorNii):
            dt2nifti = camino.DT2NIfTI(in_file = restore_tensor, header_file = dwi_masked, output_root=restoreTensorNii.strip('.nii.gz'), ignore_exception=True, terminal_output='stream')
            dt2nifti.run()
        #cleanup
        os.remove(dwi_masked)
        os.remove(dwi_camino)
        os.remove(restore_eigen)
        os.remove(restore_tensor)
    
    def bedpost(self):
        bpxDir = self.dtifitDir + '.bedpostX'
        nfibers=2
        logFile = os.path.join(self.dtifitDir, 'bedpostX.log')
        #the rotated bvecs has '.ecc' appened...copy this to be bvecs for bedpost
        shutil.copyfile(self.bvecs, self.bvecs.replace('.ecc', ''))
        if not os.path.exists(bpxDir):
            try:
                bpx = fsl.BEDPOSTX(dwi=self.data, bvecs=self.bvecs, bvals=self.bvals, mask=self.nodif_brain_mask, fibres=nfibers, bpx_directory=self.dtifitDir, ignore_exception=True)
                print bpx.cmdline
                bpx.run()
                write_log(bpx.cmdline, logFile)
            except:
                print('cannot run bedpost via nipype interface...use bedpostCartesius')
                #cartesius chokes on the previous set because nipuype was updated...use the bedpostCartesius function
                #bpx = fsl.BEDPOSTX5(logdir=self.dtifitDir, out_dir=self.dtifitDir, terminal_output='stream', ignore_exception=True)
                #bpx.run()
                #write_log(bpx.cmdline, logFile)

    def bedpostCartesius(self, **kwargs):
        overWriteBpx = False
        for i in kwargs.keys():
            if i == 'overWriteBpx':
                overWriteBpx = kwargs[i]
            if i == 'cpBvecs':
                #this is handy if you skip eddy/dtifit, etc...
                cpBvecs = kwargs[i]
        bpxDir = self.dtifitDir + '.bedpostX'
        nfibers=2
        logFile = os.path.join(self.dtifitDir, 'bedpostX.log')
        #the rotated bvecs has '.ecc' appened...copy this to be bvecs for bedpost
        if cpBvecs:
            shutil.copyfile(self.bvecs, self.bvecs.replace('.ecc', ''))
        if not os.path.exists(bpxDir) or overWriteBpx:
            optStr = ['bedpostx', self.dtifitDir, '-n', str(nfibers)]
            sp.call(optStr)
            write_log(str(optStr), logFile)


class knicrDTIreg():
    def __init__(self, dtiDir, subject, **kwargs):
        self.searchAngle = [-180, 180]
        self.bidsSes = False
        for i in kwargs.keys():
            if i == 'searchAngle':
                self.searchAngle = kwargs[i]
            if i == 'bidsSes':
                self.bidsSes = kwargs[i]
        if self.bidsSes:
            subjDir = os.path.join(dtiDir, subject, self.bidsSes)
        else:
            subjDir = os.path.join(dtiDir, subject)
        self.subjDir = os.path.join(subjDir, 'dmri')
        self.regDir = os.path.join(subjDir, 'reg')
        self.labelDir = os.path.join(subjDir, 'label')
        self.subject = subject
        if not os.path.exists(self.subjDir):
            print 'dmri folder not present: ', self.subjDir
            print 'Please run knicrDTIpp FIRST!'
            return False
        if not os.path.exists(self.regDir):
            os.makedirs(self.regDir)
        if not os.path.exists(self.labelDir):
            os.makedirs(self.labelDir)

    
    def fsreg(self, subjects_dir):
        """
        fsreg:  Tool to align the dwi data to the FreeSurfer brain.mgz using FNIRT
        
        Prereq:  knicrDTIprep is run first and autorecon-all is run on the T1.
        
        Step 1: 6dof flirt dwi -> brain.mgz
        Step 2: FNIRT dwi -> brain.mgz
        Step 3: invert warp from step 2
        """
        #take the FA map (requires knicrDTIprep is run)
        fa = os.path.join(self.subjDir, 'dtifit_FA.nii.gz')
        #where the brain.mgz will be found later
        fsdir = os.path.join(subjects_dir, self.subject)
        #data prov. 
        logFile = os.path.join(self.regDir, 'fsreg.log')
        #make sure they have a freesurfer dataset
        if not os.path.exists(fsdir):
            print 'Cannot find subjects dir: ', fsdir
            return False
        #first convert brain.mgz -> brain.nii.gz for FSL
        iFile = os.path.join(fsdir, 'mri', 'brain.mgz')
        anat = os.path.join(self.regDir, 'fsbrain.nii.gz')
        if not os.path.exists(iFile):
            print 'Subject does not have brain.mgz', iFile
            return False
        if not os.path.exists(anat):
            mri_convert = freesurfer.MRIConvert(in_file=iFile, out_file=anat)
            mri_convert.run()
            write_log(mri_convert.cmdline, logFile)
        #because the BET isn't perfect on the FA map, we have a bright ring around the brain
        #do a few erosions to get rid of it, as it messes up the FNIRT
        fa_ero = os.path.join(self.regDir, 'dtifit_FA_ero.nii.gz')
        if not os.path.exists(fa_ero):
            fslmaths = fsl.ImageMaths(in_file=fa, op_string='-ero -ero', out_file=fa_ero, output_type='NIFTI_GZ', out_data_type='float')
            fslmaths.run()
            write_log(fslmaths.cmdline, logFile)
        ##do we need to do the "flip for fsl" crap?
        #next do the flirt, FA -> brain.nii.gz
        oFile = os.path.join(self.regDir, os.path.basename(fa_ero).replace('.nii.gz', '_to_fsbrain_lin.nii.gz'))
        mat = os.path.join(self.regDir, os.path.basename(fa_ero).replace('.nii.gz', '_to_fsbrain_lin.mat'))
        if not os.path.exists(oFile):
            print 'Registering FA to fs anatomical using FLIRT.'
            search = [-180, 180]
            flirt = fsl.FLIRT(in_file=fa_ero, reference=anat, out_file=oFile, out_matrix_file=mat, dof=6, searchr_x=search, searchr_y=search, searchr_z=search)
            flirt.run()
            write_log(flirt.cmdline, logFile)
        #next do the FNIRT, FA -> brain.nii.gz
        oFile = os.path.join(self.regDir, os.path.basename(fa_ero).replace('.nii.gz', '_to_fsbrain_nonlin.nii.gz'))
        log = oFile.replace('.nii.gz', '_log.txt')
        warp = oFile.replace('.nii.gz', '_warp.nii.gz')
        config='FA_2_FMRIB58_1mm'
        if not os.path.exists(warp):
            print 'Registering FA to fs anatomical using FNIRT.'
            fnirt = fsl.FNIRT(in_file=fa_ero, fieldcoeff_file=warp, warped_file=oFile, affine_file=mat, config_file=config, ref_file=anat, log_file=log)
            fnirt.run()
            write_log(fnirt.cmdline, logFile)
        #next do the invwarp to the warp file
        oFile = os.path.join(self.regDir, 'fsbrain_to_' + os.path.basename(fa_ero).replace('.nii.gz', '_nonlin.nii.gz'))    
        inverted_warp = oFile.replace('.nii.gz', '_warp.nii.gz')
        if not os.path.exists(inverted_warp):
            ##there is an error in the nipype.interfaces.fsl.utils file....the inverse_warp input spec is set to True, and should not be.
            print 'Inverting nonlinear warp.'
            invwarp = fsl.InvWarp(reference=fa_ero, warp=warp, inverse_warp=inverted_warp)
            invwarp.run()
            write_log(invwarp.cmdline, logFile)
        #last but not least, apply the warp to the T1, so we have a t1 image in dti space (to make sure the reg worked)
        if not os.path.exists(oFile):
            applywarp = fsl.ApplyWarp(in_file=anat, field_file=inverted_warp, out_file=oFile, ref_file=fa_ero, relwarp=True)
            applywarp.run()
            write_log(applywarp.cmdline, logFile)
    
    def mapFSLabel2Diff(self, subjects_dir, **kwargs):
        """
        Map freesurfer labels to diffusion data.
        
        ***Assumes the fs_reg function has been run!
        
        Will map the wmparc, aparc and aseg labels by default.
        
        Ranges are set for indices within the freesurfer color lut. 
        
        todo: Eventually build in kwargs for the following:
            indice ranges for wmparc, aparc and aseg, to specify a sub-set of labes.
            by default, it ignores the a2009s labels
            perhaps add support for other segmentations (i.e., ba labels)
        """
        #set a default parc type
        parcTypes = ['wmparc', 'aparc', 'aseg', 'lobes']
        for i in kwargs.keys():
            if i == 'parcTypes':
                parcTypes=kwargs[i]
        #specify the warp transform from anat->diff (from fs_reg)
        warp = os.path.join(self.regDir, 'fsbrain_to_dtifit_FA_ero_nonlin_warp.nii.gz')
        #point to the fa map (only for output resolution / dims)
        fa = os.path.join(self.regDir, 'dtifit_FA_ero.nii.gz')
        #do you want to ignore the -G and -S (gyral and sulcal labels...probably yes)
        ignore_a2009s = True
        sw = ('ctx-lh-S', 'ctx-lh-G','ctx-rh-S', 'ctx-rh-G', 'wm-lh-S', 'wm-lh-G', 'wm-rh-S', 'wm-rh-G')
        #which parcellations will be run
        #make sure the warp is present...if not, must quit
        if not os.path.exists(warp):
            print 'Cannot find nonlinear warp file..must exit..', warp
            return False
        #these are the indice ranges from the freesurfer color lut.
        #0-80 roughly = aseg ; 1000-3000 is aparc, and 3000-5000 is wmparc
        parc_range={'wmparc':[3000, 4999], 'aparc':[1000, 2999], 'aseg':[1, 60], 'lobes': [1000, 5003]}
        #keep provenance info
        logFile = os.path.join(self.labelDir, 'mapFSlabel2diff.log')
        #read the color look-up-table
        #we have short names "parcTypes" for each parcellation, and these have different file names
        for parcType in parcTypes:
            fs_clut = read_fs_clut(lut=parcType)
            print 'Mapping FreeSurfer labels to diffusion space: ', parcType,
            if parcType == 'wmparc':
                parcFile = 'wmparc'
            elif parcType == 'aparc':
                parcFile = 'aparc+aseg'
            elif parcType == 'aseg':
                parcFile = 'aparc+aseg'
            elif parcType == 'lobes':
                parcFile = 'lobes'
            else:
                print 'Unknown parcellation...must skip'
                continue
            #make an output folder for the data ; one for the 1mm FS data, and one for the DTI data
            anatLabelDir = os.path.join(self.labelDir, parcType, 'anat')
            diffLabelDir = os.path.join(self.labelDir, parcType, 'diff')
            if not os.path.exists(anatLabelDir):
                os.makedirs(anatLabelDir)
            if not os.path.exists(diffLabelDir):
                os.makedirs(diffLabelDir)
            #The original fs parcellation file
            parc_anat = os.path.join(subjects_dir, self.subject, 'mri', parcFile + '.mgz')
            #the new one in Nifti format
            parc = os.path.join(anatLabelDir, parcFile + '.nii.gz')
            #if the nifti one isn't mde yet, use mri_convert to generate it
            if not os.path.exists(parc):
                mri_convert = freesurfer.MRIConvert(in_file=parc_anat, out_file=parc)
                mri_convert.run()
                write_log(mri_convert.cmdline, logFile)
            #now loop over the index ranges you specified above
            for index in range(parc_range[parcType][0], parc_range[parcType][1]):
                index = str(index)
                #make sure it's actually in our dictionary
                if not fs_clut.has_key(index):
                    continue
                #get the ROI
                roi = fs_clut[index]
                #ignore the gyral and sulcal labels
                if ignore_a2009s:
                    if roi.startswith(sw):
                        continue
                #make an output mask in FS space
                anat_label = os.path.join(anatLabelDir, roi + '.nii.gz')
                if not os.path.exists(anat_label):
                    print '.',
                    op_str = '-thr ' + index + ' -uthr ' + index + ' -bin'
                    fslmaths = fsl.ImageMaths(in_file=parc, op_string=op_str, out_file=anat_label, output_type='NIFTI_GZ')
                    fslmaths.run()
                    write_log(fslmaths.cmdline, logFile)
                #now warp the mask from the previous step to diffusion space
                diff_label = os.path.join(diffLabelDir, roi + '.nii.gz')
                if not os.path.exists(diff_label):
                    applywarp = fsl.ApplyWarp(in_file=anat_label, field_file=warp, out_file=diff_label, ref_file=fa, relwarp=True, interp='nn')
                    applywarp.run()
                    write_log(applywarp.cmdline, logFile)
            print 'Done.'
    
    def stSpaceReg(self, template):
        fa = os.path.join(self.subjDir, 'dtifit_FA.nii.gz')
        fa_ero = os.path.join(self.regDir, 'dtifit_FA_ero.nii.gz')
        logFile = os.path.join(self.regDir, os.path.basename(template).replace('.nii.gz', '_registration.log'))
        if not os.path.exists(fa_ero):
            fslmaths = fsl.ImageMaths(in_file=fa, op_string='-ero -ero', out_file=fa_ero, output_type='NIFTI_GZ', out_data_type='float')
            fslmaths.run()
            write_log(fslmaths.cmdline, logFile)
        oFile = os.path.join(self.regDir, os.path.basename(fa_ero).replace('.nii.gz', '_to_') + os.path.basename(template).replace('.nii.gz', '_lin.nii.gz'))
        mat = oFile.replace('.nii.gz', '.mat')
        if not os.path.exists(oFile):
            print 'Registering FA map to standard space with FLIRT: ', template
            #search = [-180, 180]
            flirt = fsl.FLIRT(in_file=fa_ero, reference=template, out_file=oFile, out_matrix_file=mat, dof=12, searchr_x=self.searchAngle, searchr_y=self.searchAngle, searchr_z=self.searchAngle)
            flirt.run()
            write_log(flirt.cmdline, logFile)
        oFile = oFile.replace('_lin.nii.gz', '_nonlin.nii.gz')
        log = oFile.replace('.nii.gz', '.log')
        warp = oFile.replace('.nii.gz', '_warp.nii.gz')
        config='FA_2_FMRIB58_1mm'
        if not os.path.exists(warp):
            print 'Registering FA map to standard space with FNIRT: ', template
            fnirt = fsl.FNIRT(in_file=fa_ero, fieldcoeff_file=warp, warped_file=oFile, affine_file=mat, config_file=config, ref_file=template, log_file=log)
            fnirt.run()
            write_log(fnirt.cmdline, logFile)
        oFile = os.path.join(self.regDir, os.path.basename(template).replace('.nii.gz', '_to_') + os.path.basename(fa_ero).replace('.nii.gz', '_nonlin.nii.gz'))
        inverted_warp = oFile.replace('.nii.gz', '_warp.nii.gz')
        if not os.path.exists(inverted_warp):
            print 'Inverting nonlinear warp field.'
            invwarp = fsl.InvWarp(reference=fa_ero, warp=warp, inverse_warp=inverted_warp)
            invwarp.run()
            write_log(invwarp.cmdline, logFile)
        if not os.path.exists(oFile):
            applywarp = fsl.ApplyWarp(in_file=template, field_file=inverted_warp, out_file=oFile, ref_file=fa_ero, relwarp=True)
            applywarp.run()
            write_log(applywarp.cmdline, logFile)
    


class knicrAutoPtx():
    def __init__(self, dtiDir, subject, autoptxLib):
        self.subjDir = os.path.join(dtiDir, subject)
        self.regDir = os.path.join(self.subjDir, 'reg')
        self.autoptxDir = os.path.join(self.subjDir, 'autoptx')
        self.bpxDir = os.path.join(self.subjDir, 'dmri.bedpostX')
        self.autoptxLib = autoptxLib
        if not os.path.exists(self.bpxDir):
            return False
        if not os.path.exists(self.autoptxDir):
            os.makedirs(self.autoptxDir)
    
    def autoPtx_1_preproc(self):
        """
        Preprocess the diffusion data for automated tracking.
        
        The knicr version of this is actually done with knicrDTIprep() & knicrDTIreg()
        
        Here, we only make the 1mm native space refvol.
        """
        refvol = os.path.join(self.autoptxDir, 'refvol.nii.gz')
        if not os.path.exists(refvol):
            logFile = os.path.join(self.autoptxDir, 'autoptx_1_preproc.log')
            fa = os.path.join(self.regDir, 'dtifit_FA_ero.nii.gz')
            op_str = '-applyisoxfm 1'
            flirt = fsl.FLIRT(in_file=fa, reference=fa, out_file=refvol, args=op_str)
            flirt.run()
            write_log(flirt.cmdline, logFile)
            op_str = '-mul 0'
            fslmaths = fsl.ImageMaths(in_file=refvol, op_string=op_str, out_file=refvol, output_type='NIFTI_GZ')
            fslmaths.run()
            write_log(fslmaths.cmdline, logFile)
    
    def autoPtx_2_launchTractography(self):
        overwrite = False
        tractDir = 'tracts'
        logFile = os.path.join(self.autoptxDir, 'autoptx_2_launchtractography.log')
        refvol = os.path.join(self.autoptxDir, 'refvol.nii.gz')
        def read_autoPtxStructureList():
            structureList = open(os.path.join(self.autoptxLib, 'structureList'), 'r')
            structureInfo = {}
            structureCsv = csv.reader(structureList, delimiter=' ')
            for line in structureCsv:
                rmws = True
                while rmws:
                    try:
                        line.remove('')
                    except:
                        rmws = False
                struct,factor,walltime = line
                structureInfo[struct] = factor
            structureList.close()
            return structureInfo
        
        structInfo = read_autoPtxStructureList()
        print 'Running automated tractography of subject: ', self.subjDir
        for struct in structInfo.keys():
            nSeed = 1000
            nSeed = int(nSeed * float(structInfo[struct]))
            print 'Structure: ', struct, ' using ', nSeed, ' seeds per voxel.'
            masks = os.path.join(self.autoptxLib, 'protocols', struct)
            warp = os.path.join(self.regDir, 'FMRIB58_FA_1mm_to_dtifit_FA_ero_nonlin_warp.nii.gz')
            tracts = os.path.join(self.autoptxDir, 'tracts', struct)
            invert=False
            if overwrite and os.path.exists(tracts):
                    shutil.rmtree(tracts)
            if not os.path.exists(tracts):
                os.makedirs(tracts)
                os.makedirs(os.path.join(tracts, tractDir))
                def warp_masks(iFile, oFile):
                    if not os.path.exists(oFile):
                        applywarp = fsl.ApplyWarp(in_file=iFile, field_file=warp, out_file=oFile, ref_file=refvol, datatype='float')
                        applywarp.run()
                        write_log(applywarp.cmdline, logFile)
                        op_str = '-thr 0.1 -bin'
                        fslmaths = fsl.ImageMaths(in_file=oFile, op_string=op_str, out_file=oFile, output_type='NIFTI_GZ', out_data_type='char')
                        fslmaths.run()
                        write_log(fslmaths.cmdline, logFile)
                
                mask_list = ['seed', 'target', 'stop', 'exclude']
                for m in mask_list:
                    if m == 'seed':
                        seed_mni = os.path.join(masks, 'seed.nii.gz')
                        seed = os.path.join(tracts, 'seed.nii.gz')
                        warp_masks(seed_mni, seed)
                    elif m == 'target':
                        target_mni = os.path.join(masks, 'target.nii.gz')
                        target = os.path.join(tracts, 'target.nii.gz')
                        warp_masks(target_mni, target)
                    elif m == 'exclude':
                        exclude_mni = os.path.join(masks, 'exclude.nii.gz')
                        exclude = os.path.join(tracts, 'exclude.nii.gz')
                        warp_masks(exclude_mni, exclude)
                    elif m == 'stop':
                        stop_mni = os.path.join(masks, 'stop.nii.gz')
                        stop = Undefined
                        if os.path.exists(stop_mni):
                            stop = os.path.join(tracts, 'stop.nii.gz')
                            warp_masks(stop_mni, stop)
                invert = False
                if os.path.exists(os.path.join(masks, 'invert')):
                    invert = True
                ptx = fsl.ProbTrackX2(seed=seed, waypoints=target, samples_base_name=os.path.join(self.bpxDir, 'merged'), 
                mask=os.path.join(self.subjDir, 'dmri', 'nodif_brain_mask.nii.gz'), stop_mask=stop, n_samples=nSeed, opd=True, 
                out_dir=os.path.join(tracts, tractDir), avoid_mp=exclude, loop_check=True, force_dir=True, onewaycondition=invert, terminal_output='stream')
                ptx.run()
                write_log(ptx.cmdline, logFile)
                def read_waytotal(waytotal):
                    f = open(waytotal, 'r')
                    w = f.readline()
                    f.close()
                    return w
                
                way = float(read_waytotal(os.path.join(tracts, tractDir, 'waytotal')))                  
                op_str = '-div ' + str(way) + ' -range'
                #iFile = os.path.join(tracts,'tracts', 'density.nii.gz')
                iFile = os.path.join(tracts, tractDir, 'fdt_paths.nii.gz')
                oFile = os.path.join(tracts, tractDir, 'tractsNorm.nii.gz')
                fslmaths = fsl.ImageMaths(in_file=iFile, op_string=op_str, out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
    
    def autoPtx_2_launchTractographyDistance(self):
        overwrite = False
        tractDir = 'tractsDistance'
        logFile = os.path.join(self.autoptxDir, 'autoptx_2_launchtractographyDistance.log')
        refvol = os.path.join(self.autoptxDir, 'refvol.nii.gz')
        def read_autoPtxStructureList():
            structureList = open(os.path.join(self.autoptxLib, 'structureList'), 'r')
            structureInfo = {}
            structureCsv = csv.reader(structureList, delimiter=' ')
            for line in structureCsv:
                rmws = True
                while rmws:
                    try:
                        line.remove('')
                    except:
                        rmws = False
                struct,factor,walltime = line
                structureInfo[struct] = factor
            structureList.close()
            return structureInfo
        
        structInfo = read_autoPtxStructureList()
        print 'Running automated tractography with Distance Correction of subject: ', self.subjDir
        for struct in structInfo.keys():
            nSeed = 1000
            nSeed = int(nSeed * float(structInfo[struct]))
            print 'Structure: ', struct, ' using ', nSeed, ' seeds per voxel.'
            masks = os.path.join(self.autoptxLib, 'protocols', struct)
            warp = os.path.join(self.regDir, 'FMRIB58_FA_1mm_to_dtifit_FA_ero_nonlin_warp.nii.gz')
            tracts = os.path.join(self.autoptxDir, 'tracts', struct)
            invert=False
            if overwrite and os.path.exists(tracts):
                    shutil.rmtree(tracts)
            if not os.path.exists(os.path.join(tracts, tractDir)):
                #os.makedirs(tracts)
                os.makedirs(os.path.join(tracts, tractDir))
                def warp_masks(iFile, oFile):
                    if not os.path.exists(oFile):
                        applywarp = fsl.ApplyWarp(in_file=iFile, field_file=warp, out_file=oFile, ref_file=refvol, datatype='float')
                        applywarp.run()
                        write_log(applywarp.cmdline, logFile)
                        op_str = '-thr 0.1 -bin'
                        fslmaths = fsl.ImageMaths(in_file=oFile, op_string=op_str, out_file=oFile, output_type='NIFTI_GZ', out_data_type='char')
                        fslmaths.run()
                        write_log(fslmaths.cmdline, logFile)
                
                mask_list = ['seed', 'target', 'stop', 'exclude']
                for m in mask_list:
                    if m == 'seed':
                        seed_mni = os.path.join(masks, 'seed.nii.gz')
                        seed = os.path.join(tracts, 'seed.nii.gz')
                        warp_masks(seed_mni, seed)
                    elif m == 'target':
                        target_mni = os.path.join(masks, 'target.nii.gz')
                        target = os.path.join(tracts, 'target.nii.gz')
                        warp_masks(target_mni, target)
                    elif m == 'exclude':
                        exclude_mni = os.path.join(masks, 'exclude.nii.gz')
                        exclude = os.path.join(tracts, 'exclude.nii.gz')
                        warp_masks(exclude_mni, exclude)
                    elif m == 'stop':
                        stop_mni = os.path.join(masks, 'stop.nii.gz')
                        stop = Undefined
                        if os.path.exists(stop_mni):
                            stop = os.path.join(tracts, 'stop.nii.gz')
                            warp_masks(stop_mni, stop)
                invert = False
                if os.path.exists(os.path.join(masks, 'invert')):
                    invert = True
                ptx = fsl.ProbTrackX2(seed=seed, waypoints=target, samples_base_name=os.path.join(self.bpxDir, 'merged'), 
                mask=os.path.join(self.subjDir, 'dmri', 'nodif_brain_mask.nii.gz'), stop_mask=stop, n_samples=nSeed, opd=True, 
                out_dir=os.path.join(tracts, tractDir), avoid_mp=exclude, loop_check=True, force_dir=True, onewaycondition=invert, 
                correct_path_distribution=True, terminal_output='stream')
                ptx.run()
                write_log(ptx.cmdline, logFile)
                def read_waytotal(waytotal):
                    f = open(waytotal, 'r')
                    w = f.readline()
                    f.close()
                    return w
                
                way = float(read_waytotal(os.path.join(tracts, tractDir, 'waytotal')))                  
                op_str = '-div ' + str(way) + ' -range'
                #iFile = os.path.join(tracts,'tracts', 'density.nii.gz')
                iFile = os.path.join(tracts, tractDir, 'fdt_paths.nii.gz')
                oFile = os.path.join(tracts, tractDir, 'tractsNorm.nii.gz')
                fslmaths = fsl.ImageMaths(in_file=iFile, op_string=op_str, out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
                #now, convert the distance corrected path dist to a simple distance
                #the dist corrected path dist = distance * path distribution
                #To convert it, simply divide by the original path distribution (what you get when dist correction is disabled)
                iFile = os.path.join(tracts, tractDir, 'fdt_paths.nii.gz')
                opt_str = '-div ' + os.path.join(tracts, 'tracts', 'fdt_paths.nii.gz')
                oFile = os.path.join(tracts, tractDir, 'tractsDistance.nii.gz')
                fslmaths = fsl.ImageMaths(in_file=iFile, op_string=op_str, out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
    
    def autoPtx_to_Standard(self, **kwargs):
        overwrite = False
        tractDir = 'tracts'
        for i in kwargs.keys():
            if i == 'tractDir':
                tractDir = kwargs[i]
                print 'setting tractDir: ', tractDir
            elif i == 'overwrite':
                overwrite = True
                print 'setting overwite option: ', overwrite
        logFile = os.path.join(self.autoptxDir, 'autoptx_to_standard.log')
        refvol = os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'FMRIB58_FA_1mm.nii.gz')
        warp = os.path.join(self.regDir, 'dtifit_FA_ero_to_FMRIB58_FA_1mm_nonlin_warp.nii.gz')
        def read_autoPtxStructureList():
            structureList = open(os.path.join(self.autoptxLib, 'structureList'), 'r')
            structureInfo = {}
            structureCsv = csv.reader(structureList, delimiter=' ')
            for line in structureCsv:
                rmws = True
                while rmws:
                    try:
                        line.remove('')
                    except:
                        rmws = False
                struct,factor,walltime = line
                structureInfo[struct] = factor
            structureList.close()
            return structureInfo
        
        structInfo = read_autoPtxStructureList()
        print 'Warping tracts to FMRIB58 space...', self.subjDir
        for struct in structInfo.keys():
            if tractDir == 'tracts':
                tract = os.path.join(self.autoptxDir, 'tracts', struct, tractDir, 'tractsNorm.nii.gz')
            elif tractDir == 'tractsDistance':
                tract = os.path.join(self.autoptxDir, 'tracts', struct, tractDir, 'tractsDistance.nii.gz')
            oFile = tract.replace('.nii.gz', '_to_FMRIB58_FA_1mm.nii.gz')
            applywarp = fsl.ApplyWarp(in_file=tract, field_file=warp, out_file=oFile, ref_file=refvol, datatype='float')
            applywarp.run()
            write_log(applywarp.cmdline, logFile)
    
    def subDivideAtpxTracts(self, **kwargs):
        print 'Extracting Subdivisions for subject: ', self.subjDir
        overwrite = False
        subDivisions = {'acgc_l': 'cgc_l', 'acgc_r': 'cgc_r'}
        maskDir = '/Volumes/rbraid/xnat/nii/software'
        tractDir = 'tracts'
        for i in kwargs.keys():
            if i == 'subdivisions':
                subDivisions = kwargs[i]
            elif i == 'tractDir':
                tractDir = kwargs[i]
                print 'setting tractDir: ', tractDir
            elif i == 'overwrite':
                overwrite = True
                print 'setting overwite option: ', overwrite
            elif i == 'maskDir':
                maskDir = kwargs[i]
        warp = os.path.join(self.regDir, 'FMRIB58_FA_1mm_to_dtifit_FA_ero_nonlin_warp.nii.gz')
        refvol = os.path.join(self.autoptxDir, 'refvol.nii.gz')
        logFile = os.path.join(self.autoptxDir, 'autoptx_SubDivide.log')
        for subDivision in subDivisions.keys():
            struct = subDivisions[subDivision]
            iTract = os.path.join(self.autoptxDir, 'tracts', struct, tractDir, 'tractsNorm.nii.gz')
            if not os.path.exists(iTract):
                print "input tract does not exist...must exit"
                return
            if not os.path.exists(os.path.join(self.autoptxDir, 'tracts', subDivision, tractDir)):
                os.makedirs(os.path.join(self.autoptxDir, 'tracts', subDivision, tractDir))
            if subDivision.startswith('acgc'):
                iMask = os.path.join(maskDir, 'FMRIB58_FA_1mm_acgc_block_mask.nii.gz')
            oMask = os.path.join(self.autoptxDir, 'tracts', subDivision, tractDir, subDivision + '_block_mask.nii.gz')
            if not os.path.exists(oMask):
                applywarp = fsl.ApplyWarp(in_file=iMask, field_file=warp, out_file=oMask, ref_file=refvol, datatype='float')
                applywarp.run()
                write_log(applywarp.cmdline, logFile)
            oTract = os.path.join(self.autoptxDir, 'tracts', subDivision, tractDir, 'tractsNorm.nii.gz')
            if not os.path.exists(oTract):
                op_str = '-mul ' + oMask
                fslmaths = fsl.ImageMaths(in_file=iTract, op_string=op_str, out_file=oTract, output_type='NIFTI_GZ', out_data_type='float')
                fslmaths.run()
                write_log(fslmaths.cmdline, logFile)
    


class knicrTBSS():
    """
    knicrTBSS 
    
    A class of functions to run FSLs TBSS module in a stepwise fashion
    
    ...more to come...
    """
    def __init__(self, dtiDir, tbssDir, statsDir, template, **kwargs):
        #set some general paths to be used by all functions
        self.pfix = 'idc'
        self.dtiDir = dtiDir
        self.tbssDir = tbssDir
        self.faDir = os.path.join(self.tbssDir, 'FA')
        self.template = template
        self.interp = 'trilinear'
        if not os.path.exists(self.faDir):
            os.makedirs(self.faDir)
        self.statsDir = statsDir
        if self.tbssDir == self.statsDir:
            print 'Stats output folder CANNOT be the same as the TBSS dir...'
            print 'Please change the stats output folder...'
            sys.exit(0)
        if not os.path.exists(self.statsDir):
            os.makedirs(self.statsDir)
        for i in kwargs.keys():
            if i == 'interp':
                self.interp = kwargs[i]
            
    
    def tbss_1_preproc(self, subject, **kwargs):
        """
        tbss_1_preproc
        
        An adaptation from FSLs tbss_1_preproc: $FSLDIR/bin/tbss_1_preproc
        
        OPTIONAL ARGUMENTS:
            fa = The input FA map for registration
                [NLLS, fsl, knicr, camino_restore_fmed]
                The currently function two options are:
                            knicr (uses the fsl dtifit fa map)
                            camino_restore_fmed (median filtered CAMINO RESTORE FA map)
        
        !!!Changes from the FSL version!!!
        Some zero-padding is done in FSL, presumably to speed up the nonlinear registration. This is skipped here!
        """
        self.subject = subject
        self.nodif_brain_mask = os.path.join(self.dtiDir, subject, 'dmri', 'nodif_brain_mask.nii.gz')
        fa = False
        tracula_format = False
        tbss_preproc_default = False
        #check the optional arguments
        for i in kwargs.keys():
            if i == 'fa':
                if kwargs[i] == 'NLLS':
                    fa = os.path.join(self.dtiDir, subject, 'dmri', 'NLLS', 'dti_NLLS_fa.nii.gz')
                #because tracula seems to be doing some strange masking (with the aseg/wmparc)
                #we may need to compute the FA map again ourselves, with our own mask
                #this is because of the erosion step in tbss_1. 
                #The tracula brain_mask has some holes in the middle of the brain, where there are vessels.
                #this causes major problems with the erosion.  Further, the susceptibility artifact also causes masking problems for tbss.
                elif kwargs[i] == 'fsl':
                    fa = os.path.join(self.dtiDir, subject, 'dmri', 'fsl', 'dtifit_fsl_FA.nii.gz')
                elif kwargs[i] == 'knicr':
                    fa = os.path.join(self.dtiDir, subject, 'dmri', 'dtifit_FA.nii.gz')
                elif kwargs[i] == 'camino_restore_fmed':
                    fa = os.path.join(self.dtiDir, subject, 'dmri', 'camino_restore_fmed', 'camino_restore_fmed_FA.nii.gz')
                elif os.path.exists(kwargs[i]) and kwargs[i].endswith('.nii.gz'):
                    fa = kwargs[i]
            elif i == 'pfix':
                self.pfix = kwargs[i]
            elif i == 'tbss_preproc_default':
                tbss_preproc_default = kwargs[i]
        #if you didn't specify any FA map, assume tracula format as default.
        if not fa:
            tracula_format=True
        #make sure you can even find a tracula fa map
        if tracula_format:
            fa = os.path.join(self.dtiDir, subject, 'dmri', 'dtifit_FA.nii.gz')
        #if not, exit out.
        if not os.path.exists(str(fa)):
            print 'Cannot locate FA map: '
            print 'dtiDir: ', self.dtiDir
            print 'Subject: ', subject
            print 'fa: ', fa
            return False
        #begin the tbss_1 preprocessing
        print 'Processing: ', subject
        #first we need to erode the FA map a bit, to get rid of the high signal voxels around the edges
        logFile = os.path.join(self.faDir, self.pfix + '_' + subject + '.tbss_1_preproc.log')
        self.fa = os.path.join(self.faDir, self.pfix + '_' + subject + '_FA.nii.gz')
        if not os.path.exists(self.fa):
            #RLM Edited 21July2018
            #this erosion step causes issues with zeros inside of the brain (e.g., behind the splenium)
            #fix here not to erode the FA, but rather the input mask...
            if tbss_preproc_default:
                fslmaths = fsl.ImageMaths(in_file=fa, op_string='-ero', out_file=self.fa, output_type='NIFTI_GZ', out_data_type='float')
            else:
                fslmaths = fsl.ImageMaths(in_file=self.nodif_brain_mask, op_string='-ero -mul ' + fa, out_file=self.fa, output_type='NIFTI_GZ', out_data_type='float')
            fslmaths.run()
            write_log(fslmaths.cmdline, logFile)
        #make a mask of the eroded file now
        #this actually is a boundary mask for FLIRT input weighting
        self.mask = self.fa.replace('.nii.gz', '_mask.nii.gz')
        if not os.path.exists(self.mask):
            #there are zeros inside the brain in some places...this is not good.
            #first smooth the fa map, then do the erosion, then make the mask
            if tbss_preproc_default:
                fslmaths = fsl.ImageMaths(in_file=fa, op_string='-fmean -kernel gauss 4 -ero -bin', out_file=self.mask, output_type='NIFTI_GZ', out_data_type='float')
            else:
                fslmaths = fsl.ImageMaths(in_file=self.nodif_brain_mask, op_string='-ero -bin', out_file=self.mask, output_type='NIFTI_GZ', out_data_type='float')
            fslmaths.run()
            write_log(fslmaths.cmdline, logFile)
            op_str = '-dilD -dilD -sub 1 -abs -add ' + self.mask
            fslmaths = fsl.ImageMaths(in_file=self.mask, op_string=op_str, out_file=self.mask, output_type='NIFTI_GZ', out_data_type='char')
            fslmaths.run()
            write_log(fslmaths.cmdline, logFile)
        return True
    
    def tbss_2_reg(self, **kwargs):
        """
        tbss_2_reg
        
        Tool to run the flirt and fnirt steps of tbss.
        
        Note: we skip the computation of the mean deformation (i.e., _warp.msf, for determining the best fit from the all-to-all registration)
        
        This step must always follow the tbss_1_preproc....the fa map determined in step one is passed as self.fa.
        If tbss_1 has already been run once, and it's run again, it will not overwrite...
        """
        config='FA_2_FMRIB58_1mm'
        if not os.path.exists(self.fa):
            print 'Fa map not found: ', self.fa
            print 'run tbss_1_preproc first...'
            return
        obn = os.path.join(self.faDir, os.path.basename(self.fa).replace('.nii.gz', '') + '_to_' + os.path.basename(self.template).replace('.nii.gz', '') + '_lin')
        oFile = obn + '.nii.gz'
        mat = obn + '.mat'
        logFile = os.path.join(self.faDir, 'idc_' + self.subject + '.tbss_2_reg.log')
        if not os.path.exists(oFile):
            flirt = fsl.FLIRT(in_file=self.fa, reference=self.template, out_file=oFile, in_weight=self.mask, out_matrix_file=mat)
            flirt.run()
            write_log(flirt.cmdline, logFile)
        obn = os.path.join(self.faDir, os.path.basename(self.fa).replace('.nii.gz', '') + '_to_' + os.path.basename(self.template).replace('.nii.gz', '') + '_nonlin')
        self.warpedFile = obn + '.nii.gz'
        self.warp = obn + '_warp.nii.gz'
        log = obn + '_log.txt'
        if self.interp == 'trilinear':
            interp = 'linear'
        else:
            interp = self.interp
        if not os.path.exists(self.warp):
            args = '--interp=' + interp
            fnirt = fsl.FNIRT(in_file=self.fa, fieldcoeff_file=self.warp, warped_file=self.warpedFile, affine_file=mat, config_file=config, ref_file=self.template, log_file=log, args=args)
            fnirt.run()
            write_log(fnirt.cmdline, logFile)
    
    def tbss_3_postregA(self, **kwargs):
        """
        tbss_3_postregA (A= only applywarp step)
        
        This simply does the applywarp step from tbss_3.
        This skips some main steps...namely the fsl merge and the skeletonization. Those are done in two subsequent steps...
        """
        if not os.path.exists(self.warp):
            print 'Cannot find warp file: ', self.warp
            return
        applywarp = fsl.ApplyWarp(in_file=self.fa, field_file=self.warp, out_file=self.warpedFile, ref_file=self.template, relwarp=True, interp=self.interp)
        applywarp.run()
        logFile = os.path.join(self.faDir, 'idc_' + self.subject + '.tbss_3_postreg.log')
        write_log(applywarp.cmdline, logFile)
    
    def tbss_3_postregB(self, subjList, **kwargs):
        self.threshold = 0.2 #pass a kwarg to change it?
        #first read in the subject list if it's not already in txt format.
        if not isinstance(subjList, list):
            txtopts = ('.csv', '.txt')
            if subjList.endswith(txtopts):
                subjList = read_txt(subjList)
        self.subjList = subjList
        mergeList = []
        #go through the subjlist and make a list of file names to append
        for subj in self.subjList:
            fa = os.path.join(os.path.join(self.faDir, self.pfix + '_' + subj + '_FA_to_' + os.path.basename(self.template).replace('.nii.gz', '') + '_nonlin.nii.gz'))
            if os.path.exists(fa):
                mergeList.append(fa)
            else:
                #error check....if the fa map isn't there...the whole show stops...
                print 'Subject: ', subj, '  FA MAP!!'
                print fa
                print 'Cannot continue!!!'
                return False
        merged4d = os.path.join(self.statsDir, 'all_FA.nii.gz')
        print '4d: ', merged4d
        subjinfo = os.path.join(self.statsDir, 'tbss_3_postreg_mergeInfoFiles.csv')
        write_txt(mergeList, subjinfo)
        #for some reason, nipype throws an os error for too many args
        #this seems to work though, via sub process
        opts = ['fslmerge', '-t', merged4d]
        for i in mergeList:
            opts.append(i)
        sp.call(opts)
        #fslmerge = fsl.Merge(dimension='t', in_files=mergeList, merged_file=merged4d)
        #fslmerge.run()
        logFile = os.path.join(self.statsDir, 'tbss_3_postregB.log')
        write_log(str(opts), logFile)
        subjinfo = os.path.join(self.statsDir, 'tbss_3_postreg_mergeInfoSubjects.csv')
        write_txt(subjList, subjinfo)
        #next we make the mean FA image
        op_str = '-max 0 -Tmin -bin'
        fslmaths = fsl.ImageMaths(in_file=merged4d, op_string=op_str, out_file=os.path.join(self.statsDir, 'mean_FA_mask.nii.gz'), output_type='NIFTI_GZ', out_data_type='char')
        fslmaths.run()
        write_log(fslmaths.cmdline, logFile)
        op_str = '-mas ' + os.path.join(self.statsDir, 'mean_FA_mask.nii.gz')
        fslmaths = fsl.ImageMaths(in_file=merged4d, op_string=op_str, out_file=merged4d, output_type='NIFTI_GZ')
        fslmaths.run()
        write_log(fslmaths.cmdline, logFile)
        op_str = '-Tmean'
        fslmaths = fsl.ImageMaths(in_file=merged4d, op_string=op_str, out_file=os.path.join(self.statsDir, 'mean_FA.nii.gz'), output_type='NIFTI_GZ')
        fslmaths.run()
        write_log(fslmaths.cmdline, logFile)
        #now make the initial skeleton
        tbss_skeleton = fsl.TractSkeleton(in_file=os.path.join(self.statsDir, 'mean_FA.nii.gz'), skeleton_file=os.path.join(self.statsDir, 'mean_FA_skeleton.nii.gz'))
        tbss_skeleton.run()
        write_log(tbss_skeleton.cmdline, logFile)
        self.merged4d = merged4d
        return True
    
    def tbss_4_prestats(self, **kwargs):
        threshold = 0.2 #use a kwarg to edit this?
        logFile = os.path.join(self.statsDir, 'tbss_4_prestats.log')
        print 'Creating Skeleton Mask using threshold', str(threshold)
        op_str = '-thr ' + str(threshold) + ' -bin'
        fslmaths = fsl.ImageMaths(in_file=os.path.join(self.statsDir, 'mean_FA_skeleton.nii.gz'), op_string=op_str, out_file=os.path.join(self.statsDir, 'mean_FA_skeleton_mask.nii.gz'), output_type='NIFTI_GZ')
        fslmaths.run()
        write_log(fslmaths.cmdline, logFile)
        print 'Creating skeleton distancemap (for use in projection search)'
        op_str = '-mul 1 -add 1 -add ' + os.path.join(self.statsDir, 'mean_FA_skeleton_mask.nii.gz')
        fslmaths = fsl.ImageMaths(in_file=os.path.join(self.statsDir, 'mean_FA_mask.nii.gz'), op_string=op_str, out_file=os.path.join(self.statsDir, 'mean_FA_skeleton_mask_dst.nii.gz'), output_type='NIFTI_GZ')
        fslmaths.run()
        write_log(fslmaths.cmdline, logFile)
        distancemap = fsl.DistanceMap(in_file=os.path.join(self.statsDir, 'mean_FA_skeleton_mask_dst.nii.gz'), distance_map=os.path.join(self.statsDir, 'mean_FA_skeleton_mask_dst.nii.gz'))
        distancemap.run()
        write_log(distancemap.cmdline, logFile)
        #do a quick check on the cingulum file...make sure it's the right size!
        nb_template = nb.load(self.template)
        x = abs(nb_template.get_affine()[0][0])
        y = abs(nb_template.get_affine()[1][1])
        z = abs(nb_template.get_affine()[2][2])
        cing_file1mm = os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'LowerCingulum_1mm.nii.gz')
        if x == 1 and y == 1 and z == 1:
            cing_file = cing_file1mm
        elif x == 2 and y == 2 and z == 2:
            cing_file = os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'LowerCingulum_2mm.nii.gz')
            if not os.path.exists(cing_file):
                cing_file = os.path.join(self.statsDir, 'LowerCingulum_2mm.nii.gz')
                flirt = fsl.FLIRT(in_file=cing_file1mm, reference=os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'MNI152_T1_2mm.nii.gz'), out_file=cing_file, args='-applyisoxfm 2')
                flirt.run()
                write_log(flirt.cmdline, logFile)
        print 'Projecting all FA data onto Skeleton'
        final4d = self.merged4d.replace('.nii.gz', '_skeletonized.nii.gz')
        tbss_skeleton = fsl.TractSkeleton(threshold=threshold, distance_map=os.path.join(self.statsDir, 'mean_FA_skeleton_mask_dst.nii.gz'), 
        data_file=self.merged4d, projected_data=final4d, in_file=os.path.join(self.statsDir, 'mean_FA.nii.gz'), project_data=True)
        tbss_skeleton.inputs.use_cingulum_mask = Undefined
        tbss_skeleton.inputs.search_mask_file=cing_file
        tbss_skeleton.run()
        write_log(tbss_skeleton.cmdline, logFile)
    
    def tbss_nonFA_reg(self, subject, nonFA, **kwargs):
        """
        Required arguments:
            subject -- id number of subject
                string type
            nonFA --- this is the scalar (FA, MD, RD, AD)
                string type
        
        optional arguments:
            fit_method --- default is dtifit. data reside in dir/subj/dmri/dtifit_FA.nii.gz
            
            supported fit methods are:
                dipy_nlls (nonlinear least squares)
                dtifit_wls (linear weighted least squares)
            Assume that fit method folder containing scalars is within the dmri folder:
                dir/subj/dmri/fit_method/fit_method_FA.nii.gz
            string type.
                
        """
        fit_method = 'dtifit'
        for i in kwargs.keys():
            if i == 'fit_method':
                fit_method = kwargs[i]
        if fit_method == 'dtifit':
            iFile = os.path.join(self.dtiDir, subject, 'dmri', fit_method + '_' + nonFA + '.nii.gz')
        else:
            iFile = os.path.join(self.dtiDir, subject, 'dmri', fit_method, fit_method + '_' + nonFA + '.nii.gz')
        #do the reg for each subject here
        warp = os.path.join(self.faDir, self.pfix + '_' + subject + '_FA_to_' + os.path.basename(self.template).replace('.nii.gz', '') + '_nonlin_warp.nii.gz')
        if not os.path.exists(warp):
            print 'Cannot find warp file: ', warp
            sys.exit(0)
            return False
        oFile = os.path.join(self.faDir, self.pfix + '_' + subject + '_' + fit_method + '_' + nonFA + '_to_' + os.path.basename(self.template).replace('.nii.gz', '') + '_nonlin.nii.gz')
        if not os.path.exists(oFile):
            applywarp = fsl.ApplyWarp(in_file=iFile, field_file=warp, out_file=oFile, ref_file=self.template, relwarp=True, interp=self.interp)
            applywarp.run()
            logFile = os.path.join(self.faDir, 'idc_' + subject + '.tbss_nonFA.log')
            write_log(applywarp.cmdline, logFile)
    
    def tbss_nonFA_merge(self, nonFA, subjList):
        threshold=0.2
        #first read in the subject list if it's not already in txt format.
        if not isinstance(subjList, list):
            txtopts = ('.csv', '.txt')
            if subjList.endswith(txtopts):
                subjList = read_txt(subjList)
        self.subjList = subjList
        mergeList = []
        #go through the subjlist and make a list of file names to append
        for subj in subjList:
            scalar = os.path.join(self.faDir, self.pfix + '_' + subj + '_' + nonFA + '_to_' + os.path.basename(self.template).replace('.nii.gz', '') + '_nonlin.nii.gz')
            if os.path.exists(scalar):
                mergeList.append(scalar)
            else:
                #error check....if the fa map isn't there...the whole show stops...
                print 'Subject: ', subj, '  Cannot locate nonFA MAP!!'
                print scalar
                print 'Cannot continue!!!'
                return False
        merged4d = os.path.join(self.statsDir, 'all_' + nonFA + '.nii.gz')
        print '4d: ', merged4d
        #mer
        #for some reason, nipype throws an os error for too many args
        #this seems to work though, via sub process
        opts = ['fslmerge', '-t', merged4d]
        for i in mergeList:
            opts.append(i)
        sp.call(opts)
        #fslmerge = fsl.Merge(dimension='t', in_files=mergeList, merged_file=merged4d)
        #fslmerge.run()
        logFile = os.path.join(self.statsDir, 'tbss_nonFA.log')
        #write_log(fslmerge.cmdline, logFile)
        write_log(str(opts), logFile)
        subjinfo = os.path.join(self.statsDir, 'tbss_nonFA_' + nonFA + '_mergeInfo.csv')
        write_txt(subjList, subjinfo)
        #mask it
        op_str = '-mas ' + os.path.join(self.statsDir, 'mean_FA_mask.nii.gz')
        fslmaths = fsl.ImageMaths(in_file=merged4d, op_string=op_str, out_file=merged4d, output_type='NIFTI_GZ')
        fslmaths.run()
        logFile = os.path.join(self.statsDir, 'tbss_nonFA.log')
        write_log(fslmaths.cmdline, logFile)
        #skeletenize
        nb_template = nb.load(self.template)
        x = abs(nb_template.get_affine()[0][0])
        y = abs(nb_template.get_affine()[1][1])
        z = abs(nb_template.get_affine()[2][2])
        cing_file1mm = os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'LowerCingulum_1mm.nii.gz')
        if x == 1 and y == 1 and z == 1:
            cing_file = cing_file1mm
        elif x == 2 and y == 2 and z == 2:
            cing_file = os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'LowerCingulum_2mm.nii.gz')
            if not os.path.exists(cing_file):
                cing_file = os.path.join(self.statsDir, 'LowerCingulum_2mm.nii.gz')
                flirt = fsl.FLIRT(in_file=cing_file1mm, reference=os.path.join(os.environ['FSLDIR'], 'data', 'standard', 'MNI152_T1_2mm.nii.gz'), out_file=cing_file, args='-applyisoxfm 2')
                flirt.run()
                write_log(flirt.cmdline, logFile)
        else:
            cing_file=False
        if not cing_file:
            print 'Cannot find cingulum file.'
            return
        print 'Projecting all FA data onto Skeleton'
        final4d = merged4d.replace('.nii.gz', '_skeletonized.nii.gz')
        #run the tbss skeletonize thing
        tbss_skeleton = fsl.TractSkeleton(threshold=threshold, distance_map=os.path.join(self.statsDir, 'mean_FA_skeleton_mask_dst.nii.gz'), alt_data_file=merged4d, 
        data_file=os.path.join(self.statsDir, 'all_FA.nii.gz'), projected_data=final4d, in_file=os.path.join(self.statsDir, 'mean_FA.nii.gz'), project_data=True)
        tbss_skeleton.inputs.use_cingulum_mask = Undefined
        tbss_skeleton.inputs.search_mask_file=cing_file
        tbss_skeleton.run()
        write_log(tbss_skeleton.cmdline, logFile)
    


class knicrDTIstats():
    def __init__(self, mydb, dtiDir):
        print 'Negotiating connection to mysql DB...'
        ##connect to the mysql server
        self.con = mysqldb.connect(db = mydb, read_default_file=os.path.join(os.getenv("HOME"), '.my.cnf'))
        ##set up a cursor to use
        self.cursor = self.con.cursor()
        self.dtiDir = dtiDir
    
    def extract_fs_stats(self, tblroot, subjid, **kwargs):
        """
        extract DTI values from FreeSurfer ROIs
        
        *Assumes knicrDTIreg has been run, and the labels have been mapped to native space.
        
        Inputs: tblroot, which is the prefix to the table name that will be created
                This will be, tblroot_wmparc_dti_stats
                
                subjid  the id number of the subject.
        
        Optional inputs:    parcType    ---> default is wmparc
                                        ---> string format
                                        ---> Othersupported options: aseg, aparc, lobes
                            
                            scalarType  ---> default is dtifit from fsl
                                        ---> string format
                                        ---> The scalar maps for anything else must be in this format:
                                                dtiDir/subj/dmri/scalarType/scalarType_FA.nii.gz
                            
                            scalars     ---> default = [FA, MD, RD, AD]
                                        ---> list type
                            
                            stats       ---> default = [mean]
                                        ---> list type
                                        
                            
        """
        parcType = 'wmparc'
        scalarType = 'dtifit'
        stats = ['mean', 'vol']
        scalars = ['FA', 'MD', 'AD', 'RD']
        for i in kwargs.keys():
            if i == 'parcType':
                parcType = kwargs[i]
                print 'setting parcType to: ', parcType
            elif i == 'scalarType':
                scalarType = kwargs[i]
                print 'Setting scalar type to: ', scalarType
            elif i == 'scalars':
                if isinstance(kwargs[i], list):
                    scalars = kwargs[i]
                else:
                    print 'Cannot reset: ', i
                    print 'Must pass in list format...using default list: ', scalars
            elif i == 'stats':
                if isinstance(kwargs[i], list):
                    stats = kwargs[i]
                else:
                    print 'Cannot reset: ', i
                    print 'Must pass in list format...using default list: ', stats
        tbl = tblroot + '_' + parcType + '_dti_stats'
        id_col = 'ID'
        idval = subjid
        try:
            self.cursor.execute("""create table %s (%s varchar(10))""" % (tbl, id_col))
        except (mysqldb._mysql.OperationalError):
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s=\'%s\'""" % (id_col, tbl, id_col, idval))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, id_col, idval))
        #now read in the data
        subjDir = os.path.join(self.dtiDir, subjid)
        if scalarType == 'dtifit':
            scalarDir = os.path.join(subjDir, 'dmri')
        else:
            scalarDir = os.path.join(subjDir, 'dmri', scalarType)
        #make sure the file exists
        if not os.path.exists(scalarDir):
            print 'Cannot find scalar maps: ', scalarDir
            return False
        #now find the FS masks
        labelDir = os.path.join(subjDir, 'label', parcType, 'diff')
        if os.path.exists(labelDir):
            labels = sorted(os.listdir(labelDir))
        else:
            print 'Cannot find labels: ', labelDir
            return False
        #iterate through each type of DTI data
        run_vol_stat = True
        for scalar in scalars:
            #load the map in ahead of time
            scalar_map = os.path.join(scalarDir, scalarType + '_' + scalar + '.nii.gz')
            data = nb.load(scalar_map).get_data()
            #loop over each mask
            for label in labels:
                #we need to fix the name to make it mysql-friendly
                roi = label.replace('.nii.gz', '').replace('-', '_')
                if roi.split('_')[0] == '3rd' or roi.split('_')[0] == '4th' or roi.split('_')[0] == '5th':
                    roi = roi.replace('3rd', 'Third')
                    roi = roi.replace('4th', 'Fourth')
                    roi = roi.replace('5th', 'Fifth')
                #now, we have some stat options (mean, weighted-average, etc)
                for stat in stats:
                    #compute the value/statistic
                    val = compute_stat(data, os.path.join(labelDir, label), stat=stat)
                    if not val:
                        #if the mask is empty...just skip it from here
                        continue
                    #set up a column name for hte mysql db
                    if not stat == 'vol':
                        col = roi + '_dti_' + stat + '_' + scalar 
                        try:
                            #see if that column exists...if not...make it.
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col))
                        except:
                            pass
                        #update the data base with the data value...
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col, val, id_col, idval))
                    else:
                        if not run_vol_stat:
                            continue
                        vol_types = ['vol', 'nvox']
                        for v in vol_types:
                            col = roi + '_dti_' + v
                            try:
                                #see if that column exists...if not...make it.
                                self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col))
                            except:
                                pass
                            if v == 'nvox':
                                vol_val = val[0]
                            elif v == 'vol':
                                vol_val = val[1]
                            else:
                                continue
                            self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col, vol_val, id_col, idval))
                            
                        
                        
    
    def extract_autoptx_stats(self, tblroot, subjid, **kwargs):
        """
        extract DTI values from autoptx tractss
        
        *Assumes knicrAutoPtx has been run, 
        
        Inputs: tblroot, which is the prefix to the table name that will be created
                This will be, tblroot_wmparc_dti_stats
                
                subjid  the id number of the subject.
                
        Optional inputs:   
                            scalarType  ---> default is dtifit from fsl
                                        ---> string format
                                        ---> The scalar maps for anything else must be in this format:
                                                dtiDir/subj/dmri/scalarType/scalarType_FA.nii.gz
                                                
                            scalars     ---> default = [FA, MD, RD, AD]
                                        ---> list type
                                        
                            stats       ---> default = [mean, vol], which is the average of the normalized path distr. and the volume of the normalized dist.
                                        ---> list type
                                        ---> options: wavg, mean, median, stdev, vol
                            
                            structList =    default from autoptx (April 2014)
                                            dictionary type. key = structure, value = threshold.
                                            
        """
        scalarType = 'dtifit'
        stats = ['mean', 'vol']
        scalars = ['FA', 'MD', 'AD', 'RD']
        structList = {'ar_l':False,'ar_r':False,'atr_l': 0.002,'atr_r': 0.002,'cgc_l': 0.01,'cgc_r': 0.01,'cgh_l':0.02,'cgh_r':0.02,'cst_l':0.005,'cst_r':0.005,'fma':0.005,'fmi':0.01,'ifo_l':0.01,'ifo_r':0.01,
        'ilf_l':0.005,'ilf_r': 0.005,'mcp':0.0001,'ml_l':0.005,'ml_r':0.005,'ptr_l':0.005,'ptr_r':0.005,'slf_l':0.001,'slf_r':0.001,'str_l':0.005,'str_r':0.005,'unc_l':0.01,'unc_r':0.01}
        var_sfix = ''
        for i in kwargs.keys():
            if i == 'scalarType':
                scalarType = kwargs[i]
            elif i == 'scalars':
                if isinstance(kwargs[i], list):
                    scalars = kwargs[i]
                else:
                    print 'Cannot reset: ', i
                    print 'Must pass in list format...using default list: ', scalars
            elif i == 'stats':
                if isinstance(kwargs[i], list):
                    stats = kwargs[i]
                else:
                    print 'Cannot reset: ', i
                    print 'Must pass in list format...using default list: ', stats
            elif i == 'structList':
                if isinstance(kwargs[i], dict):
                    structList = kwargs[i]
                else:
                    print 'Cannot reset: ', i
                    print 'must pass list format....default will be used: ',
            elif i == 'var_sfix':
                 var_sfix = '_' + kwargs[i]
        tbl = tblroot + '_autoPtx_dti_stats'
        id_col = 'ID'
        idval = subjid
        try:
            self.cursor.execute("""create table %s (%s varchar(10))""" % (tbl, id_col))
            #indexing the ID column seems to help with slow downs once the table grows in size substantially
            self.cursor.execute("""ALTER TABLE %s ADD INDEX (%s)""" % (tbl, id_col))
        except (mysqldb._mysql.OperationalError):
            pass
        #check the table to see if the subject exists
        self.cursor.execute("""select %s from %s where %s=\'%s\'""" % (id_col, tbl, id_col, idval))
        pn_exist = self.cursor.fetchone()
        if not pn_exist:
            self.cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, id_col, idval))
        #now read in the data
        subjDir = os.path.join(self.dtiDir, subjid)
        if scalarType == 'dtifit':
            scalarDir = os.path.join(subjDir, 'dmri')
        else:
            scalarDir = os.path.join(subjDir, 'dmri', scalarType)
        #make sure the file exists
        if not os.path.exists(scalarDir):
            print 'Cannot find scalar maps: ', scalarDir
            return False
        #now find the FS masks
        autoptxDir = os.path.join(subjDir, 'autoptx', 'tracts')
        if not os.path.exists(autoptxDir):
            print 'Cannot find labels: ', autoptxDir
            return False
        #iterate through each type of DTI data
        for scalar in scalars:
            #load the map in ahead of time
            scalar_map = os.path.join(scalarDir, scalarType + '_' + scalar + '.nii.gz')
            data = nb.load(scalar_map).get_data()
            #loop over each tract
            for tract in structList.keys():
                normPath = os.path.join(autoptxDir, tract, 'tracts', 'tractsNormNative.nii.gz')
                if not os.path.exists(normPath) and os.path.exists(normPath.replace('Native.nii.gz', '.nii.gz')):
                    #this only works if the input data are isotropic...fixed it with below
                    #flirt = fsl.FLIRT(in_file=normPath.replace('Native.nii.gz', '.nii.gz'), reference=scalar_map, out_file=normPath, args='-applyisoxfm 2')
                    isoMat = 'isoXfm.mat'
                    if not os.path.exists(isoMat):
                        xfm = open(isoMat, 'w')
                        xfmWrite = csv.writer(xfm, delimiter=' ')
                        xfmWrite.writerows([[1,  0,  0,  0], [0,  1,  0,  0], [0,  0,  1,  0],[0,  0,  0,  1]])
                        xfm.close()
                    flirt = fsl.FLIRT(in_file=normPath.replace('Native.nii.gz', '.nii.gz'), reference=scalar_map, out_file=normPath, apply_xfm=True, in_matrix_file=isoMat)
                    flirt.run()
                    logFile = os.path.dirname(normPath) + 'downsampleNormPath.log'
                    write_log(flirt.cmdline, logFile)
                if not os.path.exists(normPath):
                    continue
                #now, we have some stat options (mean, weighted-average, etc)
                for stat in stats:
                    #compute the value/statistic
                    threshold = structList[tract]
                    if stat == 'wavg':
                        #i don't think it's necessary, but I rescale the mask to be from 0->1
                        val = compute_stat(data, normPath,stat=stat, rescale=True, threshold=threshold)
                    else:
                        #threshold = structList[tract]
                        if not threshold:
                            threshold = 0.005
                        val = compute_stat(data, normPath, threshold=threshold, stat=stat)
                    if not val:
                        #if the mask is empty...just skip it from here
                        print tract, 'has an empty or nan mask...'
                        continue
                    #set up a column name for hte mysql db
                    col = False ; col_vol = False ; col_vox = False
                    if not stat == 'vol':
                        col = tract + '_dti_' + scalarType + '_' + stat + '_' + scalar + var_sfix
                    else:
                        col_vol = tract + '_dti_vol' + var_sfix
                        col_vox = tract + '_dti_vox' + var_sfix
                    #print col, ':', val 
                    try:
                        #see if that column exists...if not...make it.
                        if col:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col))
                        if col_vol:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vol))
                        if col_vox:
                            self.cursor.execute("""alter table %s add column (%s float)""" % (tbl, col_vox))        
                    except:
                        pass
                    #update the data base with the data value...
                    print col
                    print col_vol
                    print col_vox
                    print val
                    if col:
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col, val, id_col, idval))
                    if col_vol:
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vol, val[1], id_col, idval))
                    if col_vox:
                        self.cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, col_vox, val[0], id_col, idval))
    

                    





