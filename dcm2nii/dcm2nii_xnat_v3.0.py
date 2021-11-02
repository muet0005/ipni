import os as os
import sys as sys
import csv as csv
import shutil as shutil
import pyxnat
import socket as socket
import datetime as datetime
from pyxnat import Interface
from nipype.interfaces import dcm2nii
from nipype.interfaces import afni
import dicom as pd
import knicr.datamanagement.r_to_idc_f9 as kdri
import multiprocessing as mp


#set your project. Can only be one at a time.
project = 'GenR_F9_MRI'
#where are the dicom data?
dcmDir = os.path.join('/Volumes/rbraid/xnat/xnat_rb/archive', project, 'arc001')
#where you want the nifti images to be saved out.
niiDir = os.path.join('/Volumes/rbraid/xnat/nii', project)

irDir = '/Volumes/rbraid/xnat/nii/software/isRunning'

convertIdc = True
cpDicom = False
parallel = False
ncores = 5

__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "1.1"


#these data have already been converted, but have goofy names that were manually changed...
#this will ensure they won't get re-converted again with the same goofy name
subs_exclude=['CogMed3']

#the project that contains your data

#this makes it easier to initialize a connection with xnat...
#put the .pyxnat.rv.cfg file in your home folder, and chmod 700
pyxnat_cfg = os.path.join(os.environ['HOME'], '.pyxnat.rb.cfg')

#a tuple of seq names, in case there are more than 1 variant
seq_names = ('T1_weighted_FreeSurfer_BRAVO_IR_ARC2_1mm_Cor', 'DTI-Ax-35DIR-2.5mm_TR11.5_TE93_b900_PA', 'RestingState-200Volumes-36slices-TR1.76sec-Asset2','3D CUBE Eyes Loc', 'PU:T1_weighted_FreeSurfer_BRAVO_IR_ARC2_1mm_Cor', 'SIRP-verb-WkMem-TR1.76sec-Asset2', 'SIRP_visspac-WkMem-TR1.76sec-Asset2', 'T1_weighted_FreeSurfer_BRAVO_IR-SPGR_ASSET2_1mm_Cor')

seq_dict = {'T1_weighted_FreeSurfer_BRAVO_IR_ARC2_1mm_Cor':'t1', 'DTI-Ax-35DIR-2.5mm_TR11.5_TE93_b900_PA':'dti_pa', 'RestingState-200Volumes-36slices-TR1.76sec-Asset2':'rsfmri','3D CUBE Eyes Loc':'t2', 'PU:T1_weighted_FreeSurfer_BRAVO_IR_ARC2_1mm_Cor':'t1_pure', 'SIRP-verb-WkMem-TR1.76sec-Asset2': 'sirp_verb_wm', 'SIRP_visspac-WkMem-TR1.76sec-Asset2': 'sirp_vispat_wm', 'T1_weighted_FreeSurfer_BRAVO_IR-SPGR_ASSET2_1mm_Cor':'t1_asset'}

afniBrik = ('RestingState-200Volumes-36slices-TR1.76sec-Asset2')
#set to true to overwrite existing nifti images...
overwriteExisting = False

#connect to the db
central = Interface(config=pyxnat_cfg)

#set up some functions to re-use later
def get_scanData(session_id):
    #get the scan data you're interested in.  In this case, it's the series description (scan name) and the series ID number (the number on the disk)
    constraint = [('xnat:mrScanData/IMAGE_SESSION_ID', '=', session_id)]
    f = central.select('xnat:mrScanData', ['xnat:mrScanData/SERIES_DESCRIPTION', 'xnat:mrScanData/ID', 'xnat:mrScanData/QUALITY']).where(constraint)
    return f.items()


def eval_scans(scans):
    #function to determine whether or not a scan should be selected
    #first see if the scan type is an eligible type (scan_dict)
    #If there are more than 1 instance of a given scan type, it then looks at quality
    usable_scan_dict = {}
    mult_scan_dict = {}
    for series_desc, series_id, quality in scans:
	try:
		test=int(series_id)
	except:
		continue
        if series_desc.startswith(seq_names):
            if not usable_scan_dict.has_key(series_desc):
                usable_scan_dict[series_desc] = [series_id, quality, series_desc]
            else:
                ###add scan specific things here for the DTI/FMRI/T1?
                if ratings[quality] >= ratings[usable_scan_dict[series_desc][1]]:
                    if int(series_id) > int(usable_scan_dict[series_desc][0]):
                        usable_scan_dict[series_desc] = [series_id, quality, series_desc]
    for series_desc, series_id, quality in scans:
        if series_desc.startswith(seq_names):
            if not mult_scan_dict.has_key(series_desc):
                mult_scan_dict[series_desc] = [series_id]
            else:
                exSeries = mult_scan_dict[series_desc]
                exSeries.append(series_id)
                mult_scan_dict[series_desc] = exSeries
    return usable_scan_dict, mult_scan_dict


def write_log(cmd, oFile, addcmd, seriesInfo):
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
    f.write('SeriesInfo=' + str(seriesInfo) + '\n')
    if addcmd:
        f.write(cmd + '\n')
    else:
        f.write('COMMAND SUPRESSED! CHECK IDC CONVERSION LOG' + '\n')
    f.write('#######--------END OF LOG--------#######' + '\n')
    f.close()


def writeidcLog(idc, mrsessionLabel, uuid, idcLogDir):
    ts = datetime.datetime.now().strftime('%d%b%Y_%H:%M:%S')
    f = open(os.path.join(idcLogDir, '.' + idc + '_idcCvtLog_' + ts), 'w')
    f.write(idc + ',' + mrsessionLabel + ',' + uuid)


def getDcm4dcm2nii(scanDir):
    dcms = os.listdir(scanDir)
    i = 0
    dcm = False
    while not dcm and i < 50:
        if dcms[i].endswith('.dcm'):
            dcm = os.path.join(scanDir, dcms[i])
        else:
            i+=1
    return dcm


def genrDcm2nii(hq, dcm, oniiDir, series_desc, series_id, quality):
    if not os.path.exists(oniiDir):
        os.makedirs(oniiDir)
    #set the output name based on the seq_dictionary
    if hq:
        niiOutName = seq_dict[series_desc] + '_' + str(series_id) + '_best'
    else:
        niiOutName = seq_dict[series_desc] + '_' + str(series_id)
    if os.path.exists(os.path.join(oniiDir, niiOutName + '.nii.gz')) and not overwriteExisting:
        return
    #run the dcm2nii conversion
    dcmCvt = dcm2nii.Dcm2nii(source_names=dcm, terminal_output='stream', gzip_output=True, output_dir=oniiDir, reorient=False, reorient_and_crop=False)
    #print dcmCvt.cmdline
    dcmCvt.run()
    #it outputs some nonsense, so we have to rename it
    os.rename(dcmCvt.output_files[0], os.path.join(oniiDir, niiOutName + '.nii.gz'))
    if hasattr(dcmCvt, 'bvecs'):
        if len(dcmCvt.bvecs) > 0:
            os.rename(dcmCvt.bvecs[0], os.path.join(oniiDir, niiOutName + '.bvec'))
            os.rename(dcmCvt.bvals[0], os.path.join(oniiDir, niiOutName + '.bval'))
    #data prov.
    addcmd = True
    if convertIdc:
        #remove any instance of R-number...
        addcmd = False
    logFile = os.path.join(oniiDir, niiOutName + '.log')
    write_log(dcmCvt.cmdline, logFile, addcmd, [series_id, series_desc, quality])

def genrDcm2afni(hq, dcm, oniiDir, series_desc, series_id, quality, scanType):
    if not os.path.exists(oniiDir):
        os.makedirs(oniiDir)
    #set the output name based on the seq_dictionary
    if hq:
        afniOutName = seq_dict[series_desc] + '_' + str(series_id) + '_best'
    else:
        afniOutName = seq_dict[series_desc] + '_' + str(series_id)
    if os.path.exists(os.path.join(oniiDir, afniOutName)) and not overwriteExisting:
        return
    #gather some fun facts about the scan
    dcmInfo = pd.read_file(dcm, stop_before_pixels=True)
    repetitionTime = int(dcmInfo.RepetitionTime)
    numRepetitions = int(dcmInfo.NumberOfTemporalPositions)
    numImages = int(dcmInfo.ImagesInAcquisition)
    #this tag is not easy to get, so you need to first get a pydicom friendly tag
    nSlices = int(dcmInfo[pd.tag.BaseTag('2166863L')].value)
    if not nSlices == (numImages/numRepetitions):
        print('something is not right...num slices not matching num dicom images')
    #run the to3d conversion
    optStr = '-session ' + oniiDir + ' -time:zt ' + str(nSlices) + ' ' + str(numRepetitions) + ' ' + str(repetitionTime)
    to3d = afni.To3D(in_folder=os.path.dirname(dcm), out_file=afniOutName, outputtype='AFNI', args=optStr, terminal_output='stream', ignore_exception=True)
    print(to3d.cmdline)
    to3d.run()
    addcmd = True
    if convertIdc:
        #remove any instance of R-number...
        addcmd = False
    logFile = os.path.join(oniiDir, afniOutName + '.log')
    write_log(to3d.cmdline, logFile, addcmd, [series_id, series_desc, quality, scanType])    
    
    
ratings = {'Excellent': 6, 'Good':5, 'Questionable': 4, 'Poor':3, 'Unusable/Other(comment)':2, 'usable':1, 'Unrated':1}


#subjs = get_subjs(project)
project_info = central.select.project(project)
######mrsessions = project_info.experiments().get()
mrsessions = ['GenerationR_E03255']
#mrsessions = ['GenerationR_E01815']
#mrsessions = ['GenerationR_E01419']
mrsessions = ['GenerationR_E01665']

#initalize the IDC converter tool
cvtIdc = kdri.cvtRtoIDC()
#read in the idc data
cvtIdc.readIdcLink()

noIdc = []


def convertMrSession(mrsession):
    lockFile = os.path.join(irDir, mrsession + '.Dcm2Nii.IsRunning')
    if os.path.exists(lockFile):
        return
    lf = open(lockFile, 'w')
    lf.close()
    idcLogDir = False
    #get the xnat mr_session label (what is displayed in xnat)
    mrsessionLabel = central.select.experiment(mrsession).label()
    if mrsessionLabel in subs_exclude:
        return
    #for each mr session, get some more info
    print "SESSION: ", mrsessionLabel
    #get all of the mri scans/sequences within a given session
    scans = get_scanData(mrsession)
    #get additional details for each session
    usable_scan_dict, mult_scan_dict = eval_scans(scans)
    #make a copy of the dicom? Only if the IDC dingetje is turned on
    if cpDicom and convertIdc:
        cpDcmDir = os.path.join(niiDir, subj, 'dicom')
        os.makedirs(cpDcmDir)
        #    shutil.copytree(os.path.join(dcmDir, subj, 'SCANS', series_id), os.path.join(cpdcmDir, series_id))
        #add some stuffs in here to anonymize the dicom
    uuid = False
    IDC = False
    for series_desc in mult_scan_dict.keys():
        for series_id in mult_scan_dict[series_desc]:
	    try:
                test = int(series_id)
            except:
		continue
            hq = False
            quality = False
            if series_id == usable_scan_dict[series_desc][0]:
                quality = usable_scan_dict[series_desc][1]
                hq = True
            scanDir = os.path.join(dcmDir, mrsessionLabel, 'SCANS', str(series_id), 'DICOM')
            if not os.path.exists(scanDir):
                print 'PATH DOES NOT EXIST: ', scanDir
                print "Why not?????"
                continue
            #find a dicom to point to dcm2nii
            dcm = getDcm4dcm2nii(scanDir)
            if convertIdc and project == 'GenR_F9_MRI':
                IDC, dcmUUID, idcLog = cvtIdc.getIdc(mrsessionLabel, dcm)
                if not idcLogDir:
                    idcLogDir = idcLog
                if not IDC:
                    print 'Cannot yet convert session: ', mrsessionLabel
		    noIdc.append(mrsessionLabel)
                    return
                if not uuid:
                    uuid = dcmUUID
                else:
                    if not uuid == dcmUUID:
                        print 'PROBLEM!!!! Found more than 1 UUID for a given scan session!!!', mrsessionLabel
                        return
                oniiDir = os.path.join(niiDir, IDC)
                #print oniiDir
            else:
                oniiDir = os.path.join(niiDir, mrsessionLabel)
            genrDcm2nii(hq, dcm, oniiDir, series_desc, series_id, quality)
            if series_desc.startswith(afniBrik):
                ##THIS MUST BE CHANGED IF YOU WANT TO USE ANYTHING OTHER THAN RSFMRI
                scanType = 'epan'
                genrDcm2afni(hq, dcm, oniiDir, series_desc, series_id, quality, scanType)
    if convertIdc and IDC:
        #save a record of the conversion in a secure location
        writeidcLog(IDC, mrsessionLabel, uuid, idcLogDir)


if not __name__ == '__main__':
    parallel = False

if parallel:
    pool = mp.Pool(processes=ncores)
    pool.map(convertMrSession, mrsessions)
else:
    print mrsessions
    for mrsession in mrsessions:
        convertMrSession(mrsession)


print 'no IDC for: ', len(noIdc)
print 'number of IDCs:  ', len(cvtIdc.uuidDict.keys())
print 'idc link: ', cvtIdc.idcLink
