import os as os
import sys as sys
import csv as csv
import dicom as pd
import tarfile as tarfile
import codecs

__author__ = "Ryan Muetzel"
__license__ = "GPL"
__version__ = "0.1"

class cvtRtoIDC():
    def __init__(self, **kwargs):
        """
        cvtRtoIDC is a class of functions to convert GenR MRI data to IDC
        
        Optional arguments:
        idcDir = where the idc tar file lives.
        baseName = the name of the tar file and the csv file living inside of the tar archive..
        
        Usage:
        
        >>>import knicr.datamanagement.r_to_idc_f9 as kdri
        #initalize the IDC converter tool
        >>>cvtIdc = kdri.cvtRtoIDC()
        #read in the idc data
        >>>cvtIdc.readIdcLink()
        #get the IDC of a case
        >>>IDC, dcmUUID = cvtIdc.getIdc(mrsessionLabel, dcm)
        
        """
        self.idcDir = '/Users/rmuetzel/.focus_at_9/.idc/'
        self.baseName = '.Ryan_MRIConversion_March2013-November2015_Part1_rlm'
        for k in kwargs.keys():
            if k == 'idcDir':
                self.idcDir = kwargs[k]
            elif k == 'baseName':
                self.baseName = kwargs[k]
        self.idcLink = os.path.join(self.idcDir, self.baseName + '.tar')
        if not os.path.exists(self.idcLink):
            print 'Cannot Find IDC link file: ', self.idcLink
            sys.exit(9)
    
    def readIdcLink(self):
        """
        readIdcLink will read a tar file, and extract the data from the included csv.
        
        """
        #First, open the tar file into memory
        tar = tarfile.open(self.idcLink)
        f = False
        #find the csv file inside the archive (should only be 1)
        for m in tar.getmembers():
            #make sure it has the right name!
            if m.path == self.baseName + '.csv':
                #pull out that one file
                f = tar.extractfile(m)
            if not f:
                print 'Cannot find the idc link file in the tar file...'
                print 'tar file: ', self.idcLink, 
                print 'csv file: ', self.baseName + '.csv'
                print 'tar member file', m
                return
            #read the data into memory, let it know it is csv data
            content = csv.reader(f.readlines(), delimiter=',')
            #dump the data into a list 
            self.idcLinkData = []
            for i in content:
                #first fix the headers, in case spss laid a little utf8 terd. I couldn't get the codecs decoder to work here
                if i[0].startswith('\xef\xbb\xbf'):
                    i[0] = i[0].strip('\xef\xbb\xbf')
                self.idcLinkData.append(i)
            #assume the first line is the header information
            #store that info into a dictionary for use later
            self.headerDict = {}
            for i in self.idcLinkData[0]:
                self.headerDict[i] = self.idcLinkData[0].index(i)
            #now set up a dictonary based on the session UUID of the dicom
            if not self.headerDict.has_key('rnummer'):
                return False
            self.uuidDict = {}
            for i in self.idcLinkData:
                if i[self.headerDict['rnummer']] == 'rnummer':
                    continue
                self.uuidDict[i[self.headerDict['ScanUUID_DICOM']]] = i
            return True
    
    def getIdc(self, xnat_session, dcm):
        """
        getIdc will get an idc number, based on the dicom uuid, from the info read in by readIdcLink.
        
        """
        #make sure the dicom file exists
        if not os.path.exists(dcm):
            print 'Cannot find dicom: ', dcm
            return False
        #get the UUID in that dicom
        dcmUUID = pd.read_file(dcm).StudyInstanceUID
        IDC = False
        #see if that UUID is ready to be converted
        if self.uuidDict.has_key(dcmUUID):
            #if so, make sure the xnat session you think it is, is the real xnat session
            if self.uuidDict[dcmUUID][self.headerDict['xnat_session']] == xnat_session:
                #if so, get the IDC associated with that UUID
                IDC =  self.uuidDict[dcmUUID][self.headerDict['IDC']]
        return IDC, dcmUUID, os.path.join(self.idcDir, '.logs')
    


