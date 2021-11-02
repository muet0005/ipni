import os as os
import sys as sys
import csv as csv
import nibabel as nb
import numpy as np
import MySQLdb as mysqldb



subjDir = '/Volumes/rbraid/mr_data_idc/aug2013_final/rsfmri'
drDir = os.path.join(subjDir, 'dual_regression/output_knicr_25Nov2015_fslFix_metaComponents/stage2')
melodicIC = os.path.join(subjDir, 'dual_regression/output_knicr_25Nov2015_fslFix_metaComponents/components/melodic_IC.nii.gz')

#numbering starts at 0!
components = [0, 1, 2, 3, 4, 5, 7, 11, 13, 14, 15, 16, 18, 19, 21]

pFix = 'idc_'
sFix = '_25Aug2014.ica'

subjs = os.listdir(subjDir)
#subjs = ['1']

#p-values and corresponding z-value
#0.05 = 1.96
#0.01 = 2.58
#0.001 = 3.29
zThresh = 3.29


#set some mysql db stuffs
print 'Negotiating connection to mysql DB...'
##connect to the mysql server
mydb = 'rsfmri'
con = mysqldb.connect(db = mydb, read_default_file=os.path.join(os.getenv("HOME"), '.my.cnf'))
##set up a cursor to use
cursor = con.cursor()
tbl = 'dual_regression_IC_stats_z329_21jan2016'
idCol = 'ID'

#make the table if it doesn't exist
def mkMysqlTbl(cursor, tbl, idCol):
    try:
        cursor.execute("""create table %s (%s varchar(10))""" % (tbl, idCol))
        return True
    except (mysqldb._mysql.OperationalError):
        return False

def insertSubMysql(cursor, idCol, tbl, idVal):
        #check the table to see if the subject exists
        cursor.execute("""select %s from %s where %s=\'%s\'""" % (idCol, tbl, idCol, idVal))
        pn_exist = cursor.fetchone()
        if not pn_exist:
            cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl, idCol, idVal))
        return pn_exist

def mkColMysql(cursor, tbl, colName):
    try:
        #see if that column exists...if not...make it.
        cursor.execute("""alter table %s add column (%s float)""" % (tbl, colName))
    except:
        pass

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
    """
    #set up some defaults
    stat = 'mean'
    rescale = False
    maskMinThresh = 1
    maskMaxThresh = False
    threshold = False
    val = False
    mask4d = False
    #check for optional args
    for i in kwargs.keys():
        if i == 'stat':
            stat = kwargs[i]
        elif i == 'rescale':
            rescale = kwargs[i]
        elif i == 'maskMinThresh':
            if kwargs[i] == 'min':
                maskMinThresh = False
            else:
                maskMinThresh = kwargs[i]
        elif i == 'maskMaxThresh':
            maskMaxThresh = kwargs[i]
        elif i == 'mask4d':
                mask4d = kwargs[i]
    #load in the mask data. First check if it is a 4d file
    is4d = len(nb.load(mask).shape)
    if is4d == 3:
        rawMaskData = nb.load(mask).get_data()
    else:
        rawMaskData = nb.load(mask).get_data()[:,:,:,mask4d]        
    #setup a new mask based on the thresholds specified above
    if not maskMaxThresh:
        maskMaxThresh = np.max(rawMaskData)
    if not maskMinThresh:
        maskMinThresh = np.min(rawMaskData)
    #print maskMinThresh
    #print maskMaxThresh
    maskData = np.zeros(np.shape(rawMaskData))
    maskData[ (rawMaskData>=maskMinThresh) & (rawMaskData<=maskMaxThresh)] = 1
    #get rid of any zeros in the data/stat image
    maskData[iData==0] = 0
    if stat == 'wavg':
        try:
            val = np.ma.average(np.ma.masked_where(maskData!=1, iData), weights=np.ma.masked_where(maskData!=1, rawMaskData))
        except ZeroDivisionError:
            val = False
    elif stat == 'mean':
        val = np.ma.mean(np.ma.masked_where(maskData!=1, iData))
    elif stat == 'median':
        val = np.ma.median(np.ma.masked_where(maskData!=1, iData))
    elif stat == 'stdev':
        val = np.ma.std(np.ma.masked_where(maskData!=1, iData))
    elif stat == 'vol':
        nvox = np.count_nonzero(maskData==1)
        mask_hd = nb.load(mask).get_header()
        x = float(mask_hd['pixdim'][1]) ;y = float(mask_hd['pixdim'][2]) ; z = float(mask_hd['pixdim'][3])
        vol = nvox * (x * y * z)
        val = [nvox, vol]
    if  not stat == 'vol':
        if np.isnan(val):
            print 'Warning....nan detected'
            val = False
    return val


def zpInt(x):
    if x < 10:
        xZp = '000' + str(c)
    elif x < 100:
        xZp = '00' + str(c)
    elif x < 1000:
        xZp = '0' + str(c)
    else:
        xZp = str(c)
    return xZp
    




mkMysqlTbl(cursor, tbl, idCol)

for c in components:
    cZp = zpInt(c)
    colName = 'dr_ica_' + cZp + '_positive_mean' 
    mkColMysql(cursor, tbl, colName)
    colName = 'dr_ica_' + cZp + '_negative_mean' 
    mkColMysql(cursor, tbl, colName)
    print 'Extracting data from component: ', cZp
    for s in subjs:
        sDir = os.path.join(subjDir, s, pFix + s + sFix)
        if os.path.exists(sDir):
            drMap = os.path.join(drDir, 'dr_stage2_idc_' + s +  '_ic' + cZp + '.nii.gz')
            if os.path.exists(drMap):
                insertSubMysql(cursor, idCol, tbl, s)
                try:
                    drData = False
                    drData = nb.load(drMap).get_data()
                except:
                    print 'FAILED TO OPEN: '
                    print drMap
                    continue
                colName = 'dr_ica_' + cZp + '_positive_mean' 
                val = compute_stat(drData, melodicIC, stat='mean', maskMinThresh=zThresh, mask4d=c)
                cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, colName, val, idCol, s))
                colName = 'dr_ica_' + cZp + '_negative_mean'
                val = compute_stat(drData, melodicIC, stat='mean', maskMinThresh='min', maskMaxThresh=zThresh*-1, mask4d=c)
                cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl, colName, val, idCol, s))



cursor.close()
con.close()
