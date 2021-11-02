
#import the stuff you need
import knicr.structural.freesurfer_utils as kfu
import os as os

#point to this FSGD file
fsgdFile = '/Volumes/rbraid/mr_data_idc/aug2013_final/freesurfer/desi/MTX/new/glm_fitpython/crude_sleep720.txt'
#list the subjects_dir
subjectsDir = '/Volumes/rbraid/mr_data_idc/aug2013_final/freesurfer'
#set an output folder to store everythong (not analysis specific)
oDir = '/Volumes/rbraid/mr_data_idc/aug2013_final/freesurfer/desi/MTX/new/glm_fitpython'
#specify where you want the GLM results to be saved (analysis specific)
outBaseName = 'left_hemi'
glmDir = os.path.join(oDir, outBaseName)


#specify the contrast file names. THEY HAVE TO LIVE IN 'oDir'
#I am actually mostly interested in C1 and C4 contrasts
contrasts = ['C1.mtx','C2.mtx']


#if you want o run sim-z, specify the parameters
#it is type, threshold 4=0.0001, sign (abs, neg, pos), and cluster-wise p (usually 0.05)
sim = ['mc-z', 1.3, 'abs', 0.05]


#load the module
#the default is to look at thickness. you also can specify other measures using measure = 'pial.area'
#or measure = 'volume' or measure = 'pial_lgi'
#for the hemi, the default is lh. to specify hemi, use hemi = 'rh'
#use below "overwrite=True" if you want to overwrite a previous analysis
#e.g., kfu.surfaceAnalysis(subjectsDir=subjectsDir, oDir=oDir, outBaseName=outBaseName, fsgdFile=fsgdFile, overwrite=True)
sa = kfu.surfaceAnalysis(subjectsDir=subjectsDir, oDir=oDir, outBaseName=outBaseName, fsgdFile=fsgdFile,overwrite=True)
#fix the excel-based text file to be linux friendly
sa.fix_newline()
#make a design mat in Matlab format from the fsgd/design hybrid
sa.fsgd_2_mat()
#make a 4D file of everyone
sa.mrisPreproc()
#set up the full paths to contrasts
contrastPaths = []
for c in contrasts:
    contrast = os.path.join(oDir, c)
    contrastPaths.append(contrast)
#run the glm
sa.mriGlmFit(glmdir = glmDir, contrasts = contrastPaths)
#correct for multiple testing
sa.mriGlmFitSim(sim=sim, oneHemi=True)
