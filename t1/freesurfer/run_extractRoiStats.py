import os as os
import sys as sys
import knicr.structural.freesurfer_utils as kfu

subjects_dir = '/Volumes/rbraid/mr_data_idc/aug2013_final/freesurfer'
subjects = os.listdir(subjects_dir)
mydb = 'freesurfer_v51'

extractStats = kfu.extractFreeSurferStats(subjects_dir, mydb=mydb, varSfix='_f5')

swExclude = ('kn', 'fsaverage', 'lh', 'rh')
for subject in subjects:
    statDir = os.path.join(subjects_dir, subject, 'stats')
    if not os.path.exists(statDir) or subject.startswith(swExclude):
        continue
    print subject
    extractStats.dumpAparcStats(subject)
    extractStats.dumpBAStats(subject)
    extractStats.dumpAsegStats(subject)
    extractStats.dumpWmparcStats(subject)
    extractStats.dumpTbvStats(subject)