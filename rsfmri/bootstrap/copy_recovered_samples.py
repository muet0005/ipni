import os as os
import sys as sys
import shutil as shutil

DIR = '/Volumes/rbraid-1/mr_data_idc/aug2013_final/rsfmri/melodic_samples_d25_n50_s1000/matched.recovered'

umDir = os.path.join(DIR, 'unmatched')

unmatched = os.listdir(umDir)
recovered = []

for i in unmatched:
    if i.startswith('recovered'):
        recovered.append(i)

print 'found: ', len(recovered)

unique_cmpts = []
for r in recovered:
    print r
    cmpt = r.split('_UNMATCHED')[0].split('recovered_')[1]
    sample = r.split('sample')[1].split('.')[0]
    rec_cmpt = os.path.join(DIR, cmpt, cmpt + '_' + sample + '.nii.gz')
    if not os.path.exists(rec_cmpt):
        print cmpt, sample, rec_cmpt
        shutil.copyfile(os.path.join(DIR, 'unmatched', r), rec_cmpt)
    if cmpt not in unique_cmpts:
        unique_cmpts.append(cmpt)
    else:
        print 'Something is wrong...this component should not be recovered....'
        print cmpt, sample, rec_cmpt
        #sys.exit(0)
        

for i in unique_cmpts:
    print i
    