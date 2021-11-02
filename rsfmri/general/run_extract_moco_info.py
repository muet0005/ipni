import os as os
import csv as csv
import knicr.fmri.extract_mcflirt_moco_par as moco



DIR = '/Volumes/rbraid/mr_data_idc/aug2013_final/rsfmri'
idcs = os.listdir(DIR)
feat_dir_pfix = 'idc_'
feat_dir_sfix = '_27July2013.feat'

moco_summary = open(os.path.join(DIR, 'MCFLIRT_moco_summary_19Nov2013.csv'), 'w')

max_translation = 3

headers = ['idc', 'exclude_max_trans_'+str(max_translation), 'max_abs_trans', 'mean_rel_trans']
data = []
data.append(headers)

for idc in idcs:
    feat_dir = os.path.join(DIR, idc, feat_dir_pfix + idc + feat_dir_sfix)
    if os.path.exists(feat_dir):
        x = moco.eval_mcflirt_moco(feat_dir)
        x.read_mcflirt_moco()
        exclude_case,max_abs,mean_rel = x.extract_moco_summary(maxtrans=max_translation)
        print idc, exclude_case, max_abs, mean_rel
        t = [idc, exclude_case, max_abs, mean_rel]
        data.append(t)


cw = csv.writer(moco_summary, delimiter=',')
cw.writerows(data)
moco_summary.close()