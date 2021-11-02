import os as os
import knicr.fmri.fmri_utils as kfmri

#set the folder
DIR = '/home/genr/data/rsfmri'

#set n of subjects per sample
nSubs = 50
#set n of samples
nSamples = 1000
#if you already run x samples, you can start at that point again, for ease of combining data later
startSample = 0
#set the folder name
dirname = 'melodic_samples_d16_n50_s1000'

mDIR = os.path.join(DIR, dirname)

subj_list = os.path.join(mDIR, 'RLM_template_sample_n494_noClinProb_moco3mm.txt')

gbs = kfmri.gen_bootstrap_samples()

gbs.read_subj_list(subj_list)

gbs.gen_samples(startSample, nSamples, nSubs)

gbs.gen_melodic_lisa(startSample, DIR, 'idc_', '_27July2013.feat', mDIR, dirname)

