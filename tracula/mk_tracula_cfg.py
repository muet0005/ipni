import os as os
import sys as sys

###Read in the passed arguments
subjects_dir = sys.argv[1]
dtroot = sys.argv[2]
subjid = sys.argv[3]
cfg_file = sys.argv[4]

#open up the empty configuration file
cfg = open(cfg_file, 'w')

#Write in the necessary lines/options
cfg_str = 'setenv SUBJECTS_DIR ' + subjects_dir
cfg.write(cfg_str + '\n')
cfg_str = 'set dtroot = ' + dtroot
cfg.write(cfg_str + '\n')
cfg_str = 'set subjlist = (' + subjid + ')'
cfg.write(cfg_str + '\n')
cfg_str = 'set runlist = (1)'
cfg.write(cfg_str + '\n')
cfg_str = 'set dcmroot = ' + os.path.join(dtroot, subjid)
cfg.write(cfg_str + '\n')
cfg_str = 'set dcmlist = (dti_idc_' + subjid + '_2mm.nii.gz)'
cfg.write(cfg_str + '\n')
cfg_str = 'set bvalfile = ' + os.path.join(dtroot, subjid + '/dti_idc_' + subjid + '.transposed.bval')
cfg.write(cfg_str + '\n')
cfg_str = 'set bvecfile = ' + os.path.join(dtroot, subjid + '/dti_idc_' + subjid + '.transposed.bvec')
cfg.write(cfg_str + '\n')
cfg_str = 'set nb0 = 3'
cfg.write(cfg_str + '\n')
cfg_str = 'set doeddy = 1'
cfg.write(cfg_str + '\n')
cfg_str = 'set dorotbvecs = 1'
cfg.write(cfg_str + '\n')
cfg_str = 'set thrbet = 0.25'
cfg.write(cfg_str + '\n')
cfg_str = 'set doregflt = 0'
cfg.write(cfg_str + '\n')
cfg_str = 'set doregbbr = 1'
cfg.write(cfg_str + '\n')
cfg_str = 'set doregmni = 1'
cfg.write(cfg_str + '\n')
cfg_str = 'set mnitemp = $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz'
cfg.write(cfg_str + '\n')
#in v5.3, some of these options seem to have changed. The defaults are now included as "all", so these are now comented out.
#cfg_str = 'set trainfile = /home/genr/software/bitbucket/lisa/tracula/train_lists/trainlist.txt'
#cfg.write(cfg_str + '\n')
#cfg_str = 'set pathlist = (lh.cst_AS rh.cst_AS lh.ilf_AS rh.ilf_AS lh.unc_AS rh.unc_AS fmajor_PP fminor_PP lh.atr_PP rh.atr_PP lh.cab_PP rh.cab_PP lh.ccg_PP rh.ccg_PP lh.slfp_PP rh.slfp_PP lh.slft_PP rh.slft_PP)'
#cfg.write(cfg_str + '\n')
#cfg_str = 'set ncpts = 5'
#cfg.write(cfg_str + '\n')
cfg.close()
