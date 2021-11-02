import os as os
import nipype.interfaces.fsl as fsl

#specify the inputs
DIR = '/Volumes/rbraid/mr_data_idc/aug2013_final/rsfmri/melodic_samples_d25_n50_s100C/matched'
oDIR = os.path.join(DIR, 'merged')

#templates = ['cerebellum', 'DMN', 'inferior_mid_frontal', 'insula_subcortical', 'left_pf', 'mid_frontal', 'noise_ant_frontal', 'noise_lower_brainstem', 'noise_pons_vessel', 'noise_sinus', 'noise_sup_frontal', 'noise_susceptibility', 'noise_upper_brainstem', 'noise_vent_wm', 'parietal', 'right_pf', 'sensory_motor', 'superior_mid_frontal', 'visual']
#templates = ['ant_frontal', 'auditory', 'cerebellum', 'dmn', 'executive_control', 'inf_mid_frontal', 'inf_parietal_temporal', 'insula_acc', 'lat_dorsal_frontal', 'left_pf', 'mid_frontal', 'parietal_occipital', 'pcc', 'pcc_parietal', 'posterior_cerebellum', 'right_pf', 'sensorimotor', 'subcortical', 'visual', 'noise_inf_brainstem', 'noise_lower_brainstem_cerebellum', 'noise_pons_large', 'noise_pons_small', 'noise_sinus', 'noise_sup_frontal', 'noise_susceptibility', 'noise_temporal_vessel', 'noise_wm_vent', 'noise_cerebellum_occipital']
templates = ['sup_insula_acc', 'lat_dor_frontal', 'ant_frontal', 'auditory','cerebellumA', 'cerebellumB', 'executive_control', 'inf_mid_frontal', 'insula_acc', 'lang_atten', 'left_pf', 'mid_frontal', 'noise_inf_brainstem', 'noise_lower_brainstem_cerebellum', 'noise_occ_cerebellum', 'noise_pons_large', 'noise_pons_small', 'noise_sinus', 'noise_subcortical', 'noise_sup_frontal', 'noise_temporal_vessel', 'noise_wm_vent', 'parietal_occipital', 'pcc', 'pcc_parietal', 'right_pf', 'sensory', 'sup_motor', 'visual']

nsamples = 100

for template in templates:
	oFile = os.path.join(DIR, 'merged', template + '_merged.nii.gz')
	if os.path.exists(oFile):
		continue
	print template
	fList = []
	for s in range(0, nsamples):
		if s < 10:
			s_str = '000' + str(s)
		elif s < 100:
			s_str = '00' + str(s)
		elif s < 1000:
			s_str = '0' + str(s)
		else:
			s_str = str(s)
		f = os.path.join(DIR, template, template + '_' + s_str + '.nii.gz')
		if os.path.exists(f):
			fList.append(f)
	fslmerge = fsl.Merge(dimension='t', terminal_output='stream',in_files=fList, merged_file=oFile, output_type='NIFTI_GZ')
	fslmerge.run()
for template in templates:
	print template
	iFile = os.path.join(DIR, 'merged', template + '_merged.nii.gz')
	oFile = os.path.join(DIR, 'mean', template + '_merged.nii.gz')
	fslmaths = fsl.ImageMaths(in_file=iFile, op_string= '-Tmean', out_file=oFile, output_type='NIFTI_GZ', out_data_type='float')
	fslmaths.run()