import os as os
import nibabel as nb
import numpy as np
#from scipy import spatial as sp
#from scipy.spatial.distance import rogerstanimoto
#from scipy.stats.stats import pearsonr as pr
import scipy.stats.stats as spstat
import nipype.interfaces.fsl as fsl
from multiprocessing import Pool
import csv as csv

sample_sfix = 'melodic_samples_d25_n50_s100C'
#specify the inputs
DIR = os.path.join('/Volumes/rbraid/mr_data_idc/aug2013_final/rsfmri', sample_sfix)
oDIR = os.path.join(DIR, 'matched')
templateDir = '/Volumes/rbraid/mr_data_idc/aug2013_final/rsfmri/melodic_ica_templates/d25'
#templates = ['cerebellum', 'DMN', 'inferior_mid_frontal', 'insula_subcortical', 'left_pf', 'mid_frontal', 'noise_ant_frontal', 'noise_lower_brainstem', 'noise_pons_vessel', 'noise_sinus', 'noise_sup_frontal', 'noise_susceptibility', 'noise_upper_brainstem', 'noise_vent_wm', 'parietal', 'right_pf', 'sensory_motor', 'superior_mid_frontal', 'visual']
#templates = ['ant_frontal', 'auditory', 'cerebellum', 'dmn', 'executive_control', 'inf_mid_frontal', 'inf_parietal_temporal', 'insula_acc', 'lat_dorsal_frontal', 'left_pf', 'mid_frontal', 'parietal_occipital', 'pcc', 'pcc_parietal', 'posterior_cerebellum', 'right_pf', 'sensorimotor', 'subcortical', 'visual', 'noise_inf_brainstem', 'noise_lower_brainstem_cerebellum', 'noise_pons_large', 'noise_pons_small', 'noise_sinus', 'noise_sup_frontal', 'noise_susceptibility', 'noise_temporal_vessel', 'noise_wm_vent', 'noise_cerebellum_occipital']
templates = ['noise_susceptibility', 'par_occ_temp', 'pcc_frontal', 'pcc_insula','sup_insula_acc', 'lat_dor_frontal', 'ant_frontal', 'auditory','cerebellumA', 'cerebellumB', 'executive_control', 'inf_mid_frontal', 'insula_acc', 'lang_atten', 'left_pf', 'mid_frontal', 'noise_inf_brainstem', 'noise_lower_brainstem_cerebellum', 'noise_occ_cerebellum', 'noise_pons_large', 'noise_pons_small', 'noise_sinus', 'noise_subcortical', 'noise_sup_frontal', 'noise_temporal_vessel', 'noise_wm_vent', 'parietal_occipital', 'pcc', 'pcc_parietal', 'right_pf', 'sensory', 'sup_motor', 'visual']
ndim=25
nsamples = 100
ncores = 25

#for thresholding the IC maps....they are Z maps
z = 3.3

#make a function -- two args are the trg and src maps
def compute_ui(trg_data, src_data):
    #read in the nifti data for the target (template components) and the source (the experimental set that needs matching)
    #trg_nii = 
    #src_nii = 
    #pull out the data array
    #trg_data = nb.load(trg_map).get_data()
    #src_data = nb.load(src_map).get_data()
    #make an empty array to store the union over intersection values
    #this is trg# of components x src# of components
    ui = np.zeros([trg_data.shape[3], src_data.shape[3]])
    #scc = np.zeros([trg_data.shape[3], src_data.shape[3]])
    #tdis = np.zeros([trg_data.shape[3], src_data.shape[3]])
    #loop through the trg components, and the src components
    for c_trg in range(0, trg_data.shape[3]):
        #print 'checking component: ', c_trg
        for c_src in range(0, src_data.shape[3]):
            #first threshold to get the positive values in both maps
            trg_mask_pos = trg_data[:,:,:,c_trg] >= z
            src_mask_pos = src_data[:,:,:,c_src] >= z
            #next identify the negative values
            #trg_mask_neg = trg_data[:,:,:,c_trg] <= -(z)
            #src_mask_neg = src_data[:,:,:,c_src] <= -(z)
            #find the intersection of the positive map
            intersect_pos = np.sum((trg_mask_pos.astype(int) + src_mask_pos.astype(int)) >=2)
            #find the union of the positive map
            union_pos = np.sum(trg_mask_pos.astype(int)) + np.sum(src_mask_pos.astype(int)) - intersect_pos
            #same for the negative maps
            #intersect_neg = np.sum((trg_mask_neg.astype(int) + src_mask_neg.astype(int)) >=2)
            #union_neg = np.sum(trg_mask_neg.astype(int)) + np.sum(src_mask_neg.astype(int)) - intersect_neg
            #now, sum the pos/neg union/intersects
            #intersect = intersect_pos + intersect_neg
            #union = union_pos + union_neg
            intersect = intersect_pos
            union = union_pos
            #dump the union over intersection into the array with the coordinates that correspond to which component matches
            ui[c_trg, c_src] = float(intersect)/float(union)
            #scc[c_trg, c_src] = pr(trg_data[:,:,:,c_trg].ravel(), src_data[:,:,:,c_src].ravel())[0]
            #tdis[c_trg, c_src] = rogerstanimoto(trg_mask_pos.ravel(), src_mask_pos.ravel())
    #make an empty matrix to store the matches
    #this is an array that is 2 x n-components
    matched_components = np.zeros([2, trg_data.shape[3]])
    #now iterate through each target component, and find the max U/I...this is your match
    for i in range(0, trg_data.shape[3]):
        #here we fill this array with the volume where the max UI is (1st row), and also the UI value (2nd row)
        matched_components[0][i] = np.argmax(ui[i])
        matched_components[1][i] = np.max(ui[i])
    return matched_components


#a function to read in the power spectra, and select the column for a given component
def read_ftmix(ftmix, component):
    f = open(ftmix, 'r')
    csv_ft = csv.reader(f, delimiter=' ')
    ft = []
    for line in csv_ft:
        rm_ws = True
        while rm_ws:
            try:
                line.remove('')
            except:
                rm_ws = False
        ft.append(line)
    ft_component = []
    for line in ft:
        ft_component.append(line[component])
    f.close()
    return ft_component


#a function to write out a file with the power spectrum for the column of interest inthe melodic_FTmix file
def write_ftmix(ftmix_data, ftmix_file):
    f = open(ftmix_file, 'w')
    for i in ftmix_data:
        f.write(i + '\n')
    f.close()



trg_maps = {}
print 'pre-loading templates...'
for template in templates:
    print 'loading...', template,
    trg_map = os.path.join(templateDir, template + '_merged.nii.gz')
    trg_maps[template] = nb.load(trg_map).get_data()
    print '....done.'


def run_sample(sample):
    if sample < 10:
        str_sample = '000' + str(sample)
    elif sample < 100:
        str_sample = '00' + str(sample)
    elif sample < 1000:
        str_sample = '0' + str(sample)  
    else:
        str_sample = str(sample)
    #make an empty dictionary to hold the matches
    matched_components = {}
    src_map = os.path.join(DIR, 'sample.' + str(sample) + '.' + sample_sfix, 'melodic_IC.nii.gz')
    src_data = nb.load(src_map).get_data()
    for template in templates:
        print template
        #trg_map = os.path.join(DIR, 'templates', template + '_merged.nii.gz')
        #for each template --- add the matched component info
        #this is a 2 x n-component numpy array
        #The first row holds the volume number from the source/input data that had the max UI with a given tempalte
        #the second column holds the UI for that volume
        matched_components[template] = compute_ui(trg_maps[template], src_data)
    #now we make a final dictionary to hold the actual final decisions on matches
    template_matches = {}
    for template in matched_components:
        print template
        #find the mode, most commonly matched volume for a given template
        dict_key = int(spstat.mode(matched_components[template][0])[0][0])
        print 'vols: ', matched_components[template][0]
        print 'max: ', matched_components[template][1]
        print 'mode: ', spstat.mode(matched_components[template][0])[0][0]
        if int(spstat.mode(matched_components[template][0])[1][0]) == 1:
            print 'WARNING.....THERE IS NO TRUE MODE for the following: '
            print 'template: ', template
            print 'max match is: ', np.max(matched_components[template][1])
        if not template_matches.has_key(dict_key):
            template_matches[dict_key] = template
            print dict_key, template
        else:
            max_new = np.max(matched_components[template][1])
            max_existing = np.max(matched_components[template_matches[spstat.mode(matched_components[template][0])[0][0]]][1])
            if max_new > max_existing:
                print 'Dictionary key: ', dict_key, ' was reassigned from ', template_matches[dict_key], ' to: ', template
                if not os.path.exists(os.path.join(oDIR, 'unmatched')):
                    os.makedirs(os.path.join(oDIR, 'unmatched'))
                f = open(os.path.join(oDIR, 'unmatched', template_matches[dict_key] + '_sample' + str(sample) + '_' + str(dict_key) + 'reassigned.txt'), 'a')
                f.write('Dictionary key: ' +  str(dict_key) + ' was reassigned from ' +  template_matches[dict_key] +  ' to: ' +  template + '\n')
                f.write('max new:' + str(max_new) + 'max exisiting: ' + str(max_existing) + '\n')
                f.close()
                template_matches[dict_key] = template
                print max_new, max_existing
    vol_selected = []
    for vol in template_matches.keys():
        print 'Writing out components...'
        vol_selected.append(vol)
        component = template_matches[vol]
        if not os.path.exists(os.path.join(oDIR, component)):
            os.makedirs(os.path.join(oDIR, component))
        oFile = os.path.join(oDIR, component, component + '_' + str_sample + '.nii.gz')
        fslroi = fsl.ExtractROI(in_file=src_map, roi_file=oFile, t_min=vol, t_size=1)
        fslroi.run()
        #then write out the power spectrum goodness
        ft_mix = os.path.join(DIR, 'sample.' + str(sample) + '.' + sample_sfix, 'melodic_FTmix')
        ft_data = read_ftmix(ft_mix, vol)
        ft_out = os.path.join(oDIR, component, component + '_' + str_sample + '.FTmix.txt')
        write_ftmix(ft_data, ft_out)
        #and why not the time courses as well...
        ft_mix = os.path.join(DIR, 'sample.' + str(sample) + '.' + sample_sfix, 'melodic_mix')
        ft_data = read_ftmix(ft_mix, vol)
        ft_out = os.path.join(oDIR, component, component + '_' + str_sample + '.mix.txt')
        write_ftmix(ft_data, ft_out)
    #######find out which components didn't get matched and write them out as extras...
    if not os.path.exists(os.path.join(oDIR, 'unmatched')):
        os.makedirs(os.path.join(oDIR, 'unmatched'))
    for c in range(0, ndim):
        if not c in vol_selected:
            vol_sample = 'vol_' + str(c) + '_sample_' + str_sample
            recover_dict = {}
            print 'vol was not matched to a component: ', c
            print 'checking to see if there is an alternative match...'
            for template in templates:
                mcList = matched_components[template][0].tolist()
                if c in mcList:
                    recover_dict[template] = mcList.count(c)
            if not len(recover_dict.keys()) > 0:
                print 'no potential match found....'
                continue
            max_count = np.max(recover_dict.values())
            recovered = []
            for template in recover_dict.keys():
                if recover_dict[template] == max_count:
                    print 'potential match is: ', template
                    oFile = os.path.join(oDIR, template, template + '_' + str_sample + '.nii.gz')
                    if not os.path.exists(oFile):
                        print 'writing out recovered volume: ', template
                        recovered.append(vol_sample)
                        oFile = os.path.join(oDIR, 'unmatched', 'recovered_' + template + '_UNMATCHED_vol' + str(c) + '_sample' + str_sample + '.nii.gz')
                        fslroi = fsl.ExtractROI(in_file=src_map, roi_file=oFile, t_min=c, t_size=1)
                        fslroi.run()
                        break
            if not vol_sample in recovered:
                print 'Could not find match for: ', vol_sample, ' ....writing out unmatched volume.'
                oFile = os.path.join(oDIR, 'unmatched',  'unmatched_vol' + str(c) + '_sample' + str_sample + '.nii.gz')
                fslroi = fsl.ExtractROI(in_file=src_map, roi_file=oFile, t_min=c, t_size=1)
                fslroi.run()



print 'matching bootstrap sample components to templates...'
sample_queue = range(nsamples)

if __name__ == '__main__':
    pool = Pool(processes=ncores)
    pool.map(run_sample, sample_queue)





###in case the masking is needed, this preserves the values within the images...doesn't switch to bools
#set the z threshold for the positive and negative ranges
#trg_mask_pos = np.invert(trg_data[:,:,:,0] >= z)
#src_mask_pos = np.invert(src_data[:,:,:,0] >= z)
#trg_mask_neg = np.invert(trg_data[:,:,:,0] <= -(z))
#mask the arrays
#trg_data_masked_pos = np.ma.masked_array(trg_data[:,:,:,0], trg_mask_pos)
#src_data_masked_pos = np.ma.masked_array(src_data[:,:,:,0], src_mask_pos)
#trg_data_masked_neg = np.ma.masked_array(trg_data[:,:,:,0], trg_mask_neg)
