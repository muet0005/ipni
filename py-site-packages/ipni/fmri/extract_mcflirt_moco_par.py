import os as os
import sys as sys
import csv as csv
import numpy as np
import math as math
import argparse as argparse


class eval_mcflirt_moco():
    """
    Class to read in the .par motion correction file from FSLs MCFLIRT and
    return three items (see extract_moco_summary function)
    """
    def __init__(self, feat_dir):
        """
        Initialize the class.  Takes one argument, feat_dir.  THis is the folder where the feat output reside, and should contain an "mc"
        folder that holds the moco output.
        """
        self.feat_dir = feat_dir
    
    def read_mcflirt_moco(self, **kwargs):
        """
        Will read in the .par file. It will assume the name of the file is called, "prefiltered_func_data_mcf.par" and that
        it is within the feat_dir and a subfolder called "mc"
        i.e.,   feat_dir/mc/prefiltered_func_data_mcf.par
        You can specifiy par_file= to change the name of hte par file...
        """
        par_file = 'prefiltered_func_data_mcf.par'
        for i in kwargs.keys():
			if i == 'par_file':
			    par_file = kwargs[i]
        mcDir = os.path.join(self.feat_dir, 'mc')
        if not os.path.isdir(mcDir):
        	print 'No_Feat_DIR'
        	sys.exit(0)
        parFile = open(os.path.join(mcDir, par_file), 'r')
        csv_parFile = csv.reader(parFile, delimiter=' ')
        data_list = []
        for line in csv_parFile:
        	rmWS = True
        	while rmWS:
        		try:
        			line.remove('')
        		except:
        			rmWS = False
        	data_list.append(line)
        parFile.close()
        self.trs = len(data_list)
        par = len(data_list[0])
        if not par == 6:
            print "Did not find 6 motion parameters (3 rotation, 3 translation)...must exit."
            sys.exit(0)
        self.mc_pars = np.zeros(shape=(self.trs, par))
        for i in range(0, len(data_list)):
        	for j in range(0, par):
        		self.mc_pars[i] = data_list[i]
    
    def extract_moco_summary(self, **kwargs):
        """
        will evaluate the moco parameters and return 3 pieces of information:
        
        1.) True/False on whether the main exclude crieteria were met for abs translations
        2.) The max abs translation
        3.) the mean relative translation
        
        Default max translation is set to 3mm
        Pass arg maxtrans= to set a different value
        """
        maxtrans = 3
        for i in kwargs.keys():
			if i == 'maxtrans':
			    maxtrans = int(kwargs[i])
        excl_dict = {}
        moco_abs_info = []
        moco_rel_info = []
        for i in range(0, self.trs):
        	#first figure out max absolute displacement
        	# Get the x,y,z translations in mm
        	x,y,z = self.mc_pars[i][3],self.mc_pars[i][4],self.mc_pars[i][5]
        	# compute the absolute displacement
        	d_abs = np.sqrt((x**2) + (y**2) + (z**2))
        	moco_abs_info.append(d_abs)
        	if d_abs > maxtrans:
        		if excl_dict.has_key(self.feat_dir):
        			excl_dict[self.feat_dir].append(d_abs)
        		else:
        			excl_dict[self.feat_dir] = [d_abs]
        	#next, compute the relative displacement for each tr
        	if not i == 0:
        	    #first get the xyz translations in mm for the previous TR
        		x_ref,y_ref,z_ref = self.mc_pars[i-1][3],self.mc_pars[i-1][4],self.mc_pars[i-1][5]
        		#compute the relative displacement
        		d_rel = np.sqrt(((x-x_ref)**2) + ((y-y_ref)**2) + ((z-z_ref)**2))
        		moco_rel_info.append(d_rel)
        exclude_case = False
        if excl_dict.has_key(self.feat_dir):
        	exclude_case = True
        return exclude_case, np.max(moco_abs_info), np.mean(moco_rel_info)
    
    def extract_rms(self, rms_file):
        """
        extrac the rms value from either the abs or rms mean files from McFlirt
        """
        mcDir = os.path.join(self.feat_dir, 'mc')
        if not os.path.isdir(mcDir):
        	print 'No_Feat_DIR'
        	sys.exit(0)
        rmsFile = open(os.path.join(mcDir, rms_file), 'r')
        csv_parFile = csv.reader(rmsFile, delimiter=' ')
        rmsData = []
        for line in rmsFile:
            rmsData.append(line)
        if len(rmsData) == 1:
            return float(rmsData[0])
        else:
            print 'RMS file: ', rms_file, ' Has multiple lines....cannot use this function to parse this file.'
            return False
    



if __name__ == '__main__':
    #If run via the cmd line....
    #First parse the command-line arguments that we need.
    #set up the parser
    parser = argparse.ArgumentParser(description='Parse MCFLIRT moco files to extract information')
    #add the args you want parsed
    parser.add_argument("-f", "--feat_dir", help="Location of subject's feat folders", required=True, nargs=1)
    parser.add_argument("--maxtra", help="maximum amount of translations (in mm)",	required=False, nargs=1)
    #parser.add_argument("--genr", help="assume the outdir, nifti pfix and sfix are ordered in a certain way", required=False, action="store_true")
    #parse the args
    args = parser.parse_args()
    if args.maxtra:
    	maxtra = float(args.maxtra[0])
    else:
    	maxtra = 3
    feat_dir = args.feat_dir[0]
    x = eval_mcflirt_moco(feat_dir)
    x.read_mcflirt_moco()
    exclude_case,max_abs_trans,mean_rel_trans = x.extract_moco_summary()
    print "EXCLUDE CASE: ", exclude_case
    print "Max abs translation: ", max_abs_trans
    print "Mean rel translation: ", mean_rel_trans









  
#if math.degrees(np.max(np.abs(mc_pars[:,0]))) > maxrot:
#	if not excl_dict.has_key(feat_dir):
#		excl_dict[feat_dir] = [math.degrees(np.max(mc_pars[:,0])), 'xrot']
#if math.degrees(np.max(np.abs(mc_pars[:,1]))) > maxrot:
#	if not excl_dict.has_key(feat_dir):
#		excl_dict[feat_dir] = [math.degrees(np.max(mc_pars[:,1])), 'yrot']
#if math.degrees(np.max(np.abs(mc_pars[:,2]))) > maxrot:
#	if not excl_dict.has_key(feat_dir):
#		excl_dict[feat_dir] = [math.degrees(np.max(mc_pars[:,2])), 'zrot']
#if np.max(np.abs(mc_pars[:,3])) > maxtra:
#	if not excl_dict.has_key(feat_dir):
#		excl_dict[feat_dir] = [np.max(mc_pars[:,3]), 'xtrans']
#if np.max(np.abs(mc_pars[:,4])) > maxtra:
#	if not excl_dict.has_key(feat_dir):
#		excl_dict[feat_dir] = [np.max(mc_pars[:,4]), 'ytrans']
#if np.max(np.abs(mc_pars[:,5])) > maxtra:
#	if not excl_dict.has_key(feat_dir):
#		excl_dict[feat_dir] = [np.max(mc_pars[:,5]), 'ytrans']				
#print 'Max X rotation (rad, deg): ', np.max(mc_pars[:,0]), math.degrees(np.max(mc_pars[:,0]))
#print 'Max Y rotation (rad, deg): ', np.max(mc_pars[:,1]), math.degrees(np.max(mc_pars[:,1]))
#print 'Max Z rotation (rad, deg): ', np.max(mc_pars[:,2]), math.degrees(np.max(mc_pars[:,2]))
#print 'Max X translation (mm): ', np.max(mc_pars[:,3])
#print 'Max Y translation (mm): ', np.max(mc_pars[:,4])
#print 'Max Z translation (mm): ', np.max(mc_pars[:,5])