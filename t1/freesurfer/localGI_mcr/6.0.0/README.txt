Steps to make the localGI functions from the recon-all stream to work in the matlab runtime compiler (MRC)


Note: all scripts should have "#RLM" at/around where I made changes...

1.) Compile the following matlab scripts using the matlab runtime compiler:
	compute_lgi    
	find_corresponding_center_FSformat   
	make_outer_surface  
	make_roi_paths
	
	The output is simply a "run_${scriptname}.sh" file that loads the proper MRC env, and points to the compiled code, ${scriptname}.
	
	The "run" script takes the same args as the matlab code, with 1 exception....the first argument must be the MRC path (see step 3b below)
	
	**Before compiling!!!  A minor edit was made to the scripts that have non-string inputs (e.g., radius).  The arguments passed through the run script to the compiled code are all in string format as far as the matlab compiler is concerned.  So, in the scripts that needed it, the non-string variables were converted to double with:
		v = str2double(v), e.g., radius = str2double(radius)
		
	The edited scripts are included in the matlab folder of the zip archive.
	
	The compiled code is in the "distrib" folder of the zip archive.
	
2.) Edits to the recon-all scripts  (See recon-all_genr in the enclosed zip archive)
	Specifically, two points:
		a.) Add flags, -localGI_mrc, -lgi_mrc, and -lGI_mrc
		b.) associate a call to mris_compute_lgi_mrc (and edited version of mris_compute_lgi)

	Very basically, I searched for anything related to localGI or gyrifaction and added sections for the MRC env

3.) Edit the mris_compute_lgi (see step 2b in above, also included in zip archive)
	There are 4 matlab calls that need to be edited. the original calls to matlab were left in.
	Two variables are set at the top of this script:
		a.) fs_mrc  must be set (perhaps this can be added to the setup scripts?) This is where the compiled matlab scripts live.  Easiest case scenario, these are packaged with the $FREESURFER_HOME/bin folder..
	 	b.) MRC -->this is where the matlab runtime compiler lives.  Perhaps this could be added to the matlab chunk in the setup scripts...
	
	Very basically, I replaced the calls to the matlab scripts with calls more similar to what is done with plain old binaries in $FREESURFER_HOME/bin