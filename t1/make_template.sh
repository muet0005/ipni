#!/bin/bash


DIR=/Volumes/mdraid1/mr_data_idc/feb2013_datalock/t1
oDIR=${DIR}/knicr166

idcs=(66 103 205 310 341 365 445 568 725 982 1021 1438 1620 1811 1868 1883 2023 2189 2206 2258 2423 2426 2497 2579 2650 2917 3042 3064 3124 3184 3340 3608 3789 3984 4174 4336 4496 4831 4914 5000 5027 5053 5168 5207 5616 5775 6191 6336 6397 6728 7023 7303 7418 7518 7769 7809 8013 9145 9358 84 97 108 140 172 241 286 323 345 415 459 485 530 552 669 672 734 898 963 990 1104 1143 1149 1350 1401 1404 1429 1567 1578 1656 1733 1842 2000 2056 2071 2121 2806 2873 2875 2910 3104 3128 3167 3245 3390 3586 3710 3794 3828 3864 4087 4166 4169 4202 4360 4450 4505 4510 4536 4554 4654 4849 4950 5034 5180 5192 5301 5313 5550 5700 5805 6322 6455 6509 6616 6663 7056 7153 7193 7302 7313 7336 7340 7393 7650 7738 7758 7960 8049 8056 8106 8195 8228 8232 8372 8492 8586 8598 8638 8805 8989 9037 9095 9458 9485 9489 9544)

#EXCLUDED 
#180 -> image not there?
#8255 -> Huge ventricles, CC agenesis

#motion
#1085 1414 1589 2155 231 2500 2841 4154 4186 4289 4629 4824 4919 507 5251 545 5587 6273 6398 6399 7234 818 8631 888 9588


if [ ! -d $oDIR ] ; then
	mkdir $oDIR
fi

for reg in fast lin1 lin2 nonlin1 nonlin2 nonlin3 ; do
	if [ ! -d $oDIR/${reg} ] ; then
		mkdir $oDIR/${reg}
	fi
done


##do the initial 6dof reg
for idc in ${idcs[*]} ; do
	echo $idc "linear 6dof"
	if [ ! -e $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_6dof.nii.gz -a ! -e $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
		#touch an isrunning so you can run on multiple cores..
		touch $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
		#resample from 512x512 -> 256x256, 0.9mm iso
		if [ ! -e $oDIR/lin1/t1_idc_${idc}_1mm.nii.gz ] ; then
			echo "Resampling"
			mri_convert ${DIR}/${idc}/t1_idc_${idc}.nii.gz -vs 0.9 0.9 0.9 $oDIR/lin1/t1_idc_${idc}_1mm.nii.gz
		fi
		#run an initial 12dof reg
		if [ ! -e $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_12dof.nii.gz ] ; then
			echo "Inital 12dof reg"
			${FSLDIR}/bin/flirt -in $oDIR/lin1/t1_idc_${idc}_1mm.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_1mm -out $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_12dof -omat $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_12dof.mat -dof 12
		fi
		if [ ! -e $oDIR/lin1/t1_idc_${idc}_1mm_brain.nii.gz ] ; then
			#invert the 12dof xfm
			echo "inverting xfm"
			$FSLDIR/bin/convert_xfm -omat $oDIR/lin1/MNI152_T1_1mm_12dof_2_t1_idc_${idc}.mat -inverse $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_12dof.mat
			#flirt the mni brain mask back to native space....use the mask with the vents/csf filled, and dilated 2x
			echo "aligning brain mask to native space"
			$FSLDIR/bin/flirt -in $oDIR/templates/MNI152_T1_1mm_brain_mask_CSFfilled_dil2_rlm -ref $oDIR/lin1/t1_idc_${idc}_1mm -out $oDIR/lin1/MNI152_T1_1mm_12dof_brain_mask_2_t1_idc_${idc} -applyxfm -init $oDIR/lin1/MNI152_T1_1mm_12dof_2_t1_idc_${idc}.mat
			#threshold the mni brain mask to get rid of interp
			$FSLDIR/bin/fslmaths $oDIR/lin1/MNI152_T1_1mm_12dof_brain_mask_2_t1_idc_${idc} -thr 0.5 -bin $oDIR/lin1/MNI152_T1_1mm_12dof_brain_mask_2_t1_idc_${idc}
			#bet the native brain with the mask from mni space
			echo "brain extraction"
			$FSLDIR/bin/fslmaths $oDIR/lin1/t1_idc_${idc}_1mm -mul $oDIR/lin1/MNI152_T1_1mm_12dof_brain_mask_2_t1_idc_${idc} $oDIR/lin1/t1_idc_${idc}_1mm_brain
		fi
		#intensity normalize the brain image with FAST
		if [ ! -e $oDIR/fast/t1_idc_${idc}_1mm_brain_restore.nii.gz  ] ; then
			echo "intensity normalization"
			$FSLDIR/bin/fast -b -B -o $oDIR/fast/t1_idc_${idc}_1mm_brain $oDIR/lin1/t1_idc_${idc}_1mm_brain
		fi
		#now run a 6dof reg with the instensity/bias corrected brain image
		if [ ! -e $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof.nii.gz ] ; then
			echo "Inital 6dof reg"
			$FSLDIR/bin/flirt -in $oDIR/fast/t1_idc_${idc}_1mm_brain_restore -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain -out $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof -omat $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof.mat -dof 6
		fi
		#correct the head image for intensity bias using the fast bias field output....NOTE...this will not correct the skull/non-brain!!!
		$FSLDIR/bin/fslmaths $oDIR/lin1/t1_idc_${idc}_1mm -div $oDIR/fast/t1_idc_${idc}_1mm_brain_bias $oDIR/fast/t1_idc_${idc}_1mm_restore
		#apply the 6dof xfm to the whole head 
		$FSLDIR/bin/flirt -in $oDIR/fast/t1_idc_${idc}_1mm_restore -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain -out $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_6dof -applyxfm -init $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm_brain_6dof.mat -dof 6
		rm $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
	fi
done


#make the first rough average brain from the linear reg.  first merge all images into a 4d file
if [ ! -e $oDIR/templates/knicr166_all_lin1_t1_1mm_brain_6dof.nii.gz ] ; then 
	cd $oDIR/lin1
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_lin1_t1_1mm_brain_6dof t1_idc_*_2_MNI152_T1_1mm_brain_6dof.nii.gz
fi
#now make the average over time
if [ ! -e $oDIR/templates/knicr166_lin1_t1_1mm_brain_6dof.nii.gz ] ; then
	cd $oDIR/lin1
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_lin1_t1_1mm_brain_6dof -Tmean $oDIR/templates/knicr166_lin1_t1_1mm_brain_6dof
fi

#do the same for the head...
if [ ! -e $oDIR/templates/knicr166_all_lin1_t1_1mm_6dof.nii.gz ] ; then
	cd $oDIR/lin1
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_lin1_t1_1mm_6dof t1_idc_*_2_MNI152_T1_1mm_6dof.nii.gz
fi
if [ ! -e $oDIR/templates/knicr166_lin1_t1_1mm_6dof.nii.gz ] ; then
	cd $oDIR/lin1
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_lin1_t1_1mm_6dof -Tmean $oDIR/templates/knicr166_lin1_t1_1mm_6dof
fi

#begin the 2nd linear reg step....this time use 12dof
for idc in ${idcs[*]} ; do
	echo $idc "linear 12dof"
	if [ ! -e $oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.nii.gz  -a ! -e $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
		touch $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
		echo "lin2 12dof reg"
		#reg to the 6dof aligned average brain from the first step above
		$FSLDIR/bin/flirt -in $oDIR/fast/t1_idc_${idc}_1mm_brain_restore -ref $oDIR/templates/knicr166_lin1_t1_1mm_brain_6dof -out $oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof -omat $oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.mat -dof 12
		#and for the head too
		$FSLDIR/bin/flirt -in $oDIR/fast/t1_idc_${idc}_1mm_restore -ref $oDIR/templates/knicr166_lin1_t1_1mm_brain_6dof -out $oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_12dof -applyxfm -init $oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.mat -dof 12
		rm $oDIR/lin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
	fi
done


#make the 2nd linear average, this time from the 12dof fit brains...
if [ ! -e $oDIR/templates/knicr166_all_lin2_t1_1mm_brain_12dof.nii.gz ] ; then 
	cd $oDIR/lin2
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_lin2_t1_1mm_brain_12dof t1_idc_*_2_knicr166_lin1_t1_1mm_brain_12dof.nii.gz
fi

if [ ! -e $oDIR/templates/knicr166_lin2_t1_1mm_brain_12dof.nii.gz ] ; then 
	cd $oDIR/lin2
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_lin2_t1_1mm_brain_12dof -Tmean $oDIR/templates/knicr166_lin2_t1_1mm_brain_12dof
fi

#and for the head images
if [ ! -e $oDIR/templates/knicr166_all_lin2_t1_1mm_12dof.nii.gz ] ; then 
	cd $oDIR/lin2
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_lin2_t1_1mm_12dof t1_idc_*_2_knicr166_lin1_t1_1mm_12dof.nii.gz
fi

if [ ! -e $oDIR/templates/knicr166_lin2_t1_1mm_12dof.nii.gz ] ; then 
	cd $oDIR/lin2
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_lin2_t1_1mm_12dof -Tmean $oDIR/templates/knicr166_lin2_t1_1mm_12dof
fi



#begin the first iteration of nonlinear warps to the 12dof aligned knicr brain
#use the previous 12dof xfm to initialize it
for idc in ${idcs[*]} ; do
	echo $idc "nonlinear 1"
	if [ ! -e $oDIR/nonlin1/t1_idc_${idc}_2_knicr166_lin2_t1_1mm_brain_nonlin1.nii.gz  -a ! -e $oDIR/nonlin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
		touch $oDIR/nonlin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
		$FSLDIR/bin/fnirt --in=$oDIR/fast/t1_idc_${idc}_1mm_brain_restore --ref=$oDIR/templates/knicr166_lin2_t1_1mm_brain_12dof --aff=$oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.mat --cout=$oDIR/nonlin1/t1_idc_${idc}_2_knicr166_lin2_t1_1mm_brain_nonlin1_cout --iout=$oDIR/nonlin1/t1_idc_${idc}_2_knicr166_lin2_t1_1mm_brain_nonlin1 --logout=$oDIR/nonlin1/t1_idc_${idc}_2_knicr166_lin2_t1_1mm_brain_nonlin1.log
		rm $oDIR/nonlin1/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
	fi
done

		
#make the first nonlinear average
if [ ! -e $oDIR/templates/knicr166_all_nonlin1_t1_1mm_brain.nii.gz ] ; then 
	cd $oDIR/nonlin1
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_nonlin1_t1_1mm_brain t1_idc_*_2_knicr166_lin2_t1_1mm_brain_nonlin1.nii.gz
fi

if [ ! -e $oDIR/templates/knicr166_nonlin1_t1_1mm_brain.nii.gz ] ; then 
	cd $oDIR/nonlin1
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_nonlin1_t1_1mm_brain -Tmean $oDIR/templates/knicr166_nonlin1_t1_1mm_brain
fi


#begin the second iteration of nonlinear warps...this time to the first nonlin1 aligned knicr brain
#use the previous 12dof xfm to initialize it
for idc in ${idcs[*]} ; do
	echo $idc "nonlinear2"
	if [ ! -e $oDIR/nonlin2/t1_idc_${idc}_2_knicr166_nonlin1_t1_1mm_brain_nonlin2.nii.gz  -a ! -e $oDIR/nonlin2/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
		touch $oDIR/nonlin2/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
		$FSLDIR/bin/fnirt --in=$oDIR/fast/t1_idc_${idc}_1mm_brain_restore --ref=$oDIR/templates/knicr166_nonlin1_t1_1mm_brain --aff=$oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.mat --cout=$oDIR/nonlin2/t1_idc_${idc}_2_knicr166_nonlin1_t1_1mm_brain_nonlin2_cout --iout=$oDIR/nonlin2/t1_idc_${idc}_2_knicr166_nonlin1_t1_1mm_brain_nonlin2 --logout=$oDIR/nonlin2/t1_idc_${idc}_2_knicr166_nonlin1_t1_1mm_brain_nonlin2.log
		rm $oDIR/nonlin2/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
	fi
done



#make the second nonlinear average
if [ ! -e $oDIR/templates/knicr166_all_nonlin2_t1_1mm_brain.nii.gz ] ; then 
	cd $oDIR/nonlin2
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_nonlin2_t1_1mm_brain t1_idc_*_2_knicr166_nonlin1_t1_1mm_brain_nonlin2.nii.gz
fi

if [ ! -e $oDIR/templates/knicr166_nonlin2_t1_1mm_brain.nii.gz  ] ; then 
	cd $oDIR/nonlin2
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_nonlin2_t1_1mm_brain -Tmean $oDIR/templates/knicr166_nonlin2_t1_1mm_brain
fi		



#begin the third iteration of nonlinear warps...this time to the second nonlin2 aligned knicr brain
#use the previous 12dof xfm to initialize it
for idc in ${idcs[*]} ; do
	echo $idc "nonlinear3"
	if [ ! -e $oDIR/nonlin3/t1_idc_${idc}_2_knicr166_nonlin2_t1_1mm_brain_nonlin3.nii.gz  -a ! -e $oDIR/nonlin3/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
		touch $oDIR/nonlin3/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
		$FSLDIR/bin/fnirt --in=$oDIR/fast/t1_idc_${idc}_1mm_brain_restore --ref=$oDIR/templates/knicr166_nonlin2_t1_1mm_brain --aff=$oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.mat --cout=$oDIR/nonlin3/t1_idc_${idc}_2_knicr166_nonlin2_t1_1mm_brain_nonlin3_cout --iout=$oDIR/nonlin3/t1_idc_${idc}_2_knicr166_nonlin2_t1_1mm_brain_nonlin3 --logout=$oDIR/nonlin3/t1_idc_${idc}_2_knicr166_nonlin2_t1_1mm_brain_nonlin3.log
		rm $oDIR/nonlin3/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
	fi
done


#make the third nonlinear average
if [ ! -e $oDIR/templates/knicr166_all_nonlin3_t1_1mm_brain.nii.gz ] ; then 
	cd $oDIR/nonlin3
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_nonlin3_t1_1mm_brain t1_idc_*_2_knicr166_nonlin2_t1_1mm_brain_nonlin3.nii.gz
fi

if [ ! -e $oDIR/templates/knicr166_nonlin3_t1_1mm_brain.nii.gz  ] ; then 
	cd $oDIR/nonlin3
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_nonlin3_t1_1mm_brain -Tmean $oDIR/templates/knicr166_nonlin3_t1_1mm_brain
fi	