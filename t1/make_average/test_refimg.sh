#!/bin/bash


DIR=/Volumes/rbraid/mr_data_idc/feb2013_datalock/t1
oDIR=${DIR}/knicr166

idcs=(66 103 205 310 341 365 445 568 725 982 1021 1438 1620 1811 1868 1883 2023 2189 2206 2258 2423 2426 2497 2579 2650 2917 3042 3064 3124 3184 3340 3608 3789 3984 4174 4336 4496 4831 4914 5000 5027 5053 5168 5207 5616 5775 6191 6336 6397 6728 7023 7303 7418 7518 7769 7809 8013 9145 9358 84 97 108 140 172 241 286 323 345 415 459 485 530 552 669 672 734 898 963 990 1104 1143 1149 1350 1401 1404 1429 1567 1578 1656 1733 1842 2000 2056 2071 2121 2806 2873 2875 2910 3104 3128 3167 3245 3390 3586 3710 3794 3828 3864 4087 4166 4169 4202 4360 4450 4505 4510 4536 4554 4654 4849 4950 5034 5180 5192 5301 5313 5550 5700 5805 6322 6455 6509 6616 6663 7056 7153 7193 7302 7313 7336 7340 7393 7650 7738 7758 7960 8049 8056 8106 8195 8228 8232 8372 8492 8586 8598 8638 8805 8989 9037 9095 9458 9485 9489 9544)


#make the third nonlinear average
if [ ! -e $oDIR/templates/knicr166_all_nonlin3_t1_1mm.nii.gz ] ; then 
	cd $oDIR/nonlin3
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_nonlin3_t1_1mm t1_idc_*_2_knicr166_nonlin2_t1_1mm_nonlin3.nii.gz
fi
if [ ! -e $oDIR/templates/knicr166_nonlin3_t1_1mm.nii.gz  ] ; then 
	cd $oDIR/nonlin3
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_nonlin3_t1_1mm -Tmean $oDIR/templates/knicr166_nonlin3_t1_1mm
fi


for idc in ${idcs[*]} ; do
	if [ ! -e $oDIR/nonlin4/t1_idc_${idc}_2_knicr166_nonlin3_t1_1mm_nonlin4.nii.gz  -a ! -e $oDIR/nonlin4/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
		touch $oDIR/nonlin4/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
		echo $idc nonlin4 head
		$FSLDIR/bin/fnirt --in=$oDIR/fast/t1_idc_${idc}_1mm_restore --ref=$oDIR/templates/knicr166_nonlin3_t1_1mm --aff=$oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.mat --cout=$oDIR/nonlin4/t1_idc_${idc}_2_knicr166_nonlin4_t1_1mm_nonlin4_cout --iout=$oDIR/nonlin4/t1_idc_${idc}_2_knicr166_nonlin3_t1_1mm_nonlin4 --logout=$oDIR/nonlin4/t1_idc_${idc}_2_knicr166_nonlin3_t1_1mm_nonlin4.log --refmask=$oDIR/templates/knicr166_nonlin3_t1_1mm_brain_mask_rlm_edits.nii.gz
		rm $oDIR/nonlin4/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
	fi
done

#make the third nonlinear average
if [ ! -e $oDIR/templates/knicr166_all_nonlin4_t1_1mm.nii.gz ] ; then 
	cd $oDIR/nonlin4
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_nonlin4_t1_1mm t1_idc_*_2_knicr166_nonlin3_t1_1mm_nonlin4.nii.gz
fi
if [ ! -e $oDIR/templates/knicr166_nonlin4_t1_1mm.nii.gz  ] ; then 
	cd $oDIR/nonlin4
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_nonlin4_t1_1mm -Tmean $oDIR/templates/knicr166_nonlin4_t1_1mm
fi


for idc in ${idcs[*]} ; do
	if [ ! -e $oDIR/nonlin5/t1_idc_${idc}_2_knicr166_nonlin4_t1_1mm_nonlin5.nii.gz  -a ! -e $oDIR/nonlin5/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning ] ; then
		touch $oDIR/nonlin5/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
		echo $idc nonlin5 head
		$FSLDIR/bin/fnirt --in=$oDIR/fast/t1_idc_${idc}_1mm_restore --ref=$oDIR/templates/knicr166_nonlin4_t1_1mm --aff=$oDIR/lin2/t1_idc_${idc}_2_knicr166_lin1_t1_1mm_brain_12dof.mat --cout=$oDIR/nonlin5/t1_idc_${idc}_2_knicr166_nonlin4_t1_1mm_nonlin5_cout --iout=$oDIR/nonlin5/t1_idc_${idc}_2_knicr166_nonlin4_t1_1mm_nonlin5 --logout=$oDIR/nonlin5/t1_idc_${idc}_2_knicr166_nonlin4_t1_1mm_nonlin5.log --refmask=$oDIR/templates/knicr166_nonlin3_t1_1mm_brain_mask_rlm_edits.nii.gz
		rm $oDIR/nonlin5/t1_idc_${idc}_2_MNI152_T1_1mm.isrunning
	fi
done




#make the fifth nonlinear average
if [ ! -e $oDIR/templates/knicr166_all_nonlin5_t1_1mm.nii.gz ] ; then 
	cd $oDIR/nonlin5
	${FSLDIR}/bin/fslmerge -t $oDIR/templates/knicr166_all_nonlin5_t1_1mm t1_idc_*_2_knicr166_nonlin4_t1_1mm_nonlin5.nii.gz
fi
if [ ! -e $oDIR/templates/knicr166_nonlin5_t1_1mm.nii.gz  ] ; then 
	cd $oDIR/nonlin5
	${FSLDIR}/bin/fslmaths $oDIR/templates/knicr166_all_nonlin5_t1_1mm -Tmean $oDIR/templates/knicr166_nonlin5_t1_1mm
fi
