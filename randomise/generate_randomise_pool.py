import os as os
import subprocess as sp

poolname='list.pool'
design_basename = 'design_mainfx'

in_imgs = ['dr_stage2_merged_pe_ic0001', 'dr_stage2_merged_pe_ic0002','dr_stage2_merged_pe_ic0003','dr_stage2_merged_pe_ic0004',
'dr_stage2_merged_pe_ic0005','dr_stage2_merged_pe_ic0006','dr_stage2_merged_pe_ic0007','dr_stage2_merged_pe_ic0008','dr_stage2_merged_pe_ic0016',
'dr_stage2_merged_pe_ic0017','dr_stage2_merged_pe_ic0018','dr_stage2_merged_pe_ic0019','dr_stage2_merged_pe_ic0020','dr_stage2_merged_pe_ic0021']

nperms=10000
ncores=100


poolFile = os.path.join(os.getcwd(), poolname + '.data')

nperms_per_core = int((nperms / ncores) + 1)

poolData = open(poolFile, 'w')

for i in in_imgs:
    for n in range(1, ncores+1):
        line = i + ' ' + design_basename + ' ' +str(n) + ' ' + str(nperms_per_core)
        poolData.write(line + '\n')
poolData.close()

sp.call(["stopos purge -p " + poolname], shell=True)
sp.call(["stopos create -p " +poolname], shell=True)
sp.call(["stopos add -p " + poolname +  " " + poolFile], shell=True)
