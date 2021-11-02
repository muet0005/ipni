import os as os
import sys as sys
import csv as csv


btbl = open(sys.argv[1], 'r')
transp_btbl = open(sys.argv[2], 'w')


bvecs_data = []
bvecs = csv.reader(btbl, delimiter=' ')
for line in bvecs:
	bvecs_data.append(line)
ncols = len(bvecs_data[0])
nrows = len(bvecs_data)
transp_bvecs_data = []
for c in range(0, ncols):
	tmp = []
	for r in range(0, nrows):
		try:
			float(bvecs_data[r][c])
			tmp.append(bvecs_data[r][c])
		except:
			pass
	transp_bvecs_data.append(tmp)
transp_bvecs_data = filter(None, transp_bvecs_data)
transp_btbl_writer = csv.writer(transp_btbl, delimiter=' ', lineterminator='\n')
transp_btbl_writer.writerows(transp_bvecs_data)



btbl.close()
transp_btbl.close()