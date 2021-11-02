import knicr.structural.freesurfer_utils as fsu
import MySQLdb as mysqldb
import os as os
from multiprocessing import Pool



#where the data live
SUBJECTS_DIR = '/Volumes/rbraid/mr_data_idc/aug2013_final/freesurfer'

#the pull path to your label (without the ".label" at the end)
labeldir = '/Volumes/rbraid/mr_data_idc/aug2013_final/freesurfer/laura/dys'
labels = ['LH_VIS_rosmidfront', 'LH_VIS_parsop', 'LH_VIS_suptemp', 'LH_VIS_prec', 'LH_VIS_latorbfront', 'RH_VIS_prec', 'LH_MEM_prec', 'LH_MEM_precentral', 'LH_MEM_suppar', 'RH_MEM_prec', 'RH_MEM_postcentral', 'LH_ATT_superior_par']

#specify the measure
ctxmeasures=['thickness']

#specify table name
tbl = 'NEPSY_thickness'

#specify the subjects to run on (default is all)
#subject_list=['1', '100']

#your microsection number
mydb = '070006'

#do you want to run in parallel:
#True or False
parallel = True
#how many cores would you like to use?
#max is 20, but that's only if no one else is using the system.
#probably 5-10 is best...
ncores = 6

#below here, you shouldn't need to edit anything, unless you want to change the selection of subjects etc.
# if you type help(fsu) you will get all options of changes that can be made.


#call the class
s = fsu.qdec_clusters(SUBJECTS_DIR)

#connect to the mysql database
con = mysqldb.connect(db = mydb, read_default_file=os.path.join(os.getenv("HOME"), '.my.cnf'))
cursor = con.cursor()


def run_extract_stats(l):
	lhsw=('lh', 'LH', 'Lh')
	rhsw=('rh', 'RH', 'Rh')
	if l.startswith(lhsw):
		hemi='lh'
	elif l.startswith(rhsw):
		hemi='rh'
	else:
		hemi=False
	if not hemi:
		print 'labelname should start with lh or rh'			
		break
	label=os.path.join(labeldir, l)
	#map the labels to native space
	s.mapQdecLabel2Subjects(label, hemi)
	for ctxmeasure in ctxmeasures:
		#make the .stats files (if you want to specify specific subjects add: , subjects=subject_list)
		s.computeQdecLabelStats(label, hemi, measure=ctxmeasure)
		#dump the data in the .stats files into the database
		s.dumpQdecLabelStatsMysql(cursor, label, hemi, measure=ctxmeasure, table=tbl)


if not __name__ == '__main__' and parallel:
    parallel = False

if not parallel:
    #just run one at a time
    for l in labels:
        run_extract_stats(l)
else:
    #run on as many cores as specified
    pool = Pool(processes=ncores)
    pool.map(run_extract_stats, labels)


