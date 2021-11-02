import pyxnat
from pyxnat import Interface
import sys as sys
import csv as csv

project = 'GenR_F9_MRI'
cfg_celeste = '/Users/rmuetzel/.pyxnat.cfg'
cfg_rietveld = '/Users/rmuetzel/.pyxnat.rv.cfg'
seq_names = ('T1_weighted')

def get_rs(project):
	central = Interface(config=cfg)
	constraint = [('xnat:SubjectData/PROJECT', '=', project)]
	f = central.select('xnat:SubjectData', ['xnat:SubjectData/SUBJECT_LABEL']).where(constraint)
	return f.items()

	
def get_mrsubjectid(r):
	central = Interface(config=cfg)
	constraint = [('xnat:SubjectData/SUBJECT_LABEL', '=', r)]
	f = central.select('xnat:SubjectData', ['xnat:SubjectData/SUBJECT_ID']).where(constraint)
	return f.items()[0][0]


def get_mrsessionid(subject_id):
	central = Interface(config=cfg)
	constraint = [('xnat:mrSessionData/SUBJECT_ID', '=', subject_id)]
	f = central.select('xnat:mrSessionData', ['xnat:mrSessionData/SESSION_ID']).where(constraint)
	return f.items()[0][0]

	
def get_scandata(session_id):
	central = Interface(config=cfg)
	constraint = [('xnat:mrScanData/IMAGE_SESSION_ID', '=', session_id)]
	f = central.select('xnat:mrScanData', ['xnat:mrScanData/SERIES_DESCRIPTION', 'xnat:mrScanData/ID', 'xnat:mrScanData/QUALITY', 'xnat:mrSessionData/DATE']).where(constraint)
	return f.items()


rs = get_rs(project)
csv_data = []
headers = ['rnumber', 'series_id', 'sequence', 'quality_rating', 'scan_date']
csv_data.append(headers)

for r in rs:
	r = r[0]
	subject_id = get_mrsubjectid(r)
	session_id = get_mrsessionid(subject_id)
	data = get_scandata(session_id)
	print 'RNumber: ', r,  ' Subject ID: ', subject_id, ' Session ID: ', session_id
	for i,j,k,l in data:
		if i.startswith(seq_names):
			tmp = [r, j, i, k,l]
			csv_data.append(tmp)
			print "Series ID: ", j, ' Series Description: ', i, ' Quality Rating: ', k, ' Date: ', l

csv_writer = csv.writer(open('/Users/rmuetzel/Desktop/xnat_qc_july11_2013.csv', 'w'), delimiter='|')
csv_writer.writerows(csv_data)