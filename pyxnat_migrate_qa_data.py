import pyxnat
from pyxnat import Interface
import psycopg2

import sys as sys
import csv as csv

project = 'GenR_F9_MRI'
cfg_celeste = '/Users/rmuetzel/.pyxnat.cfg'
cfg_rietveld = '/Users/rmuetzel/.pyxnat.rv.cfg'
seq_names = ('T1_weighted')
conn = psycopg2.connect("dbname=xnat user=xnat01 password=BrainScan1")
cur = conn.cursor()

def get_rs(project, cfg):
	central = Interface(config=cfg)
	constraint = [('xnat:SubjectData/PROJECT', '=', project)]
	f = central.select('xnat:SubjectData', ['xnat:SubjectData/SUBJECT_LABEL']).where(constraint)
	return f.items()

	
def get_mrsubjectid(r, cfg):
	central = Interface(config=cfg)
	constraint = [('xnat:SubjectData/SUBJECT_LABEL', '=', r)]
	f = central.select('xnat:SubjectData', ['xnat:SubjectData/SUBJECT_ID']).where(constraint)
	return f.items()[0][0]


def get_mrsessionid(subject_id, cfg):
	central = Interface(config=cfg)
	constraint = [('xnat:mrSessionData/SUBJECT_ID', '=', subject_id)]
	f = central.select('xnat:mrSessionData', ['xnat:mrSessionData/SESSION_ID']).where(constraint)
	return f.items()[0][0]

	
def get_scandata(session_id, cfg):
	central = Interface(config=cfg)
	constraint = [('xnat:mrScanData/IMAGE_SESSION_ID', '=', session_id)]
	f = central.select('xnat:mrScanData', ['xnat:mrScanData/SERIES_DESCRIPTION', 'xnat:mrScanData/ID', 'xnat:mrScanData/QUALITY']).where(constraint)
	return f.items()


def update_qa(session_id, series_id, qa):
	cur.execute("update xnat_imagescandata set quality=\'" + qa + "\' where image_session_id=\'" + session_id + "\' and id=\'" + series_id + "\'")
	conn.commit()


rs = get_rs(project, cfg_celeste)

#for r in rs:
r = rs[1][0]
print r
subject_id = get_mrsubjectid(r, cfg_celeste)
subject_id_rv = get_mrsubjectid(r, cfg_rietveld)
session_id = get_mrsessionid(subject_id, cfg_celeste)
session_id_rv = get_mrsessionid(subject_id_rv, cfg_rietveld)
print subject_id, subject_id_rv, session_id, session_id_rv
data = get_scandata(session_id, cfg_celeste)
print 'RNumber: ', r,  ' Subject ID: ', subject_id, ' Session ID: ', session_id
for i,j,k in data:
	print "Series ID: ", j, ' Series Description: ', i, ' Quality Rating: ', k
	update_qa(session_id, j, 'holy_hell')
	#insert_qa(session_id_rv, j, k, cfg_rietveld)
