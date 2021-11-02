import os as os
import csv as csv
import copy as copy 

xnatFile = open('/Users/rmuetzel/Desktop/xnat_dump_forRonald_new.csv', 'r')
oFile = open('/Users/rmuetzel/Desktop/GenR_F9_MRI_xnat_mrsession_FORMATTED_30Nov2018_new.csv', 'w')

print 'Reading in CSV file from XNAT: ', xnatFile
csv_reader = csv.reader(xnatFile, delimiter=',')

xnat_mrsession_data = []

seq_dict = {}
r_dict = {}

for line in csv_reader:
	xnat_mrsession_data.append(line)
xnatFile.close()

print 'building sequence dictionary...'
for line in xnat_mrsession_data:
	if not line[0] == 'MR ID':
		scans=line[5].split('  ')
		print(line)
		rm_ws = True
		while rm_ws:
			try:
				scans.remove('')
			except:
				rm_ws = False
		for scan in range(0, len(scans)):
			if not scans[scan] == '':
				pfix = ('Exp', 'Insp', ' Exp', ' Insp' )
				if not scans[scan].startswith(pfix):
					seq = scans[scan]
					if scan+1 < len(scans):
						if scans[scan+1].startswith(pfix):
							seq = scans[scan] + scans[scan+1]
				seq = seq.replace('(1)', '').replace('(2)', '').replace('(3)', '')
				if not seq_dict.has_key(seq):
					seq_dict[seq] = 0
seq_dict['date'] = False

print 'Found the following sequences:'
print seq_dict

print 'populating r-dictionary with sequence information...'
for line in xnat_mrsession_data:
	if not line[0] == 'MR ID':
		r = line[0]
		date = line[1]
		print '.',
		r_dict[r] = copy.deepcopy(seq_dict)
		r_dict[r]['date'] = date
		scans=line[5].split('  ')
		rm_ws = True
		while rm_ws:
			try:
				scans.remove('')
			except:
				rm_ws = False
		for scan in range(0, len(scans)):
			if not scans[scan] == '':
				pfix = ('Exp', 'Insp', ' Exp', ' Insp' )
				if not scans[scan].startswith(pfix):
					seq = scans[scan]
					if scan+1 < len(scans):
						if scans[scan+1].startswith(pfix):
							seq = scans[scan] + scans[scan+1]
				seq2 = seq.replace('(1)', '').replace('(2)', '').replace('(3)', '')
				if seq.endswith('(1)'):
					count = 1
				elif seq.endswith('(2)'):
					count = 2
				elif seq.endswith('(3)'):
					count = 3
				else:
					count = 4
				r_dict[r][seq2] = count

print 'Found: ', len(r_dict.keys()), '  Rnumbers in XNAT...'
csv_data = []
headers = ['rnumber']
for i in seq_dict.keys():
	headers.append(i)
csv_data.append(headers)
for r in r_dict.keys():
	tmp = [r]
	for i in seq_dict.keys():
		tmp.append(r_dict[r][i])
	csv_data.append(tmp)

print 'Writing out CSV file: ', oFile
csv_writer = csv.writer(oFile, delimiter=',')
csv_writer.writerows(csv_data)
oFile.close()


print 'Finished.....Press any key to close this window.'
raw_input()
