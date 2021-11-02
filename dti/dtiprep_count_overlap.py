import os as os
import sys as sys
import MySQLdb as mysqldb
import csv as csv
from collections import Counter

#mysql table name
tbl_in = 'dtiprep_qa_06dec2014'
tbl_out = 'dtiprep_qa_06dec2014_overlap'
#mysql database name
mydb = 'dti'

##connect to the mysql server
con = mysqldb.connect(db = mydb, read_default_file=os.path.join(os.getenv("HOME"), '.my.cnf'))
##set up a cursor to use
cursor = con.cursor()

#tell the code which column is the id column in the db
id_col = 'idc'

#select all of the subjects from the db
cursor.execute("""select %s from %s""" % (id_col, tbl_in))
rs = cursor.fetchall()

#tell it how many DTI dwis you have
nvols = 38

overlap_data = {}

#rs = ['R114204', 'R114114']
for r in rs:
    r = r[0]
    print r
    r_slice_data = {}
    #sys.exit(0)
    #make the zeropadded vol numbers
    for v in range(1, nvols+1):
        if v < 10:
            vstr = '000' + str(v)
        elif v < 100:
            vstr = '00' + str(v)
        elif v < 1000:
            vstr = '0' + str(v)
        else:
            vstr = str(v)
        col = 'GradDir_' + vstr + '_BadSliceNums'
        #pull the mysql data for each one that exists
        try:
            cursor.execute("""select %s from %s where %s=\'%s\'""" % (col, tbl_in, id_col, r))
            d = cursor.fetchone()
            #print col, d
            r_slice_data[col] = d[0].split(',')
        except:
            pass
    #get a list of all the dictionary keys for this subject
    dkeys = r_slice_data.keys()
    #count the number of them
    ndkeys = len(dkeys)
    overlap = []
    #loop over all unique pairwise combinations of volumes
    for a in range(0, ndkeys):
        for b in range(0, ndkeys):
            b += a+1
            if not b < ndkeys:
                continue
            #compute overlap between the lists
            o = set(r_slice_data[dkeys[a]]) & set(r_slice_data[dkeys[b]])
            print a,b, '    ', o
            if len(o) > 0:
                for x in o:
                    #append overlapping slice numbers into the overlap list
                    overlap.append(x)
    overlap = Counter(overlap)
    final_overlap = {}
    #now make a simple count of how many slices had overlap
    #1 = overlap between two volumes
    #2 = overlap between 3 volumes and so on...
    for x in overlap.keys():
        if final_overlap.has_key(overlap[x]):
            final_overlap[overlap[x]] += 1
        else:
            final_overlap[overlap[x]] = 1
    overlap_data[r] = final_overlap
    #last, dump it into a mysql db
    try:
        cursor.execute("""create table %s (%s varchar(10))""" % (tbl_out, id_col))
    except (mysqldb._mysql.OperationalError):
        pass
    #check the table to see if the subject exists
    idval = r
    cursor.execute("""select %s from %s where %s=\'%s\'""" % (id_col, tbl_out, id_col, idval))
    pn_exist = cursor.fetchone()
    if not pn_exist:
        cursor.execute("""insert into %s (%s) values (\'%s\')""" % (tbl_out, id_col, idval))
    for x in final_overlap.keys():
        col = 'slice_volume_overlap_' + str(x)
        val = final_overlap[x]
        try:
            #see if that column exists...if not...make it.
            cursor.execute("""alter table %s add column (%s float)""" % (tbl_out, col))
        except:
            pass
        #update the data base with the data value...
        cursor.execute("""update %s set %s=%s where %s=\'%s\'""" % (tbl_out, col, val, id_col, idval))
    print final_overlap
    print overlap
    for i in r_slice_data.keys():
        print i, '   ', r_slice_data[i]
        