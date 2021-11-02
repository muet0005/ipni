import os as os
import sys as sys
import knicr.dti.dti_utils as kdti

dtiDir = sys.argv[1]
subj = sys.argv[2]


kdti_prep = kdti.knicrDTIprep(dtiDir, subj)
kdti_prep.bedpost()