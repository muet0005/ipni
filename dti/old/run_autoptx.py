import os as os
import sys as sys
import knicr.dti.dti_utils as kdti

dtiDir = sys.argv[1]
subj = sys.argv[2]

autoPtxLib = '/home/genr/software/autptx'

kdti_aptx = kdti.knicrAutoPtx(dtiDir, subj, autoPtxLib)
kdti_aptx.autoPtx_1_preproc()
kdti_aptx.autoPtx_2_launchTractography()
