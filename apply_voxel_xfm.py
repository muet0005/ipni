import os as os
import sys as sys
import numpy as np


xfm =   np.array([[1.121683074, 4.887416731e-05, -0.001616224127, -10.84774627],
        [0.00021499525, 1.055417576, -0.03495203452, -3.695648498],
        [0.004571149443, 0.004456808645, 1.117034189, -11.02037863],
        [0,  0,  0,  1]])
        
        
x = 45 ; y = 45 ; z = 45

X = np.array([x, y, z, 1])

Y = X*xfm

y1 = Y[1, 1]
x1 = Y[2, 1]
z1 = Y[3, 1]





xfm_inv = np.linalg.inv(xfm)

Xsolv = xfm_inv*Y