#!/usr/bin/env python

from pyFDM import *
from pyFDM import NeumanBoundary

#import sys,os
#if os.environ.get('NUMPYARRAY','') == 'Numeric' or sys.modules.has_key('Numeric'):
#    from Numeric import *
#else:

from numarray import *

from Field import *
from Grid import *

def bc(x):
    if x[0] < 0.5:
        return 1
    else:
        return 0

def nc(x):
    if x[0] < 0.5:
        return 4*pi
    else:
        return 3.5*pi
    
def initc(*x):
    return cos(x[0] * pi * 0.5) + sin(x[0] * 4 * pi)

t=0.0
T=0.005
ts=150
dt = (T - t)/ts

l=0.0
L=1.0
n=50

