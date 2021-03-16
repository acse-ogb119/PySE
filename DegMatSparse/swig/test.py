#
#     (c) Copyright 2005
#     Author: Ola Skavhaug
#     Simula Research Laboratory AS
#     
#     This file is part of DegMatSparse.
#
#     DegMatSparse is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     DegMatSparse is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with DegMatSparse; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


from DP import *
from Numeric import *
from MatSparse import *
import time
from math import sqrt

nloops = 10
matrixfile = '../../../data/A.m'
#matrixfile = 'f.m'

#-----------------------------------------------------------------------
# Pretty print the timing results
#-----------------------------------------------------------------------
def report(name="Dummy", t0=-1, t1=-1, norm=-1):
    pstr = "Testing %s" % (name)
    print "%s\n%s\n%s" % (len(pstr)*'-', pstr, len(pstr)*'-')
    if (t0 > 0):
        print "The time of %d A*x (with allocation) is: %4.3e" % (nloops, t0)
    if (t1 > 0):
        print "The time of %d A*x (without allocation) is: %4.3e" % (nloops, t1)
    if (norm > 0):
        print "... and the norm is %f " % norm
    print ""
    
def anorm(x):
    n = 0.0
    for i in range(len(x)): n += x[i]*x[i]
    return sqrt(n)


#-----------------------------------------------------------------------
# Use a Diffpack::MatSparse to read the file (slow) 
#-----------------------------------------------------------------------


s = MatSparse_double()
s.load(matrixfile)
nrows = s.getNoRows()
x_dp = Vec_double(nrows)
x_dp.fill(1.0)
b_dp = Vec_double(nrows)
b_dp.fill(0.0)

A=s

t_dp = -time.time(); 
for i in range(nloops): A.prod(x_dp,b_dp)
t_dp += time.time()

report(name="MatSparse_double", t1=t_dp)



#-----------------------------------------------------------------------
# Testing the FastMatSparse implementation (pointer based)
#-----------------------------------------------------------------------

# Construct the matrix and vectors
p = s.pattern()
m2i = arange(nrows)
A = DegMatSparse()
A.values = s.getValues()
A.row = p.getIrow()
A.m2i = m2i
A.col = p.getJcol()
A.setOffset(1)

x = arange(nrows, typecode='d')
b = arange(nrows, typecode='d')
x *= 0; x += 1; b *= 0;

# Time the matrix-vector products
t0 = -time.time()
for i in range(nloops): b = A*x
t0 += time.time()

t1 = -time.time(); 
for i in range(nloops): A.prod2(x,b)
t1 += time.time()

report(name="DegMatSparse", t0=t0, t1=t1, norm=anorm(b))

#-----------------------------------------------------------------------
# Testing the FastMatSparse implementation (pointer based)
#-----------------------------------------------------------------------

B = FastMatSparse(nrows)
B.values = s.getValues()
B.row = p.getIrow()
B.col = p.getJcol()

B.setOffset(1)

t = -time.time()
for i in range(nloops):
    B.prod(x,b)
t += time.time()


report(name="FastMatSparse", t1=t, norm=anorm(b))

#-----------------------------------------------------------------------
# Changing the DP::MatSparse_double, affecting FastMatSparse
#-----------------------------------------------------------------------


print """
Change the values in the DP-matrix (by using the Numeric.array typemap
pointer). Compute the matrix-vector product by using FastMatSparse. This
changes the norm of the result vector. I.e., the pointers are set up properly.
"""

vals = s.getValues()
vals *= 2.0
B.prod(x,b)
report(name="FastMatSparse, take II", norm=anorm(b))
vals *= 0.5;
B.prod(x,b)
report(name="FastMatSparse, take III", norm=anorm(b))

print """Conclusion: 
The run time reduction (t/t_dp) is %f""" % (t/t_dp)

