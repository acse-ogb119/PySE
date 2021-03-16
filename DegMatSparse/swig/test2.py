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
from MatSparse import *
from Numeric import arange
import time
from math import sqrt

mmat = MapMatSparse()
mmat.load('../../../data/A.m')

A = DegMatSparse(mmat)
B = FastMatSparse(mmat)

print "Done reading stuff"

n = B.n


print n

x = arange(n, typecode='d')
b = arange(n, typecode='d')

x *= 0; x += 1
b *= 0; b += 2

t = -time.time()
#for i in range(10):
#    x=A*x
#t += time.time()
#print t

print "Start timer"

t = -time.time()
for i in range(10):
    B.prod(x,b)
t += time.time()
print t


norm = 0.0
for i in range(n):
    norm += b[i]*b[i]


print sqrt(norm)
