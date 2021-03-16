#!/usr/bin/env python
# -*- coding: iso8859-1 -*-
#
# (c) Copyright 2005
#     Author: Åsmund Ødegård
#     Simula Research Laboratory AS
#     
#     This file is part of PyFDM.
#
#     PyFDM is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     PyFDM is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with PyFDM; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#     Example usage of PyFDM.

from pyFDM import *
import getopt, sys
from math import exp,sin,cos,pi,sqrt
import time

def runstandalone():
    T=0.0025
    n=20
    m=40

    parasize = pypar.size()
    myrank = pypar.rank()

    # If we are parallel, and are using mpich MPI, we need to distribute the values
    # for the options, because they are only read on the first cpu when using mpich.
    
    if myrank == 0:
        opts = sys.argv[1:]
        numopts = len(opts)
        i = 0
        while i < numopts:
            opt = opts[i]
            if opt == "-t":
                T = float(opts[i+1])
                i += 1
            if opt == "-n":
                n = int(opts[i+1])
                i += 1
            if opt == "-m":
                m = int(opts[i+1])
                i += 1
            i += 1
    
    if parasize > 1:
        allopts = Numeric.zeros(3,typecode='d')
        if myrank == 0:
            allopts = Numeric.array((T,n,m))
        pypar.broadcast(allopts,0)
        if myrank != 0:
            T = float(allopts[0])
            n = int(allopts[1])
            m = int(allopts[2])

    # Solve, return the solution and a reference to the grid
    solution,g,T = solve(T,n,m)

    # The analytical solution
    analytical = lambda x,y: exp(-T)*numarray.sin(pi*x)*numarray.cos(pi*y)

    u_a = Field(g)
    u_a.fill_vec(analytical)

    # The error
    error = u_a - solution
    #error.plot(title='Error')
    #solution.plot(title='Solution')
    #u_a.plot(title='Analytical solution')

    errornorm = sqrt(g.dx*g.dy*error.inner(error))
    if g.myrank == 0 or g.myrank == -1:
        print "The norm of the error is: ",errornorm

    if solution.isParallel:
        pypar.finalize()

def usage(T,n,m):
    print "-t T: solve from 0 to T (default %s)" % (T)
    print "-n steps: Split time in so many steps (default %s)" % (n)
    print "-m div: Split each direction in so many nods (default %s)" % (m)


def solve(T=1.0, n=1000, m=100):
    # solve the problem from 0 to T, dt = 1/1000, dx = 1/100
    myrank = pypar.rank()

    if myrank <= 0:
        print "Solve the problem from 0 to %s in %s steps, %s grid-nodes" % (T,n,m)

    # Running time:
    rt = 0.0

    # compute dt:
    dt = T/(1.0*n)

    # Create a grid
    g = Grid((m,m),([0.0,1.0],[0.0,1.0]))

    dx = g.dx
    h = dt/(dx*dx)

    # Functions for source, coefficients, boundary, initial condition
    dt_f = lambda x,y: dt*exp(-rt)*(-pi*numarray.cos(pi*x)*numarray.cos(pi*y) \
             + pi*numarray.sin(pi*x)*numarray.sin(pi*y) \
             + (2.0*(x + y)*pi*pi - 1.0)*numarray.sin(pi*x)*numarray.cos(pi*y))
    # need only two k-functions, as dx is the same in both directions!
    k_p = lambda x,y: x + y + 0.5*dx
    k_m = lambda x,y: x + y - 0.5*dx
    # boundary fu:
    bf = lambda x,y: exp(-rt)*Numeric.sin(pi*x)*Numeric.cos(pi*y)
    # initical condition:
    initf = lambda x,y: numarray.sin(pi*x)*numarray.cos(pi*y)


    # build the stencil
    s = Stencil(nsd=2,varcoeff=True)
    s.addNode((0,-1),[lambda *x: h*k_p(*x)])
    s.addNode((-1,0),[lambda *x: h*k_p(*x)])
    s.addNode((0,0),[lambda *x: 1.0 - 2.0*h*(k_m(*x) + k_p(*x))])
    s.addNode((1,0),[lambda *x: h*k_m(*x)])
    s.addNode((0,1),[lambda *x: h*k_m(*x)])
    

    idstencil = Stencil(nsd=2,nodes={(0,0): 1.0})
    A = StencilList(g)
    A.addStencil(s,g.innerPoints())

    B = StencilSet(g)
    B.addStencil(idstencil,g.boundary())

    S = StencilSet(g)
    S.addStencil(idstencil,g.innerPoints())

    # Go parallel if we are on parallel stuff
    g.partition(A)
    B.doInitParallel()
    S.doInitParallel()

    u = Field(g)
    u.fill_vec(initf)
    d = Numeric.array(u.data)

    # Source and boundary fields
    F = Field(g)
    dF = Numeric.array(F.data)
    bF = Field(g)
    dbF = Numeric.array(bF.data)

    # build data structures in StencilSets
    A.buildMatrixOperator(u)
    B.buildMatrixOperator(u)
    S.buildMatrixOperator(u)

    # reset to start-time again
    rt = dt

    if myrank <= 0:
        print "start..."
    ct = -time.clock()
    laps = -time.time()

    while rt < T:
        # update source field
        F.fill_vec(dt_f)
        bF.fill_vec(bf)

        dF[:] = F.data
        dbF[:] = bF.data

        d[:] = A.direct_matvec(d) + S.direct_matvec(dF) + B.direct_matvec(dbF)

        # this should really go away... (Numeric vs. numarray problem):
        if u.isParallel:
            u.data[:] = d
            u.updateField()
            d[:] = u.data[:]

        rt += dt

    ct += time.clock()
    laps += time.time()
    if myrank <= 0:
        print "stop..."
        print "Used cputime in comp.loop: ",ct
        print "Elapsed time in comp.loop: ",laps

    if not u.isParallel:
        u.data[:] = d
    return u,g,rt

if (__name__ == '__main__') : runstandalone()
