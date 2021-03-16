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


# We want to solve u_tt = c^2\nabla^2 u in \Omega
#              
#                  u(x,0) = F(x)
#                  u_t(x,0) = 0.0
#                  
#                  and some various boundary conditions.
#
# options to the solver: 
#  -t : end time
#  -n : number of time steps
#  -m : division in x-direction and y-direction

import getopt, sys
from pyFDM import *
from math import exp,sin,cos,pi,sqrt

def runstandalone():
    # some defaults.
    T=1.0
    m=8
    dt = 1.0/m
    n=T/dt

    parasize = pypar.size()
    myrank = pypar.rank()

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
        import Numeric
        allopts = Numeric.zeros(4,typecode='d')
        if myrank == 0:
            allopts = Numeric.array((T,n,m))
        pypar.broadcast(allopts,0)
        if myrank != 0:
            T = float(allopts[0])
            n = int(allopts[1])
            m = int(allopts[2])

    solution,g = solve(T,n,m)
    
    if solution.isParallel:
        pypar.finalize()

def usage(T,n,m):
    print "-t T: solve from 0 to T (default %s)" % (T)
    print "-n steps: Split time in so many steps (default %s)" % (n)
    print "-m div: Split x/y-direction in this many nods (default %s)" % (m)


def F(x,y):
    if x>=0.375 and x <= 0.625:
        return 1.0
    else:
        return 0.0

def solve(T, n, m):
    print "Solve the problem from 0 to %s in %s steps, %s x %s grid-nodes" % (T,n,m,m)

    g = Grid(domain=([0,1],[0,1]), div=(m,m))
   
    #some constants
    dt = T/(1.0*n)
    # ensure dt according to CFL
    dt = min(dt,g.dx)
    print "Use dt: ",dt

    hx = (dt**2)/(g.dx**2)
    hy = (dt**2)/(g.dy**2)

    # All the stencils:
    s_pre = Stencil(nsd=2, varcoeff=False, nodes={(0,0): -1.0})
    s_lap = Stencil(nsd=2, varcoeff=False,\
            nodes={(0,0): 2*(1.0-hx-hy), (-1,0): hx, (1,0): hx, (0,-1): hy, (0,1): hy})
    
    # Dirichlet boundary based on bc1_fu function
#    bc1 = DirichletBoundary(2,bc1_fu)
#
#    # Dirichlet boundary, constant value:
    bc_diri = DirichletBoundary(2,0.0)


    # StencilSets, 
    # A is the constant -1.0 (for -1.0*u^{n-1})
    A = StencilSet(g)
    A.addStencil(s_pre,g.allPoints())

    # B holds the main stencils (for u^n).
    B = StencilSet(g)
    B.addStencil(s_lap,g.innerPoints())

#    bc_neum_l = createNeumanBoundary(s_lap,g,0.0,region=([0+g.dx/2.,1-g.dx/2.],[0,0]))
#    bc_neum_r = createNeumanBoundary(s_lap,g,0.0,region=([0+g.dx/2.,1-g.dx/2.],[1,1]))
#
#    B.addStencil(bc_diri, g.boundary(region=([0,0],[0,1])))
#    B.addStencil(bc_diri, g.boundary(region=([1,1],[0,1])))
#    B += bc_neum_l
#    B += bc_neum_r

    B += createNeumanBoundary(s_lap,g,0.0)

    uprev = Field(g)

    # Set initial condition. uprev should be filled with the function F,
    # u should then be -dt*A(G) + 0.5*B(F), where A(G) is A applied to a field 
    # filled with the function G, and B(F) is B applied to a field filled
    # with F (which is uprev...)
    uprev.fill(F)
    u = 0.5*B(uprev)
    
    rt = 0.0
    contflag=True
    while rt < T:
        unew = B(u) + A(uprev)

        # shift
        uprev = u
        u = unew
        #print u.shapedata
        rt += dt
        plot(field=u,movie='on',title='Wave propagation, t = %s' % rt)
        if contflag:
            print ">>return to continue - g to go on<<"
            input = sys.stdin.readline()[:-1]
            if input == "g":
                contflag = False

    return (u,g)

if (__name__ == '__main__') : runstandalone()
