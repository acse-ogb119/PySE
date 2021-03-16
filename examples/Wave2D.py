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


# We want to solve u_tt = \nable \cdot (f(x) \nabla u) in \Omege_1
#              and u_tt = c^2\nabla^2 u in \Omega_2
#              
#                  u(x,0) = F(x)
#                  u_t(x,0) = G(x)
#                  
#                  and some various boundary conditions.
#
# options to the solver: 
#  -t : end time
#  -n : number of time steps
#  -mx : division in x-direction
#  -my : division in y-direction

import getopt, sys
from pyFDM import *
from math import exp,sin,cos,pi,sqrt

def runstandalone():
    T=0.01
    n=10
    mx=40
    my=80

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
            if opt == "-mx":
                mx = int(opts[i+1])
                i += 1
            if opt == "-my":
                my = int(opts[i+1])
                i += 1
            i += 1
    
    if parasize > 1:
        import Numeric
        allopts = Numeric.zeros(4,typecode='d')
        if myrank == 0:
            allopts = Numeric.array((T,n,mx,my))
        pypar.broadcast(allopts,0)
        if myrank != 0:
            T = float(allopts[0])
            n = int(allopts[1])
            mx = int(allopts[2])
            my = int(allopts[3])

    solution,g = solve(T,n,mx,my)
    
    if solution.isParallel:
        pypar.finalize()

def usage(T,n,mx,my):
    print "-t T: solve from 0 to T (default %s)" % (T)
    print "-n steps: Split time in so many steps (default %s)" % (n)
    print "-mx div: Split x-direction in this many nods (default %s)" % (mx)
    print "-my div: Split y-direction in this many nods (default %s)" % (my)


def F(x,y):
    Au = 0.4
    s_ux = 0.05
    s_uy = 0.1
    xc_u = 0.5
    yc_u = 1.0
    #return Au*exp(-0.5*((x - xc_u)/s_ux)**2 - 0.5*((y - yc_u)/s_uy)**2)
    return 0.0

def G(x,y):
    return 0.0

def H(x,y):
    Au = 0.5
    s_ux = 0.1
    s_uy = 0.1
    xc_u = 0.7
    yc_u = 1.3
    return 1 - Au*exp(-0.5*((x - xc_u)/s_ux)**2 - 0.5*((y - yc_u)/s_uy)**2)

#def H(x,y):
#    a = 0.4
#    b = 6
#    d = -0.01
#    k = 3
#    return x + a*sin(b*pi*(x-1.0)) 
#    #return x + a*sin(b*pi*(x-1.0)) + d*(x-1.0)*sin(k*pi*y)
#    #return 1.0

def bc1_fu(x,y,t):
    p = 0.5
    val = -p*(x*x - x)*sin(pi*x*2.)*sin(4.*pi*t)
    return  val

def nc1_fu(x,y):
    return 0.0

def nc2_fu(x,y):
    return 0.05

def nc3_fu(x,y):
    return -0.1

def solve(T, n, mx, my):
    print "Solve the problem from 0 to %s in %s steps, %s x %s grid-nodes" % (T,n,mx,my)

    g = Grid(domain=([0,1],[0,2]), div=(mx,my))
   
    #some constants
    dt = T/(1.0*n)

    # ensure dt according to CFL
    dt = min(dt,sqrt(1.0/1.4)*sqrt(1.0/((1.0/g.dx**2) + (1.0/g.dy**2))))
    rt = 0.0
    print "Use dt: ",dt

    hx = (dt**2)/(g.dx**2)
    hy = (dt**2)/(g.dy**2)

    # A few constants for stencils and region-defs (epsilons of 1/2 gridcells)
    ex = 0.5*g.dx
    ey = 0.5*g.dy

    # variants of H-call:
    Hxm = lambda x,y: hx*H(x-ex,y)
    Hxp = lambda x,y: hx*H(x+ey,y)
    Hym = lambda x,y: hy*H(x,y-ey)
    Hyp = lambda x,y: hy*H(x,y+ey)
    center = lambda x,y: 2.0 - (Hxp(x,y) + Hxm(x,y)) - (Hyp(x,y) + Hym(x,y))

    # All the stencils:
    s1 = Stencil(nsd=2, varcoeff=False, nodes={(0,0): -1.0})

    s2 = Stencil(nsd=2, varcoeff=True,\
            nodes={(0,0): center, (-1,0): Hxp, (1,0): Hxm, (0,-1): Hyp, (0,1): Hym})
    s3 = Stencil(nsd=2, varcoeff=False,\
            nodes={(0,0): 2*(1.0-hx-hy), (-1,0): hx, (1,0): hx, (0,-1): hy, (0,1): hy})

    
    # Dirichlet boundary, constant value:

    bc1_fu_call = lambda x,y: bc1_fu(x,y,rt)
    bc1 = DirichletBoundary(2, bc1_fu_call)
    bc3 = DirichletBoundary(2,0.05)
    bc4 = DirichletBoundary(2,-0.05)


    # StencilSets, 
    # A is the constant -1.0 (for -1.0*u^{n-1})
    A = StencilSet(g)
    A.addStencil(s1,g.allPoints())

    # B holds the main stencils (for u^n).
    B = StencilSet(g)
    B.addStencil(s2,g.innerPoints(region=([0.4+ex,1],[1+ey,2])))
    B.addStencil(s2,g.innerPoints(region=([0,0.4],[1+ey,1.6])))
    B.addStencil(s2,g.innerPoints(region=([0,0.1],[1.6+ey,2])))
    B.addStencil(s2,g.innerPoints(region=([0.1+ex,0.4],[1.9,2])))

    B.addStencil(s3,g.innerPoints(region=([0,0.5],[0,1])))
    B.addStencil(s3,g.innerPoints(region=([0.5+ex,1],[0,0.2-ey])))
    B.addStencil(s3,g.innerPoints(region=([0.5+ex,1],[0.7+ey,1])))
    B.addStencil(s3,g.innerPoints(region=([0.75,1],[0.2,0.7])))

#    B.addStencil(s2,g.innerPoints())

    # Neuman Conditions
#    # for region with s3 as main stencil
#    nc1 = createNeumanBoundary(s3,g,nc1_fu,region=([0.5+ex,1],[0,1]))
#    nc2 = createNeumanBoundary(s3,g,nc1_fu,region=([0,0],[0,1]))
#    nc5 = createNeumanBoundary(s3,g,0.0,region=([0,0.5],[0,0]))
#    # for region with s2 as main stencil
#    nc3 = createNeumanBoundary(s2,g,nc1_fu,region=([0,1-ex],[1+ey,2]))
#    nc4 = createNeumanBoundary(s2,g,nc2_fu,region=([1,1],[1+ey,2]))
#
##    nczero1 = createNeumanBoundary(s2,g,nc1_fu,region=([0,1],[1+ey,2]))
##    nczero2 = createNeumanBoundary(s3,g,nc1_fu,region=([0,1],[0,1]))
    nczero1 = createNeumanBoundary(s2,g,0.0,region=([0,1],[1+ey,2]))
    #nczero4 = createNeumanBoundary(s2,g,0.0,region=([1,1],[1+ey,1.5]))

    nczero2 = createNeumanBoundary(s3,g,0.0,region=([0.5+ex,1],[0,1]))
    nczero3 = createNeumanBoundary(s3,g,0.0,region=([0,0],[0+ey,1]))
#    nczero2 = createNeumanBoundary(s2,g,0.0,region=([0,1],[0,2]))


    # Add boundary-conditions in B (Dirichlet can be added anywhere, but 
    # Neuman conditions must be in B. And for use to set initial-condition, all 
    # boundary conditions must be in B).
    # Dirichlet boundaries:
    # Neyman boundaries:
#    B += nc1
#    B += nc2
#    B += nc3
#    B += nc4
#    B += nc5
    B += nczero1
    B += nczero2
    B += nczero3
    #B += nczero4
#    B.addStencil(bc2, g.boundary(region=([0.5+ex,1],[0,1.5])))
#    B.addStencil(bc2, g.boundary(region=([0,1],[ey,2])))
    B.addStencil(bc1, g.boundary(region=([0,0.5],[0,0])))
    #B.addStencil(bc4, g.boundary(region=([1,1],[1.5+ey,2])))

    uprev = Field(g)
    # no need to create u, it will be created automatically...
    #u = Field(g) 
    Gfield = Field(g) 

    # Set initial condition. uprev should be filled with the function F,
    # u should then be -dt*A(G) + 0.5*B(F), where A(G) is A applied to a field 
    # filled with the function G, and B(F) is B applied to a field filled
    # with F (which is uprev...)
    uprev.fill(F)
    Gfield.fill(G)
    u = -dt*A(Gfield) + 0.5*B(uprev)
    
    #plot(field=u, title='Initical condition')
    #contflag=True
    while rt < T:
        unew = B(u) + A(uprev)

        # shift
        uprev = u
        u = unew
        #print u.shapedata
        rt += dt
        plot(field=u,movie='on',title='Wave propagation, t = %s' % rt)
        B.updateSourceDataStructures()
#        if contflag:
#            print ">>return to continue - g to go on<<"
#            input = sys.stdin.readline()[:-1]
#            if input == "g":
#                contflag = False

    return (u,g)

if (__name__ == '__main__') : runstandalone()
