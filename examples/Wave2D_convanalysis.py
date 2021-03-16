#!/usr/bin/env python
# -*- coding: iso8859-1 -*-
#
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
from math import exp,sin,cos,pi,sqrt,ceil

from pyFDM import *

def runstandalone():
    # some defaults.
    T=1./(sqrt(2.)*2.)
    m=50
    dt = 1./m
    n=ceil(T/dt)
    dt = float(T)/n
    s = 9
    variant = 1
    solutionalg = 1
    tolerance = 1.0e-05
    halfwave = False
    waves = 1.0

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
            if opt == "-e":
                tolerance = float(opts[i+1])
                i += 1
            if opt == "-n":
                n = int(opts[i+1])
                i += 1
            if opt == "-m":
                m = int(opts[i+1])
                i += 1
            if opt == "-s":
                s = int(opts[i+1])
                i += 1
            if opt == "-a":
                solutionalg = int(opts[i+1])
                i += 1
            if opt == "-v":
                variant = int(opts[i+1])
                i += 1
            if opt == "-w":
                waves = float(opts[i+1])
                i += 1
                if waves == 0.5:
                    halfwave = True
            if opt == "-h":
                usage(T,n,m,s,variant,solutionalg,tolerance)
                sys.exit(0)
            i += 1
    
    if parasize > 1:
        import Numeric
        allopts = Numeric.zeros(9,typecode='d')
        if myrank == 0:
            allopts = Numeric.array((T,n,m,s,variant,solutionalg,int(halfwave),tolerance,waves))
        pypar.broadcast(allopts,0)
        if myrank != 0:
            T = float(allopts[0])
            n = int(allopts[1])
            m = int(allopts[2])
            s = int(allopts[3])
            variant = int(allopts[4])
            solutionalg = int(allopts[5])
            halfwave = bool(int(allopts[6]))
            tolerance = float(allopts[7])
            waves = float(allopts[8])

    if s == 5:
        scheme="5pt"
    elif s == 9:
        scheme="9pt"
    else:
        raise "WaveError","Unknown scheme: %s" % (s)

    if variant == 1:
        if halfwave:
            AnSol = AnSol_I_half
        else:
            AnSol = lambda x,y,t: AnSol_I(x,y,t,waves)
    elif variant == 2:
        if halfwave:
            AnSol = AnSol_II_half
        else:
            AnSol = lambda x,y,t: AnSol_II(x,y,t,waves)
    else:
        raise "WaveError","Unknown variant: %s" % (variant)

    solution,g,l2norm,T = solve(T,n,m,scheme,variant,solutionalg,halfwave,tolerance=tolerance,waves=waves)
    
    analytical = Field(g)
    analytical_fu = lambda x,y,t=T: AnSol(x,y,t)
    analytical.fill(analytical_fu)

    error = analytical - solution
    #plot(field=error,title="Error")
    #plot(field=analytical,title="Analytical solution at t=%s" % (T))
    plot(field=solution,title="Approx. solution at t=%s" % (T))

    errornorm = error.normL2()
    if myrank == 0:
        print "The norm of the error now is: ",errornorm
        print "The error in time and space is: ",l2norm
        print "g.dx = ",g.dx
#
#
#    rt = 0.0
#    dt = T/(1.*n)
#    while rt < T:
#        analytical_fu = lambda x,y,t=rt: AnSol(x,y,t)
#        analytical.fill(analytical_fu)
#        plot(field=analytical,movie='on',title='analytical solution, t = %s' %(rt))
#        rt += dt
#

    if solution.isParallel:
        pypar.finalize()

def usage(T,n,m,s,variant,solutionalg,tolerance):
    print "-t T: solve from 0 to T (default %s)" % (T)
    print "-n steps: Split time in so many steps (default %s)" % (n)
    print "-m div: Split x/y-direction in this many nods (default %s)" % (m)
    print "-s scheme: Use 5 or 9 point scheme, specify s as 5 or 9 (default %s)" % (s)
    print "-v variant: 1 or 2. 1 use neumann boundary, 2 use dirichelet (default %s)." % (variant)
    print "-a solutionalgorithm: 1, 2, or 3: 1 use centered diff in time, 2 use Runge Kutta of 4. order, 3 use a modified laplace (default %s)" % (solutionalg)
    print "-w choose another init function and analytical solution, half wave"
    print "-e tol: tolerance for jacobi-solver in sol.alg. 3. (defaults %s)" % (tolerance)


def F_I(x,y,a=1.0):
    return cos(a*pi*x)*cos(a*pi*y) 

def F_I_half(x,y):
    return cos(0.5*pi + 0.5*pi*x)*cos(0.5*pi + 0.5*pi*y) 

def F_II(x,y,a=1.0):
    return sin(a*pi*x)*sin(a*pi*y) 

def F_II_half(x,y):
    return sin(0.5*pi+0.5*x*pi)*sin(0.5*pi+0.5*y*pi)

def AnSol_I(x,y,t,a=1.0):
    return cos(a*pi*sqrt(2.0)*t)*cos(a*pi*x)*cos(a*pi*y)

def AnSol_I_half(x,y,t):
    return cos(0.5*pi*sqrt(2.0)*t)*cos(0.5*pi + 0.5*pi*x)*cos(0.5*pi + 0.5*pi*y)

def AnSol_II(x,y,t,a=1.0):
    return cos(a*pi*sqrt(2.0)*t)*sin(a*pi*x)*sin(a*pi*y)

def AnSol_II_half(x,y,t):
    return cos(0.5*pi*sqrt(2.0)*t)*sin(0.5*pi + 0.5*pi*x)*sin(0.5*pi + 0.5*pi*y)

def solve(T, n, m, scheme="9pt",variant=1,solutionalg=1,halfwave=False,tolerance=1.0e-07,waves=1.0):
    myrank = pypar.rank()
    if myrank == 0:
        print "Solve the problem from 0 to %s in %s steps, %s x %s grid-nodes" % (T,n,m,m)

    g = Grid(domain=([-1,1],[-1,1]), div=(m,m))
   
    #some constants
    dt = T/(1.*n)
    # ensure dt according to CFL
    dt = min(dt,g.dx)
    if myrank == 0:
        print "Use dt: ",dt

    h = (dt**2)/(g.dx**2)

    # All the stencils:
    s_prepre = Stencil(nsd=2, varcoeff=False, nodes={(0,0): -1.0})
    s_pre = Stencil(nsd=2,varcoeff=False, nodes={(0,0): 2.0})

    s_lap_5pt = Stencil(nsd=2, varcoeff=False,\
            nodes={            (0,1):  1.,\
                   (-1,0): 1., (0,0):  -4., (1,0): 1.,\
                               (0,-1): 1.})

    s_lap_9pt = Stencil(nsd=2, varcoeff=False,\
            nodes={(-1,1):  1./6, (0,1):  2./3,   (1,1):  1./6, \
                   (-1,0):  2./3, (0,0):  -10./3, (1,0):  2./3, \
                   (-1,-1): 1./6, (0,-1): 2./3,   (1,-1): 1./6})
   
    s_mass = Stencil(nsd=2,varcoeff=False,\
            nodes={               (0,1): 1./12,\
                   (-1,0): 1./12, (0,0): 2./3,  (1,0): 1./12,\
                                  (0,-1): 1./12})

    # StencilSets, 
    # A is the constant -1.0 (for -1.0*u^{n-1})
    A = StencilSet(g)
    A.addStencil(s_prepre,g.allPoints())

    # B holds the main stencils (for u^n).
    B = StencilSet(g)
    # L holds the laplace operator for Runge Kutta based solutionalg.
    L = StencilSet(g)
    if scheme == "5pt":
        s_inner = s_pre + h*s_lap_5pt
        s_lap = (1./g.dx**2)*s_lap_5pt
    elif scheme == "9pt":
        s_inner = s_pre + h*s_lap_9pt
        s_lap = (1./g.dx**2)*s_lap_9pt
    else:
        raise "WaveError", "Unknown scheme specified"

    B.addStencil(s_inner,g.innerPoints())
    L.addStencil(s_lap,g.innerPoints())
    
    M = StencilSet(g)
    M.addStencil(s_mass,g.innerPoints())

#    bc_neum_l = createNeumanBoundary(s_lap,g,0.0,region=([0+g.dx/2.,1-g.dx/2.],[0,0]))
#    bc_neum_r = createNeumanBoundary(s_lap,g,0.0,region=([0+g.dx/2.,1-g.dx/2.],[1,1]))
#
#    B.addStencil(bc_diri, g.boundary(region=([0,0],[0,1])))
#    B.addStencil(bc_diri, g.boundary(region=([1,1],[0,1])))
#    B += bc_neum_l
#    B += bc_neum_r

   

    if variant == 1:
        B += createNeumanBoundary(s_inner,g,0.0)
        L += createNeumanBoundary(s_lap,g,0.0)
        M += createNeumanBoundary(s_mass,g,0.0)
        if halfwave:
            AnSol = AnSol_I_half
            F = F_I_half
        else:
            AnSol = lambda x,y,t: AnSol_I(x,y,t,waves)
            F = lambda x,y: F_I(x,y,waves)
    elif variant == 2:
        bc_diri = DirichletBoundary(2,0.0)
        B.addStencil(bc_diri,g.boundary())
        L.addStencil(bc_diri,g.boundary())
        M.addStencil(bc_diri,g.boundary())
        if halfwave:
            AnSol = AnSol_II_half
            F = F_II_half
        else:
            AnSol = lambda x,y,t: AnSol_II(x,y,t,waves)
            F = lambda x,y: F_II(x,y,waves)
    else: 
        raise "Error","Unknown variant"

    g.partition(B)
    A.doInitParallel()
    L.doInitParallel()
    M.doInitParallel()

    uprev = Field(g)
    unew = Field(g)
    error = Field(g)

    #initial conditions:
    uprev.fill(F)
    #analytical_fu = lambda x,y,t=0.0: AnSol(x,y,t)
    #uprev.fill(analytical_fu)
    
    dprev = Numeric.array(uprev.data)

    rt = 0.0
    #contflag=True

    L2norm = 0.0

    #analytical = []
    analytical = Field(g)
    if solutionalg == 1:
        # std. explicit scheme.
        # but this is a first-order approximation.
        #u = 0.5*B(uprev)
        #instead, fill in the analytical solution at t=dt:
        u = Field(g)
        analytical_fu = lambda x,y,t=dt: AnSol(x,y,t)
        u.fill(analytical_fu)
        # Trigger building of datastructures in A
        u_del = A(uprev)
        u_del = B(u)
        d = Numeric.array(u.data)

        rt = dt
        while rt < T:
            #print "rt: ",rt
            #unew = B(u) + A(uprev)
            dnew = B.direct_matvec(d) + A.direct_matvec(dprev)
     
            # inc rt to get analytical solution at right timestep
            rt += dt
            analytical_fu = lambda x,y,t=rt: AnSol(x,y,t)
            analytical.fill(analytical_fu)

            unew.data[:] = dnew
            if unew.isParallel:
                unew.updateField()
                # copy back updated data:
                dnew[:] = unew.data
            error = analytical - unew
            L2norm += (error.normL2()**2)
            #print "local error: ",error.normL2()

            # shift
            #uprev = u
            #u = unew
            dprev[:] = d
            d[:] = dnew
            #print dnew
            #plot(unew,title='approx solution, t=%s' % (rt),movie='on')
            #print u.shapedata
            

            #plot(field=u,movie='on',title='Wave propagation, t = %s' % rt)
            #if contflag:
            #    print ">>return to continue - g to go on<<"
            #    input = sys.stdin.readline()[:-1]
            #    if input == "g":
            #        contflag = False

    elif solutionalg == 2:
        # 4.order Runge Kutta based scheme
        # A Nyström variant.

        # set ddt = 0 initially.
        ddt = Numeric.array(dprev)
        ddt[:] = 0.0
        # Trigger building of datastructures in L
        d_del = L(uprev)

        Lddt = Numeric.array(ddt)
        # if parallel, need to update k1, k2 (L used recursively)
        k1field = Field(g)
        k1 = Numeric.array(k1field.data)
        k2field = Field(g)
        k2 = Numeric.array(k2field.data)

        k3 = Numeric.array(dprev)
        dnew = Numeric.array(dprev)
        udtnew = Field(g)
        ddtnew = Numeric.array(udtnew.data)

        while rt < T:
            #print "rt: ",rt

            Lddt[:] = L.direct_matvec(ddt)
            k1[:] = L.direct_matvec(dprev)
            if k1field.isParallel:
                k1field.data[:] = k1
                k1field.updateField()
                k1[:] = k1field.data
            k2[:] = k1 + 0.5*dt*Lddt + (dt**2)*(1./8)*(L.direct_matvec(k1))
            if k2field.isParallel:
                k2field.data[:] = k2
                k2field.updateField()
                k2[:] = k2field.data
            k3[:] = k1 + 1.*dt*Lddt + (dt**2)*0.5*(L.direct_matvec(k2))

            dnew[:] = dprev + dt*ddt + (dt**2)*((1./6)*k1 + (1./3)*k2)
            ddtnew[:] = ddt + dt*((1./6)*k1 + (4./6)*k2 + (1./6)*k3)

            unew.data[:] = dnew
            udtnew.data[:] = ddtnew
            if unew.isParallel:
                unew.updateField()
                dnew[:] = unew.data
                udtnew.updateField()
                ddtnew[:] = udtnew.data

            # swap
            dprev[:] = dnew
            #print dnew
            ddt[:] = ddtnew



            #plot(unew,title='approx solution, t=%s' % (rt),movie='on')

            # inc rt to get analytical solution at right timestep
            rt += dt

            # error estimation
            analytical_fu = lambda x,y,t=rt: AnSol(x,y,t)
            analytical.fill(analytical_fu)
            error.data[:] = analytical.data - unew.data
            L2norm += (error.normL2()**2)

    elif solutionalg == 3:
        # Use the modification of source for 9-point laplace. 
        # Then we need to solve a linear system on each timestep.
        # We assume that the modification can be used with out problem also for 
        # 5-point laplace, as it is a O(dx^4) approximation of the "source"
        # Remark that the "source" in our case is u_tt.


        rt = dt

#        u = Field(g)
#        analytical_fu = lambda x,y,t=rt: AnSol(x,y,t)
#        u.fill(analytical_fu)

#        # approx. u^1: Solve Mu^1 = 0.5 * [ 2*M*u^0 + dt*L*u^0 ]
        b = 0.5*(2*(M*uprev) + (dt**2)*(L*uprev))
        u = conjgrad(M, uprev, b, tolerance=tolerance, relativeconv=True)
     

        while rt < T:
            b = 2*(M*u) - M*uprev + (dt**2)*(L*u)
            #unew = jacobi(M, u, b, tolerance=tolerance, relativeconv=True)
            unew = conjgrad(M, u, b, tolerance=tolerance, relativeconv=False)

            #plot(unew,title='approx solution, t=%s' % (rt),movie='on')

            rt += dt

            # estimate error
            analytical_fu = lambda x,y,t=rt: AnSol(x,y,t)
            analytical.fill(analytical_fu)
            error = analytical - unew
            L2norm += error.normL2()**2
            #print "local error: ",error.normL2()

            # swap fields:
            uprev.data[:] = u.data
            u.data[:] = unew.data

    elif solutionalg == 4:
        # 4.order Runge Kutta based scheme
        # A Nyström variant.
        # Modified for weighted "source" for 4.order convergence in time.
        
        # set ddt = 0 initially.
        ddt = Numeric.array(dprev)
        ddt[:] = 0.0
        udt = Field(g)
        udt.data[:] = ddt

        # Trigger building of datastructures in L and M
        d_del = L(uprev)
        d_del = M(uprev)

        Lddt = Numeric.array(ddt)
        # if parallel, need to update k1, k2 (L used recursively)
        k1field = Field(g)
        k1 = Numeric.array(k1field.data)
        k2field = Field(g)
        k2 = Numeric.array(k2field.data)

        k3 = Numeric.array(dprev)
        dnew = Numeric.array(dprev)
        udtnew = Field(g)
        ddtnew = Numeric.array(udtnew.data)

        bu = Field(g)
        budt = Field(g)

        while rt < T:
            #print "rt: ",rt

            Lddt[:] = L.direct_matvec(ddt)
            k1[:] = L.direct_matvec(dprev)
            if k1field.isParallel:
                k1field.data[:] = k1
                k1field.updateField()
                k1[:] = k1field.data
            k2[:] = k1 + 0.5*dt*Lddt + (dt**2)*(1./8)*(L.direct_matvec(k1))
            if k2field.isParallel:
                k2field.data[:] = k2
                k2field.updateField()
                k2[:] = k2field.data
            k3[:] = k1 + 1.*dt*Lddt + (dt**2)*0.5*(L.direct_matvec(k2))

            bu.data[:] = M*dprev + dt*(M*ddt) + (dt**2)*((1./6)*k1 + (1./3)*k2)
            budt.data[:] = M*ddt + dt*((1./6)*k1 + (4./6)*k2 + (1./6)*k3)
            unew = conjgrad(M, uprev, bu, tolerance=tolerance, relativeconv=True)
            udtnew = conjgrad(M, udt, budt, tolerance=tolerance, relativeconv=True)

            dnew[:] = unew.data
            ddtnew[:] = udtnew.data

            # swap
            dprev[:] = dnew
            uprev.data[:] = dnew
            #print dnew
            ddt[:] = ddtnew
            udt.data[:] = ddtnew


            #plot(unew,title='approx solution, t=%s' % (rt),movie='on')

            # inc rt to get analytical solution at right timestep
            rt += dt

            # error estimation
            analytical_fu = lambda x,y,t=rt: AnSol(x,y,t)
            analytical.fill(analytical_fu)
            error.data[:] = analytical.data - unew.data
            L2norm += (error.normL2()**2)



    L2norm = sqrt(L2norm*dt)
    return (unew,g,L2norm,rt)

if (__name__ == '__main__') : runstandalone()
