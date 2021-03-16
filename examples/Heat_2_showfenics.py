#!/usr/bin/env python

import getopt, sys
from pyFDM import *
from math import exp,sin,cos,pi,sqrt
import time

def runstandalone():
    T=0.0025
    n=20
    m=40

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
    analytical = lambda *x: exp(-T)*sin(pi*x[0])*cos(pi*x[1])

    u_a = Field(g)
    u_a.fill(analytical)

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
    # defaults: solve the problem from 0 to T, dt = 1/1000, dx = 1/100
    if pyFDMconf['parallel_rank'] == 0:
        print "Solve the problem from 0 to %s in %s steps, %s grid-nodes" % (T,n,m)

    # Running time:
    rt = 0.0

    # compute dt:
    dt = T/(1.0*n)

    # Create a grid
    g = Grid((m,m),([0.0,1.0],[0.0,1.0]))

    dx = g.dx
    h = dt/(dx*dx)

    # Functions I need (use rt, dt, and dx)
    dt_f = lambda *x: dt*(-exp(-rt)*sin(pi*x[0])*cos(pi*x[1]) - pi*exp(-rt)*cos(pi*x[0])*cos(pi*x[1]) + pi*exp(-rt)*sin(pi*x[0])*sin(pi*x[1]) + 2*(x[0] + x[1])*pi*pi*exp(-rt)*sin(pi*x[0])*cos(pi*x[1]))
    # need only two k-functions, as dx is the same in both directions!
    k_p = lambda *x: x[0] + x[1] + 0.5*dx
    k_m = lambda *x: x[0] + x[1] - 0.5*dx
    # boundary fu:
    bf = lambda *x: exp(-rt)*sin(pi*x[0])*cos(pi*x[1])
    # initical condition:
    initf = lambda *x: sin(pi*x[0])*cos(pi*x[1])


    # build the stencil
    s = Stencil(nsd=2,varcoeff=True,source=dt_f)
    s.addNode((0,-1),[lambda *x: h*k_p(*x)])
    s.addNode((-1,0),[lambda *x: h*k_p(*x)])
    s.addNode((0,0),[lambda *x: 1.0 - 2.0*h*(k_m(*x) + k_p(*x))])
    s.addNode((1,0),[lambda *x: h*k_m(*x)])
    s.addNode((0,1),[lambda *x: h*k_m(*x)])
    
    bs = DirichletBoundary(2,bf) 

    A = StencilList(g)
    A.addStencil(s,g.innerPoints())
    A.addStencil(bs,g.boundary())

    # Go parallel if we are on parallel stuff
    g.partition(A)

    u = Field(g)
    u.fill(initf)

    # reset to start-time again
    rt = dt 

    #import Numeric
    u_data = Numeric.array(u.data)
    u_del = A(u)

    if pyFDMconf['parallel_rank'] == 0:
        print "start..."
    ct = -time.clock()
    laps = -time.time()

    while rt < T:
        #u = A(u)
        u_data[:] = A.direct_matvec(u_data)
        rt += dt

        # This is only necessary when u uses numarray
        # should go away in the future...
        if u.isParallel:
            u.data[:] = u_data
            u.updateField()
            u_data[:] = u.data[:]

        #A.updateSourceDataStructures()
        A.buildMatrixOperator_updateSource(u)


    ct += time.clock()
    laps += time.time()
    u.data[:] = u_data[:]
    if pyFDMconf['parallel_rank'] == 0:
        print "stop..."
        print "Used cputime in comp.loop: ",ct
        print "Elapsed time in comp.loop: ",laps

    return u,g,rt

if (__name__ == '__main__') : runstandalone()
