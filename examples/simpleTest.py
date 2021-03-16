#!/usr/bin/env python

from tools import *

grid = Grid((n),([l,L],))

lap = Laplace(grid)
id = Identity(grid.nsd)
inners = id + dt*lap

# add a source - function!
# What if I want time in the source-function?
#
# init time:
tc = t
# source 
sourcefu = lambda *x: 2.0*x[0]*T
# use time in the function:
#sourcefu = lambda *x: 2.0*x[0]*tc
inners.convert_to_varcoeff()
inners.addSource(sourcefu)


sl = StencilList(grid)
innerind = sl.addStencil(inners,grid.innerPoints())

#bs = DirichletBoundary(grid.nsd,bc)
#bsind = sl.addStencil(bs,grid.boundary())

def neumannc(*x):
    if x[0] < 0.5:
        return 4*pi
    else:
        return 3.5*pi

sl += NeumanBoundary.create(sl[innerind],grid.range(),grid,neumannc)

# partition if more than 1 cpu
grid.partition(sl)


u = Field(grid)
u.fill(initc)

#sl.buildMatrixOperator(u)

while tc < T:
    u = sl(u)
    print u.myrank,"> ",u.norm()
    #u.data = dot(sl.A,u.data)
    #u.data += sl.f
    tc += dt

    # In this case we have a time-dep. source:
    #sl.buildMatrixOperator_updateSource(u)
    #sl.updateSourceDataStructures()
    #u.plot(movie='on',title='Heat transfer in 1D, t = %e' % (tc))

#u.plot(title='Finish, Heat transfer in 1D')

if u.isParallel:
    pypar.finalize()
