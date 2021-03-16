#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# (c) Copyright 2003, 2004, 2005
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
#
# Heat; A Heat simulator, utilizing the FDM Stencil framework
#

# this example uses old-style import, see Wave2D_convanalysis.py for a better example.
import sys

from pyFDM import StencilList as _sl 
from pyFDM import Grid as _g
from pyFDM import Field as _f
from pyFDM import ImplementedStencils as _is
from pyFDM import Utils as _u
from pyFDM.Utils import _n

from pyFDM import NeumanBoundary
from copy import copy

from pyFDM import *

HeatError = "Error in Heat"

class Heat:
    """ Members:
         description: Just a description string, may be empty
         grid:          The grid we want to simulate on
         g:              A function, the boundary condition
         f:              A function, the initial condition

         Default Boundary type is dirichlet. Use setBoundaryType to switch.
    """
    def __init__(self, descr="empty", tstart=0.0, tstop=0.02, nsteps=2, implicit=0, strategy='ConjGrad'):
        self.description = descr
        self.grid = None
        self.u = None
        self.g = None
        self.f = None
        self.setTime(tstart,tstop,nsteps)

        self.dirichlet = 1
        self.neumann = 2
        self.boundaryType = self.dirichlet
        self.boundarycond = None

        # control explicit vs. implicit schema.
        if implicit:
            self.implicit = implicit
        else:
            self.implicit = 0

        self.strategy = strategy

        self.plotintimeloop = 0
        # Use vectorized filling? Default = No
        self.vecfill = False
                     
    def setGrid(self,grid):
        self.grid = grid
        self.u = _f.Field(self.grid)

    def setBoundaryCondition(self,g):
        """Set the function (pointer) g as the boundary condition for the Heat
            equation"""
        self.boundarycond = g

    def setBoundaryType(self,type):
        """Set the type of boundary-condition to use"""

        self.boundaryType = type

    def setInitialCondition(self,f):
        """Set the function (pointer) f as the initial condition for the Heat
            equation"""
        self.initialcond = f

    def setImplicit(self, implicit=1):
        """Set the self.implicit variable, to control whether implicit or explicit 
           method is used to solve the problem
        """
        if implicit:
            self.implicit = implicit
        else:
            self.implicit = 0

    def setImplicitStrategy(self, strategy='ConjGrad'):
        """Set implicit solution strategy"""

        self.strategy = strategy


    def setTime(self,tstart,tstop,nsteps=10):
        self.t0 = tstart
        self.ts = tstop
        self.nsteps = nsteps
        self.dt = (self.ts - self.t0) / nsteps
        self.curtime = self.t0


    def buildStencilCollection(self, implicit=0):
        """Build the Stencil Collection for the heat equation.
           If the keyword implicit is set to 1 or self.implicit is set to 1 before 
           buildStencilCollection is called, a colletion for using this as an 
           implicit operator in build, else for explicit schema. That is:
             Explicit: I + dt*Laplace
             Implicit: I - dt*Laplace

             Boundary-conditions are the same for explicit and implicit
        """
        
        if implicit:
            self.implicit = implicit

        if not self.grid:
            raise HeatError, "Heat is not properly initialized, grid not set"

        if self.boundarycond == None:
            raise HeatError, "Boundary condition functions is not set"

        self.stencillist = _sl.StencilList(self.grid)
        nsd = self.grid.nsd

        laplace = _is.Laplace(self.grid)
        idStencil = _is.Identity(nsd)
        if self.implicit:
            innerStencil = idStencil - self.dt*laplace
        else:
            innerStencil = idStencil + self.dt*laplace

        #Add inner stencil
        self.innerindex = self.stencillist.addStencil(innerStencil,self.grid.innerPoints())

        if self.boundaryType == self.dirichlet:
            self.setDirichletBoundary()
        elif self.boundaryType == self.neumann:
            self.setNeumannBoundary()
        else:
            raise HeatError, "Unknown boundaryTyped assigned"


    def setDirichletBoundary(self):
        """Proivde a Dirichlet boundary for the problem. 
            The boundary condition must be set as a method upfront.
        """

        print "Using a dirichlet boundary"

        boundaryStencil = _is.DirichletBoundary(self.grid.nsd,self.boundarycond)
        index = self.stencillist.addStencil(boundaryStencil,self.grid.boundary())

    def setNeumannBoundary(self):
        """Provide a Neumann boundary for the problem
            The boundary condition must be set as a method upfront.
        """

        # The create method in NeumanBoundary return a new stencillist, so 
        # we need to add that one togheter with the one we already have here.
        self.stencillist += NeumanBoundary.create(self.stencillist[self.innerindex],self.grid.range(),self.grid,self.boundarycond)

    def initField(self):
        """Fill the solution-field with the initial condition:"""
        if self.vecfill:
            self.u.fill_vec(self.initialcond)
        else:
            self.u.fill(self.initialcond)

    def plot(self):
        """Wrapper for plotting the field in the heat solver"""
        self.u.plot()

    def solve(self, tolerance=1.0E-5):
        """solve just run the timeloop, included for convenience
           tolerance is for implicit solver, passed through to timeloop.
        """
        self.timeloop(tolerance)

    def plot_in_timeloop(self,flag):
        if flag==1:
            self.plotintimeloop = 1
        else:
            self.plotintimeloop = 0

    def set_vecfill(self,flag):
        """set whether vectorized fill should be used or not
           Parameter: True/False
        """
        self.vecfill = flag

    def timeloop(self, tolerance=1.0E-5):
        """timeloop: 
           Solve the problem in the given timerange.  
            - If explicit solver is used, just a simple iteration applying the 
              stencilcollection at each timestep 
            - If an implicit solver is asked for, a "double-loop" is performed.
           tolerance is only used in the implicit case.
        """
        if not self.stencillist:
            raise HeatError, "Stencillist is not initialized"
        while self.curtime < self.ts:
            #print "time is now %e" % ( self.curtime )
            if self.implicit:
                """solve using an implicit scheme"""
                
                # TODO: create a function pointer to the implicit solver instead, 
                # that will avoid the if-test here.
                if self.strategy == "Jacobi":
                    self.u = jacobi(self.stencillist, self.u, self.u, tolerance, True)
                elif self.strategy == "ConjGrad":
                    self.u = conjgrad(self.stencillist, self.u, self.u, tolerance, True)
            else:
                self.u = self.stencillist(self.u)
            self.curtime += self.dt
            if self.plotintimeloop == 1:
                self.u.plot(movie='on',title='Heat transfer')
            
def neumanBoundary1(x):
    if x[1] < 0.1:
        return -0.5
    elif x[1] > 0.9: 
        return 0.5
    else:
        return 0.0

def runstandalone():
    try: 
        t_end = float(sys.argv[1])
        t_steps = int(sys.argv[2])
    except:
        t_end = 0.0001
        t_steps = 10

    try: 
        partition_n = int(sys.argv[3])
    except: 
        try:
            partition_n = int(sys.argv[1])
        except:
            partition_n = 10

    try: 
        implicit = int(sys.argv[4])
    except:
        implicit = 1
    

    heat = Heat("Heat solver",0.0,t_end,t_steps)
    grid = _g.Grid(div=(partition_n,partition_n),domain=([0.0,1.0],[0.0,1.0]))
    #grid = _g.Grid((40,40),2,([0.0,1.0],[0.0,1.0]))
    #grid = _g.Grid((50,50),2,([0.0,1.0],[0.0,1.0]))
    #grid = _g.Grid((40,40),2,([0.0,1.0],[0.0,1.0]))
    #grid = _g.Grid((40,),1,([0.0,1.0],))



    heat.setGrid(grid)
    heat.setInitialCondition(lambda *x: 2*_n.sin(2.0*_n.pi*x[0])*_n.sin(2.0*_n.pi*x[1]))
    #heat.setInitialCondition(lambda x: _n.sin(2.0*_n.pi*x[0]))

    heat.setBoundaryCondition(lambda *x: (x[1]==0 and x[0]*x[0] or 0.))
    #heat.setBoundaryCondition(neumanBoundary1)
    #heat.setBoundaryCondition(1.0)

    heat.setBoundaryType(heat.neumann)
    #heat.setBoundaryType(heat.dirichlet)
    heat.setImplicit(implicit)
    heat.setImplicitStrategy('ConjGrad')
    heat.buildStencilCollection();

    # go parallel:
    grid.partition(heat.stencillist)

    #heat.initField()
    heat.set_vecfill(False)
    heat.initField()

#    allind = _u.tupleIterator(2,(0,0),(partition_n,partition_n))
#    for i in allind:
#        print "%s: %e" % (i,heat.u[i])

    # turn on plotting
    heat.plot_in_timeloop(1)
    heat.solve(1.0E-3)


    # no final plot either, have to deal with this differently on chilo!
    #heat.plot()
#    allind.reset()
#    for i in allind:
#        print "%s: %e" % (i,heat.u[i])

    #print heat.u.data

    # should find some other way to handle this!
    if heat.u.isParallel:
    	_u.pypar.finalize()

    return heat

if (__name__ == '__main__') : runstandalone()                     
