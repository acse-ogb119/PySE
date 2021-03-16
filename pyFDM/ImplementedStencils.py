# -*- coding: utf-8 -*-
#
# (c) Copyright 2003, 2005
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
# ImplementedStencils; Implementation of some often used stencils.
#

from Stencil import *

class Laplace(Stencil):
    """Laplace stencil: An extension of the plain Stencil type which implement a
        Laplace stencil with fixed coefficients in the given number of space
        dimensions"""

    def __init__(self,grid):
        Stencil.__init__(self,grid.nsd,0,0.0)
        
        self.grid = grid

        # Create a basis_index (must be converted to tuple before used in
        # stencil)
        basis_index = map(lambda x: 0,range(self.nsd)) 
        b_i = basis_index
        dx = grid.division
        for i in xrange(self.nsd):
            # copy the basic_index
            cbi = list(b_i)
            # add "left" node
            cbi[i]=-1
            ti=tuple(cbi)
            value = 1.0/(dx[i]*dx[i])
            self.addNode(ti,value)
            # add "right" node
            cbi[i]=1
            ti=tuple(cbi)
            self.addNode(ti,value)
            # add "center"
            cbi[i]=0
            ti=tuple(cbi)
            value=-2.0/(dx[i]*dx[i])
            self.addNode(ti,value)

class LaplaceJ(Stencil):
    """LaplaceJ stencil: An extension of the plain Stencil type which implement a
        Laplace stencil with fixed coefficients in the given number of space
        dimensions, and with the center node removed - for use in e.g Jacobi method.
        """

    def __init__(self,grid):
        Stencil.__init__(self,grid.nsd,0,0.0)
        
        self.grid = grid
        # compute the necessary paramater when partition in different directions differ:
        self.divparm = 0.0
        dx = grid.division
        for i in xrange(self.nsd):
            self.divparm += 2.0/(dx[i]*dx[i])
        self.divparm = 1.0 / self.divparm

        # Create a basis_index (must be converted to tuple before used in
        # stencil)
        basis_index = map(lambda x: 0,range(self.nsd)) 
        b_i = basis_index
        dx = grid.division # (isn't this already set...)
        for i in xrange(self.nsd):
            # copy the basic_index
            cbi = list(b_i)
            # add "left" node
            cbi[i]=-1
            ti=tuple(cbi)
            value = self.divparm * 1.0/(dx[i]*dx[i])
            self.addNode(ti,value)
            # add "right" node
            cbi[i]=1
            ti=tuple(cbi)
            self.addNode(ti,value)

class DirichletBoundary(Stencil):
    """Dirichlet Boundary stencil: Just the simplest stencil there is - just the
        source. We assume that the source is a function, and define this as a
        var.coeff stencil"""

    def __init__(self,nsd,source):
        """Arguments: no of space dimensions (nsd) and the source-function"""
        if isinstance(source,(float,int)):
            Stencil.__init__(self,nsd,0,source)
        else:
            Stencil.__init__(self,nsd,1,source)


class Identity(Stencil):
    """Identity stencil: Implement the stencil with 1.0 in the center
    node and nothing else"""

    def __init__(self,nsd):
        Stencil.__init__(self,nsd,0,0.0)

        basic_index = tuple(map(lambda x: 0,range(self.nsd)))
        self.addNode(basic_index,1.0)

# Next step: implement a neuman condition. Some stencil must be provided, and
# some boundaryrange in the grid. For each edge-part in the boundaryrange, the
# provided stencil must be modified according to the neumann-condition. Open
# questions: Should we provide a stencillist where the modified neumannstencils
# with the right ranges add, or should we just build a list and return so that
# it can be added to stencillists elsewhere.

# edgeIteratorGenerator must be used to dived the boundary-range into the
# seperate edgeparts. Used the basetu and edgtu information ( elem. 1 and 2 in
# the tuple return by edgeIteratorGenerator to find out whether the fixed
# elements are fixed at max or min value. Here, we can either assumed fixed in
# min if the edgetu entry is 0, but we might as well check with the grid.

#class NeumanBoundary(Stencil):
#    """Neuman Boundary Stencil: This class implement a Neuman boundary
#        condition. To create this, a stencil, some boundary and the neuman
#        condition must be given.
#
#        This can't be something derived from Stencil"""
