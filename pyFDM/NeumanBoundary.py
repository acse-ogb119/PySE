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
# NeumanBoundary; Implementation of necesasry utils for building a 
#                 Neuman Boundary condition
#

import Utils as _u
import StencilList as _sl
import copy as _c


NeumanError = "Error in the Neuman condition creator"

# The neuman condition: Based on a Stencil on a given grid, a boundary
# description (may be a part of the boundary) and the exact neuman condition, 
# a stencillist with the necessary stencils for the neuman condition will 
# be build. This stencillist should be merged with the applications stencillist.


# Algorithm: For each boundarypart, make a copy of the stencil, modify it 
# according to grid-information and the neuman condition, and add it to the 
# stencillist. At last, return the stencillist.

def createNeumanBoundary(stencil,grid,condition,bounds=None, region=None,\
        type='box', direction='in', center=None, radius=None):
   """Create a neumann stencil list for given boundary and grid.
      
   Build a stencillist for a neuman boundary condition from the supplied
   argumets: A stencil, a grid, a boundary specification and a neuman
   condition. The boundary specification should be given as
   ([x0,y0,z0,...],[x1,y1,z1,...]), wher the box-region defined by this range
   is assumed to be a continous boundary for the grid.
   
   Region specify where the condition applies in terms of coordinates. 

   If neither region nor bounds is given, assume the whole grid.

   If both region and bounds are given, region will be used.
   """
   if not region and not bounds:
       bounds = grid.range()

   if region:
       bounds = grid.regionToBounds(region)

   nsd  = grid.nsd

   # the local stencil list
   neumannStencilList = _sl.StencilList(grid)

   # Check the boundary spec against nsd:
   if ( len(bounds[0]) != nsd and len(bounds[1]) != nsd ):
      raise NeumanError, "nsd in grid and boundaryspec does not match"

   # Create a controllist for edgeGenerator, to indicate which ones are 
   # real (and hence wanted) boundaries:
   cl = list(_u.getZeroTuple(bounds[0]))
   gr = grid.range()
   # if flag is False after this, bounds do not includes any boundary at all.
   flag = False
   for i in xrange(grid.nsd):
       # "left" side:
       if bounds[0][i] == gr[0][i]:
           cl[i] += 1
           flag = True
       # "right" side:
       if bounds[1][i] == gr[1][i]:
           cl[i] += 2
           flag = True
   
   if not flag:
       # there is no boundary, return empty stencillist.
       return neumannStencilList

   # create the edgeIterator:
   #edgeit = _u.edgeIteratorGenerator(nsd,bounds[0],bounds[1])
   edgeit = _u.edgeGenerator(nsd,bounds[0],bounds[1],controlList=cl,yieldAllInfo=True)

   for edge in edgeit:
      # Each edge is a list with three entries: the edgeIterator, a base-tuple
      # with ones at the boundary-indices and a tuple with the index-values for
      # the boundary indices. These must be checked against the grid to check
      # which nodes in the stencil that are ghost on this particular boundary.

      # Start with at fresh copy of the stencil
      sc = _c.copy(stencil)

      # Modify stencil according to condition.
      # edge[1] shows the fixed part, edge[2] show the actual fixed-values
      #print "edgestuff: ",edge[1],", ",edge[2]
      sc.modForNeumann(edge[1],edge[2],condition,grid)

      # Add stencil to the neumannList
      # edge[0] is a iterator, and hence the 'where' part.
      #
      # If something different than a box-type region is specified,
      # the iterator must be converted to the requested type (and possibly 
      # restricted in terms of indicies), before added to the stencillist.
      # Remark that the set of indicies may be empty, in that case we can leave 
      # the iterator out.
      if type == 'box':
          neumannStencilList.addStencil(sc,edge[0])
      elif type == 'circle':

          if not center:
              raise NeumanError, "Center must be given for a circle region"
          if not radius:
              raise NeumanError, "Radius must be given for a circle region"

          #print "Create circleIterator for center: ",center,", radius: ",radius,", direction: ",direction
          cit = _u.circleIterator(edge[0],grid,center,radius,direction)

          # Only add iterator if there are any indicies in the list.
          if len(cit.indicies) > 0:
          #    print "Add circleIterator for neumann boundary, ",cit.getBounds()
              neumannStencilList.addStencil(sc, cit)
      else:
          raise NeumanError, "unknown boundary region type specified (%s)." % (type)

   return neumannStencilList


# set an alias for backward compatibility
def create(s, b, g, c):
    return createNeumanBoundary(stencil=s, grid=g, condition=c, bounds=b)
