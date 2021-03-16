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
# StencilList; A class for collections of FDM stencils defined over some
# grid to get the full operator on a Grid. 
#

import Stencil as _s 

import Field as _f 

import Grid as _g

import Utils as _u

# import Numeric/numarray!
from Utils import _n

import MatSparse
import Numeric

StencilSetError="Error in StencilSet"

class StencilSet:
   """
      A class for collections of Stencils defined over some Grid.
      Explicit case: One iteration will be u = S(u) where u is the unknown
      field over the Grid, and S is the StencilList

      Datastructure:

      collection: Collection of Stencils. This is a list, where the Stencil
                  with given index should be used in nodes with the same
                  stencilindicator.

      where:      Collection of iterators stored in the same order as the
                  stencils. The where-iterator return a set of tuples where the
                  stencil should be applied.
   """

   def __init__(self, grid):
      # attach grid:
      self.grid = grid
      self.collection = []
      self.where = []
      self.call = self.call_buildmatrix
#      self.call_matrix = self.call_matrix_numarray
#      self.buildMatrixOperator = self.buildMatrixOperator_numarray
      self.call_matrix = self.call_matrix_fastmatsparse
      self.buildMatrixOperator = self.buildMatrixOperator_fastmatsparse

   def call_matrix_numarray(self,u):
       ureturn = _f.Field(u.grid)
       ureturn.data[:] = _n.dot(self.A,u.data)
       ureturn.data += self.f
       ureturn.updateField()
       return ureturn

   def call_matrix_fastmatsparse(self,u):
       ureturn = _f.Field(u.grid)
       self.tmpud[:] = u.data[:]
       #self.tmpd[:] = self.A*self.tmpud
       self.A.prod(self.tmpud,self.tmpd)
       self.tmpd += self.f
       ureturn.data[:] = self.tmpd[:]
       ureturn.updateField()
       return ureturn

   def direct_matvec(self,u):
       """Assume Numeric array u as input.
          The user have to handle stuff the right way...
          Does not work in parallel.
       """
       
       r_u = self.A*u
       r_u += self.f
       return r_u

   def call_buildmatrix(self,u):
       self.buildMatrixOperator(u)
       self.call = self.call_matrix
       return self.call_matrix(u)

   def call_updatesource(self,u):
       self.buildMatrixOperator_updateSource(u)
       self.call = self.call_matrix
       return self.call_matrix(u)

   def updateDataStructures(self):
       self.call = self.call_buildmatrix

   def updateSourceDataStructures(self):
       self.call = self.call_updatesource

   def useIteratorMultiplication(self):
       """A*u and A(u) performed with the basic iterators on the Grid"""
       self.call = self.old_call_algorithm

   def useMatrixMultiplication(self):
       """A*u and A(u) performed by matrix-vector multiplication"""
       self.call = self.call_buildmatrix

   def __call__(self,u):
      """Apply this StencilList to the supplied field u"""
      return self.call(u)

   def old_call_algorithm(self,u):
      # For each stencil in the collection, apply the stencil to the
      # corresponding nodes in the field given in the where object.

      # create a new field for the return-values(we can not update u directly,
      # that may violate the method)
      
      #TODO: make it possible to use a predef. field for the result, should 
      # avoid a few Field creations... ?

      ureturn = _f.Field(u.grid)

      # It should only be necessary to do this loop once, at least if 
      # coefficients functions are non-variable in time. 
      # For time-dependent problems we need to do this once for each timestep.
#      if not self.usematrix:
#          # fill self.A as the operator-matrix
#          # fill self.f as the souce-part. This is a vector/field. 
#      else:
#          # use self.A as operator instead of looping.
      
      for i in xrange(len(self.collection)):
         stencil  = self.collection[i]
         where_item    = self.where[i]
         # where_item is some tupleIterator, make sure it is reset (for timeloops)
         where_item.reset()
         #print "Running stencil ",i," with bounds: ",where_item.getBounds()
         # TODO: reimplement this with a map?
         # something like this probably:
         # ureturn = map(lambda i: stencil(u,i), where_item)
         # but I am not sure an iterator can be used in map, so stick 
         # to for-loop for now
         for ind in where_item:
            ureturn[ind] = stencil(u,ind)

      # before we return Field, make sure ghosts are up-to-date in case of 
      # a parallel machine:
      ureturn.updateField()
      return ureturn

   def buildMatrixOperator_fastmatsparse(self,u):
       ## TODO: is it possible that the same index can be in multiple itereators 
       ## The effect is that there is multiple stencils defined for one point. 
       ## If so, they should probably be added up. But the original version just 
       ## overwrite (see above in call), so I do that here, too.

       # Create a double matrix.
       #self.A = _n.zeros((u.datasize,u.datasize),'d')
       Ama  = MatSparse.MapMatSparse()
       # and a double vector for the source.
       self.f = Numeric.zeros((u.datasize),'d')
#       print u.myrank,">Datasize: ",u.datasize
#       print u.myrank,">basenode: ",u.basenode
#       print u.myrank,">shapemul: ",u.shapemultiplicator
       
       # Create storage for temporary vectors during computation.
       self.tmpd = Numeric.zeros((u.datasize),'d')
       self.tmpud = Numeric.zeros((u.datasize),'d')

       # Make sure that there is at least one zero in each row:
       for i in xrange(u.datasize):
           Ama[i,i] = 0.0

       self.allinds = {}
       for i in xrange(len(self.collection)):
           s = self.collection[i]
           w = self.where[i]
           w.reset()

#           print u.myrank,">>> Bounds: ",w.getBounds()
#           for ind in w:
#               print u.myrank,">>> ind: ",ind
#           w.reset()
           for ind in w:
#               try:
               (d,s_v) = s.effect(u,ind)
#               except:
#                   print u.myrank,">>> Bounds: ",w.getBounds()
#                   print u.myrank,">>> name:   ",type(w).__name__
#                   w.reset
#                   for tmpind in w:
#                       print tmpind
#                   raise StencilSetError, "%s>>> Can not get effect for given index, %s" % (u.myrank,ind)

               dataindex = u.indexToDataIndex(ind)
               # store all inds. and stencil-ref. for source-update:
               #self.allinds[ind] = (dataindex,s)
               self.allinds[dataindex] = (ind,s)
               for key_di in d:
                   #print "i: ",dataindex, "j: ",key_di,", val: ",d[key_di]
                   Ama[dataindex,key_di] = d[key_di]
               #print u.myrank,">f: ind, i: ",ind,dataindex, ", val: ",s_v
#               try: 
               self.f[dataindex] = s_v
#               except:
#                   dilocal = u.localIndexToDataIndex(ind)
#                   diglob  = u.globIndexToDataIndex(ind)
#                   print u.myrank,">Datasize: ",u.datasize
#                   print u.myrank,">basenode: ",u.basenode
#                   print u.myrank,">shapemul: ",u.shapemultiplicator
#                   print u.myrank,">grid: ",u.grid.sizespec, u.grid.geometry
#                   raise StencilSetError, "%s> unable to set f[%s] for ind %s with value %s, %s or %s" % (u.myrank,dataindex,ind,s_v,dilocal,diglob) 

#       self.alldi_inds = self.allinds.keys()
#       self.alldi_inds.sort()
#       self.allindicies = []
#       self.allstencils = []
#       for i in self.alldi_inds:
#           self.allindicies.append(self.allinds[i][0])
#           self.allstencils.append(self.allinds[i][1])
       self.A = MatSparse.FastMatSparse(Ama)

   def buildMatrixOperator_numarray(self,u):
       ## TODO: is it possible that the same index can be in multiple itereators 
       ## The effect is that there is multiple stencils defined for one point. 
       ## If so, they should probably be added up. But the original version just 
       ## overwrite (see above in call), so I do that here, too.

       # Create a double matrix.
       self.A = _n.zeros((u.datasize,u.datasize),'d')
       # and a double vector for the source.
       self.f = _n.zeros((u.datasize),'d')

       for i in xrange(len(self.collection)):
           s = self.collection[i]
           w = self.where[i]
           w.reset()
           for ind in w:
               (d,s_v) = s.effect(u,ind)
               dataindex = u.indexToDataIndex(ind)
               for key_di in d:
                   self.A[dataindex,key_di] = d[key_di]
               self.f[dataindex] = s_v

   def buildMatrixOperator_updateSource(self,u):

#       # use stored indecis and stencil-refs.
#       #for ind in self.allinds:
       for di in self.allinds:
           #(di,s) = self.allinds[ind]
           (ind,s) = self.allinds[di]
           s_v = s.effect_source(u,ind)
           self.f[di] = s_v

 
#       #self.f[self.alldi_inds] 
#       newvals = _n.zeros(len(self.f),'d')
#       newvals[self.alldi_inds] = [s.effect_source(u,ind) for (s,ind) in zip(self.allstencils,self.allindicies)]
#       self.f[:] = newvals
#       for i in xrange(len(self.alldi_inds)):
#           self.f[self.alldi_inds[i]] = self.allstencils[i].effect_source(u,self.allindicies[i])
#       print self.f
#       self.f[:] = newvals
#       self.f[self.alldi_inds] = [s.effect_source(u,ind) for (s,ind) in zip(self.allstencils,self.allindicies)]
#       for i in xrange(len(self.collection)):
#           s = self.collection[i]
#           w = self.where[i]
#           w.reset()
#           for ind in w:
#               s_v = s.effect_source(u,ind)
#               dataindex = u.indexToDataIndex(ind)
#               self.f[dataindex] = s_v

   def doInitParallel(self):
       """Make stencillist parallel
          This method check for a parallel grid, and in that case, the stencillist
          is adjusted. Most important, the 'where' entries are modified to represent
          nodes present on the subgrid for the current process. Kind of...

          We assume that proper iterators with the right methods, like getBounds 
          are used.
       """
        
       if self.grid.isParallel:
           grid_bounds = self.grid.myGlobalNodalBounds
           # TODO: Remove debug print
           myrank = self.grid.myrank
           #print myrank,"> grid bounds: ", grid_bounds

           # build a removal list, pop these indices later on
           remove_list = []
           for i in xrange(len(self.where)):
               w = self.where[i]
               b = w.getBounds()
               #print myrank,"> Got bounds: ",b

               # TODO: Remove debug print
               #print " Got bounds: ",b
               # check that grid_bounds and b have non-zero intersection:
               if not _u.boundsIntersect(b,grid_bounds):
                   remove_list.append(i)
                   # TODO: Remove debug print
                   #print "  Remove b"
               else:
                   # make a copy of the bounds, listify it for assignment:
                   #nb = list(_u.deepcopy(b))
                   nb = []
                   for tu in b:
                       nb.append(list(tu))
                   #for j in xrange(len(b)): # should be just j=0 and j=1...
                   for i in xrange(len(b[0])): # same as i=0,...,nsd-1
                       if b[0][i] < grid_bounds[0][i]:
                           nb[0][i] = grid_bounds[0][i]
                       if b[1][i] > grid_bounds[1][i]:
                           nb[1][i] = grid_bounds[1][i]
                   # TODO: implement this setBounds method
                   # set new Bounds in the iterator.
                   # The iterator have to handle this correctly. For instance, if w is
                   # a corner or boundary iterator, the new bounds may not make much 
                   # sense, as the corner og boundary iterator only should relate
                   # to the physical boundaries of the grid.
                   # TODO: Remove debug print
                   #print myrank,"> Set new bounds: ",nb
                   w.setNewBounds(nb,iOB=True)

           remove_list.sort()
           remove_list.reverse()
           for item in remove_list:
               # TODO: Remove debug print
               #print "Remove item ",item
               self.collection.pop(item)
               self.where.pop(item)

       #TODO: Remove debug print
       #print "All stencils and wheres in stencillist"
       #for i in xrange(len(self.collection)):
       #    print "Stencil: ",self.collection[i].debugPrint()
       #    print "Where: ",self.where[i].getBounds()
                   


   def debugPrint(self):
       for i in xrange(len(self.collection)):
           print "Stencil: ",self.collection[i].debugPrint()
           print "Where: ",self.where[i].getBounds()

   def __mul__(self,v):
      """If stencillist is used as a discrete operator, like in L_d u = f, it is 
         convenient to use the multiplication notation, instead of the 
         call-notation. mul is therefor just another way of calling!
      """

      if isinstance(v,_f.Field):
          return self(v)
      elif isinstance(v,Numeric.ArrayType):
          return self.direct_matvec(v)
      elif isinstance(v,_n.ArrayType):
          return _n.array(self.direct_matvec(v))
      else:
          raise StencilSetError, "Multiplication with unknown operand"
      

   def addStencil(self,s,w):
      # attach the stencil to the collection object
      self.collection.append(s)
      # attach the where:
      self.where.append(w)
      # return the index for this stencil:
      return len(self.collection)-1
      

   def __iadd__(self,other):
      """Operator: add another stencillist into this."""

      # TODO: Check that the two stencillists are created for the same grid?
      #
      # probably self.grid == other.grid is ok, but id(self.grid) == id(other.grid)
      # can also be used.
      self.collection += other.collection
      self.where += other.where

      return self

   def __getitem__(self,i):
      """Get the stencil at index i from the StencilList"""

      return self.collection[i]

   def __setitem__(self,i,s):
      """Set the entry i in the stencillist to a new Stencil s

      The index i must be in the range of the already built stencillist"""

      if i < len(self.collection)-1:
         self.collection[i] = s

   def stencilListWidth (self):
       """Compute the width in each direction of this stencillist
          Return: a tuple of tuples, one tuple for each direction, containing two 
          values, (l,r), where l is the maximum negative offset in that direction 
          and r is the maximum positive offset in that direction.
       """

       l_master = [0 for x in xrange(self.grid.nsd)]
       r_master = [0 for x in xrange(self.grid.nsd)]

       for s in self.collection:
           (l,r) = s.stencilWidth()
           for i in xrange(len(l)):
               l_master[i] = min(l_master[i],l[i])
               r_master[i] = max(r_master[i],r[i])
       return (tuple(l_master),tuple(r_master))

   def signature(self):
       """Build the complete signature for the whole stencillist
          All keys (offsets) are fetched from the attached stencils, not 
          taking the where-part into account, and build a new list of keys,
          hence describing the offset pattern for the complete list of 
          stencils. 

          This is to be used for describing communication patterns, and is hence 
          by no means and optimal solution, as the where-parts should preferably 
          be considerd when describing communication on each subpart.
       """

       signature = {}
       for s in self.collection:
           for k in s.stencil_coeff:
               if not signature.has_key(k):
                   signature[k] = 1

       return signature



# create an alias: StencilList (for backward compatibility)
StencilList = StencilSet
