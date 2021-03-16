# -*- coding: iso-8859-1 -*-
#
# (c) Copyright 2003, 2004, 2005
#     Authors: Åsmund Ødegård, Ola Skavhaug
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
# Stencil; A stencil class for Finite Difference Methods in n space
# dimensions.
#



StencilError = "Error in Stencil"

import operator
_o = operator
import copy
_c=copy
import Utils as _U
import inspect

class Stencil(object):
   """A stencil class"""

   def __init__(self,nsd=2,varcoeff=0,source=0.0,nodes=None):
      """Initialize Stencil with from zero to 3 arguments:
       nsd, varcoeff (a boolean) and source. The stencil is stored as a
       dictionary. 
       The source term should be used as a field independent value -
       when applied for an implicit scheme, the source-term goes into the right
       hand side, while it should be added togehter with everything else as the
       last thing in an explicit scheme. It may be space-dependent. If so, it
       is a function
       
       Datastructure.

       nsd (int):       
           Contain the number of space dimensions this Stencil is used for

       varcoeff (bool):  
           True if the stencil have variable coefficients either in the nodes
           or in the source. 

       source (value or lambda function): 
           Contain either the source value or the source function-pointer.

       nodes : If given, a dictionary with stencil points and values

       __init__ creates:
       stencil_coeff (dictionary of lists of functions or values):
           The heart of the stencil. All offsets for the node, relative to the
           center is stored as the key. If we have fixed coefficients in the
           stencil, the value for each key (offset) is simply a number. Else
           the value for each key is a list of functions. When we apply the
           stencil to some field, the point is used as argument for each
           function and the resulting values are added.
       
      """ 
      self.nsd = nsd
      self.varcoeff = varcoeff
      if self.varcoeff:
         # source may be a function-list (apply args over list to get result).
         #if type(source) == type(1) or type(source) == type(1.0):
         if isinstance(source,(int,float)):
            #source is a number, but we need a function:
            sourcefu = constfunc(source)
            self.source = [sourcefu]
         elif isinstance(source,(tuple,list)):
             self.source = source
         elif type(source).__name__ == 'function':
             self.source = [source]
         else:
             raise StencilError, "Variable coefficients, but try to use source of unknown type in constructor, type(source)=%s" % (type(source))

      else:
         self.source = source

      self.stencil_coeff = {}

      if self.varcoeff:
         self.addNode = self.addNode_vc
         self.addSource = self.addSource_vc
         self.call = self.call_vc
         self.effect = self.effect_vc
         self.effect_source = self.effect_source_vc
      else:
         self.addNode = self.addNode_nvc
         self.addSource = self.addSource_nvc
         self.call = self.call_nvc
         self.effect = self.effect_nvc
         self.effect_source = self.effect_source_nvc
      
      if nodes:
          for n in nodes:
              self.addNode(n, nodes[n])

   def debugPrint(self):
      """Print usefull information about this stencil"""

      print "nsd: %d, variable. coefficients: %d" % (self.nsd,self.varcoeff)
      if self.varcoeff:
          print "Eval in origo"
          args = _U.getZeroTuple(range(self.nsd))
          for key in self.stencil_coeff:
              print "key: %s => %s" % (key,reduce(_o.add,map(lambda i: i(*args),self.stencil_coeff[key])))
          print "source: %s" % (reduce(_o.add,map(lambda i: i(*args),self.source)))
          #print "source: %s" % (reduce(_o.add,map(lambda i: i(*args),self.source)))
      else:
          for key in self.stencil_coeff.keys():
              print "key: %s => %s" % (key,self.stencil_coeff[key])
          print "source: %s" % (self.source)


   def addNode_nvc(self,key,value):
      """Add a coefficient element to the stencil. If element is already
         defined, just add up the value"""

      #check length of key:
      if  type(key) == type(1) and self.nsd != 1:
         raise StencilError, "Length of key incompatible with nsd"
      elif  type(key) != type(1) and len(key) != self.nsd:
         raise StencilError, "Length of key incompatible with nsd"

      if self.stencil_coeff.has_key(key):
         self.stencil_coeff[key] += value
      else:
         self.stencil_coeff[key] = value

   def addSource_nvc(self,source):
      """Add a value to the source-entry in the stencil"""
      self.source += source


   def addNode_vc(self,key,funclist):
      """Add a variable coefficient element to the stencil. If the element is
         already defined, add the functions
         
         stencil_coeff[key] will in this case be a list of functions."""
      
      # check length of key:
      if type(key) == type(1) and self.nsd != 1:
         raise StencilError, "Length of key incompatible with nsd"
      elif type(key) != type(1) and len(key) != self.nsd:
         raise StencilError, "Length of key incompatible with nsd"

      # check that funclist is actually a list:
      # TODO: we should be able to handle the following cases:
      #  1) funclist is no list at all, just a single function. Create list
      #     and add apropriately
      #  2) funclist is no function either, just a number. If so, create a 
      #     constanct function, put in a list, and add.
           
      #if type(funclist) != type([1]):
      #   raise StencilError, "Funclist must be a list"
      if isinstance(funclist,(int,float)):
          # stencil coeff. given as a value. 
          # convert to a constant function.
          funclist = [constfunc(funclist)]

      if not isinstance(funclist,(list,tuple)):
          funclist = [funclist]

      if self.stencil_coeff.has_key(key):
         self.stencil_coeff[key] += funclist
      else: 
         self.stencil_coeff[key] = []
         self.stencil_coeff[key] += funclist

   def addSource_vc(self,funclist):
      """add functions to the source list. Same as for nodal values."""
#      if type(funclist) != type([1]):
#          raise StencilError, "Funclist to source must be a list"

      if inspect.isfunction(funclist):
          self.source += [funclist]
      if type(funclist).__name__ == 'function':
          self.source += [funclist]
      elif isinstance(funclist,(list,tuple)):
          self.source += funclist
      elif isinstance(funclist,(int,float)):
          self.source += constfunc(funclist)
      elif type(funclist).__name__ == 'Field':
          # TODO:
          # funclist is given as a (pyFDM) Field. 
          # We should be able to use the values
          # There are a lot of things to sort out regarding this, so the 
          # functionality is just dummy for now.
          self.source += [funclist]
          raise StencilError, "Funclist specified as a Field. This is still on the TODO list"
      else:
          raise StencilError, "Funclist to source must be a list, tuple or scalar. Unknown type used: %s" % (type(funclist))


   def resetNode(self,key,value):
      del self.stencil_coeff[key]
      self.addNode(key,value)

   
   def applyPointwise(self,key,args=None):
     
      # Remark: We now unroll the args when calling!
      if not self.stencil_coeff.has_key(key):
         raise StencilError, "node is not set"

      if self.varcoeff:
         if not args:
            raise StencilError, "Need function arguments for variable coeff"
         # apply the arguments to the functionlist:
         if len(self.stencil_coeff[key]) > 1:
            if type(args) == type(1):
               args = (args,)
            return reduce(_o.add,map(lambda i: i(*args),self.stencil_coeff[key]))

         else: 
            return self.stencil_coeff[key][0](*args)
      else:
         return self.stencil_coeff[key]

   def applySource(self,args=None):
      
      if self.varcoeff:
         if not args:
            raise StencilError, "Need function arguments for variable coeff"

         if len(self.source) > 1:
            if type(args) == type(1):
               args=(args,)
            #return reduce(_o.add,map(lambda x,a=args: apply(x,(a,)),self.source))
            # Before we started rolling out arguments, such that user sees x,y,z
            #return reduce(_o.add,map(lambda i: i(args,),self.source))
            # After (i.e., with rolling out args)
            return reduce(_o.add,map(lambda i: i(*args),self.source))

         else:
            #return self.source[0](args)
            return self.source[0](*args)
      else:
         return self.source


#   def start(self):
#      """set the iterator to the beginngin of the stencil"""
#      self.iterator=-1
#
#   def next(self):
#      """get the next element in the stencil"""
#      self.iterator += 1
#      i = self.iterator
#      if i < len(self.stencil_coeff):
#         return self.stencil_coeff[i]

   def __add__(self,other):
      """add two stencils. State of varcoeff are checked in both
         
         If both stencils are non-varcoeff, a new non-varcoeff is returned.
         Else, a new, varcoeff-stencil is returned. Values in non-varcoeff
         stencils are converted to constant functions.
      """

      if self.nsd != other.nsd:  
         raise StencilError, "Incompatible number of space dimensions"

      # not variable coefficients:
      if not self.varcoeff and not other.varcoeff:
         #create a new Stencil:
         adds = Stencil(self.nsd,self.varcoeff,self.source)
         # add all elements in self to new:
         for key in self.stencil_coeff.keys(): 
            adds.addNode_nvc(key,self.stencil_coeff[key])
         # add other
         for key in other.stencil_coeff.keys():
            adds.addNode_nvc(key,other.stencil_coeff[key])

         # at last, add source-fields:
         adds.addSource(other.source)

         return adds
      # with variable coefficients in both of the stencils:
      elif self.varcoeff and other.varcoeff:
         # create a new Stencil
         adds = Stencil(self.nsd,1,self.source)


         # add all elements in self to new:
         for key in self.stencil_coeff.keys():
            adds.addNode_vc(key,self.stencil_coeff[key])

         # add other:
         for key in other.stencil_coeff.keys():
            adds.addNode_vc(key,other.stencil_coeff[key])

         # at last, add source-field
         adds.addSource(other.source)
         return adds
      # one of the stencil has variable coefficients:
      elif self.varcoeff or other.varcoeff:
         # create a new Stencil with var. coeff
         adds = Stencil(self.nsd,1)
         # which one has var. coeff:
         vcpointer = self
         nvcpointer = other
         if not vcpointer.varcoeff: 
            vcpointer = other
            nvcpointer = self

         # add the variable coeff
         for key in vcpointer.stencil_coeff.keys():
            adds.addNode_vc(key,vcpointer.stencil_coeff[key])

         # for the non-variable coeff, create a constant function and add:
         for key in nvcpointer.stencil_coeff.keys():
            value = nvcpointer.stencil_coeff[key]
            addfunc = constfunc(value)
            adds.addNode_vc(key,[addfunc])

         # add both sources. The const.coeff one must be made variable:
         nvcsource = [constfunc(nvcpointer.source)]
         adds.addSource(nvcsource)
         adds.addSource(vcpointer.source)
         return adds
      else:
         raise StencilError, "Unimplemented: add varcoeff and non-varcoeff stencil"
         
   def __sub__(self, other):
      """subtract two stencils. 
         We do this by first multiplying other with -1, then calling add.
      """
      return self + other.__rmul__(-1.0)
  
   def __rmul__(self,other):
      """
         Multiply a stencil with some scalar valye. Return the freshly created 
         stencil. That is, newStencil = 4.0*stencil
      
         Original behavior: 
         Multiply a stencil with some scalar value. Return self; i.e. change
         current stencil instance
           - This origianl behavior is not pythonic. I can not modify self here.
             The inplace operator should be implemented for that. 
      """

      tmps = _c.copy(self)

      # Constant coefficients
      if not tmps.varcoeff:
          for key in tmps.stencil_coeff.keys():
             tmps.stencil_coeff[key] *= other
          tmps.source *= other
      # Variable coefficients
      else:
          #raise StencilError, "scalar*stencil not implemented for var. coeff"
          for key in tmps.stencil_coeff:
              tmps.stencil_coeff[key] = [_U.func_mulScalarFunclist(other, tmps.stencil_coeff[key])]
          tmps.source = [_U.func_mulScalarFunclist(other, tmps.source)]

      return tmps
         
   def __mul__(self,other):
       """
          Same as __rmul__, only the other way around: 
          newStencil = stencil*4.0

          Only int/float is supported as 'other', multiplication of stencils
          are not supported as I'm not sure what the semantic of that should be.
       """ 
       
       if isinstance(other,(int,float)):
           return self.__rmul__(other)
       else:
           raise StencilError, "Stencil only supports multiplication with scalars"
      
   def __iadd__(self,other):
      """add other stencil to this one"""
      if self.nsd != other.nsd:
         raise StencilError, "Incompatible number of space dimensions"
      
      if not self.varcoeff and not other.varcoeff:
         for key in other.stencil_coeff.keys():
            self.addNode_nvc(key,other.stencil_coeff[key])
         self.addSource(other.source)
         return self

      elif self.varcoeff and other.varcoeff:
         for key in other.stencil_coeff.keys():
            self.addNode_vc(key,other.stencil_coeff[key])
         self.addSource(other.source)
         return self
      # if self has var.coeff and other not, generate constfunc. objects for
      # keys in other and add to self. return self
      elif self.varcoeff and not other.varcoeff:
         for key in other.stencil_coeff.keys():
            value = other.stencil_coeff[key]
            addfunc = constfunc(value)
            self.addNode_vc(key,[addfunc])
         nvcsource = [constfunc(other.source)]
         self.addSource(nvcsource)
         return self
      # if self has not var.coeff while other has, generate constfunc objects
      # for keys in self and add to other. return other.
      elif not self.varcoeff and other.varcoeff:
         for key in self.stencil_coeff.keys():
            value = self.stencil_coeff[key]
            addfunc = constfunc(value)
            other.addNode_vc(key,[addfunc])
         nvcsource = [constfunc(self.source)]
         other.addSource(nvcsource)
         return other
            
   def __call__(self,u,i):
      return self.call(u,i)

   def effect_vc(self,u,i):
       # return a dict. with key: dataindex, value: stencilvalue.
       # and the source-value.
       d = {}
       for key in self.stencil_coeff:
           ti = map(_o.add,i,key)
           x = u.grid.indexToPoint(ti)
           d[u.globIndexToDataIndex(ti)] = self.applyPointwise(key,x)
       v = self.applySource(u.grid.indexToPoint(i))
       return (d,v)

   def effect_source_vc(self,u,i):
       return self.applySource(u.grid.indexToPoint(i))

   def call_vc(self,u,i):

      # May try to remove the for-loop with clever map/reduce stuff.
      # Starting-point: 
      # map(lambda k: tuple(map(lambda x,y: x+y,i,k)),map(lambda key: key,s.stencil_coeff))
      # maybe something like this:
#      result = reduce(_o.add,map(lambda ti,kk: u[ti]*self.applyPointwise(kk,u.grid.indexToPoint(ti)), map(lambda k: map(lambda x,y: x+y,i,k),self.stencil_coeff),self.stencil_coeff))
#
#      Only problem: fails if stencil does not have any keys. (e.g. a dirichlet boundary)
#
#      Maybe we can fix by add the source at the end of the reduce:
#      
#      result = reduce(_o.add,map(lambda ti,kk: u[ti]*self.applyPointwise(kk,u.grid.indexToPoint(ti)), map(lambda k: map(lambda x,y: x+y,i,k),self.stencil_coeff),self.stencil_coeff),self.applySource(u.grid.indexToPoint(i)))
#
#      it's maybe clever, but actually performs worse. So, until we find something better, we leave the
#      old implementation in place:
#
#      old implementation:
      result = 0.0
      for key in self.stencil_coeff:
         # offset copy of i relative to key:
         ti = map(_o.add,i,key)

         # get the point-values for this index:
         x = u.grid.indexToPoint(ti)

         # eval the stencil_coeff function here, and multiply with u:
         result += u[ti] * self.applyPointwise(key,x) 

      # add the source, evaluted in the center:
      #x = u.grid.indexToPoint(i)
      result += self.applySource(u.grid.indexToPoint(i))

##      # Try 3: not better either...
#      keys = [ key for key in  self.stencil_coeff ]
#      tis = map(lambda key: map(_o.add,i,key),keys)
#      xis = u.grid.indicesToPoints(tis)
##      values = map(lambda ti,key,x: u[ti]*self.applyPointwise(key,x),tis,keys,xis)
#
##      result = reduce(_o.add,values,self.applySource(u.grid.indexToPoint(i)))
#      result = reduce(_o.add,map(lambda ti,key,x: u[ti]*self.applyPointwise(key,x),tis,keys,xis),self.applySource(u.grid.indexToPoint(i)))
#
#      Try 4: This one is better, but the old impl. is still best :-(
#      result = reduce(_o.add,[u[ti]*self.applyPointwise(key,u.grid.indexToPoint(ti)) for (ti,key) in [(map(_o.add,i,key),key) for key in self.stencil_coeff]],self.applySource(u.grid.indexToPoint(i)))


      return result

   def effect_nvc(self,u,i):
       # return a dict. with key: dataindex, value: stencilvalue.
       # and the source-value.
       d = {}
       for key in self.stencil_coeff:
           ti = map(_o.add,i,key)
           d[u.globIndexToDataIndex(ti)] = self.stencil_coeff[key]
       v = self.source
       return (d,v)

   def effect_source_nvc(self,u,i):
       return self.source

   def call_nvc(self,u,i):
      """call for nonvariable coefficients. 
         Arguments: a field and a index-tuple
         Return: this stencil applied to the field in the given index.
      """

      result=0.0
      for key in self.stencil_coeff.keys():
         # offset copy of i relative to key:
         ti = self.addtuples(i,key)
         result += u[ti] * self.stencil_coeff[key]
      result += self.source

      return result
         
   def addtuples(self,tn,tm):
      ti = []
      #assume equal lengths:
      for i in range(len(tn)):
         ti.append(tn[i] + tm[i])
      return ti

   def __copy__(self):
      c=Stencil(self.nsd,self.varcoeff,_c.copy(self.source))
      # Do I really need deepcopy here? It is a problem for functional
      # values...
      #c.stencil_coeff = _c.deepcopy(self.stencil_coeff)
      #c.stencil_coeff = _c.copy(self.stencil_coeff)
      for key in self.stencil_coeff:
          # make sure that a new list is created.
          #c.stencil_coeff[key] = []
          c.stencil_coeff[key] = _c.copy(self.stencil_coeff[key])
      return c

   def modForNeumann(self,edges,indices,condition,grid):
      """Modify the current stencil according to the given condition. 
         
         Input:
         edges is a tuple of same dimension as the grid, with ones at the
         indices wich are actually on some edge, 
         
         indices are the actual index-values corresponding to the ones, we can
         check against the grid to see whether we are 'left' or 'right'.

         condition is the neumann-condition. This might be a number or a
         function. If it is a function, the arguments must be a point-tuple (of
         the same length as grid.nsd).

         grid is needed to grab the deltas in the grid.
         
         Output: Nothing.  This function transform the current stencil to an
         edge-stencil

         We assume that condition is given according to the varcoeff flag; if
         varcoeff=0, condition is a number, if varcoeff!=0, condition is a
         function.
      """

      edgeinds  = []
      leftright = []
      tmpi      = 0
      gr = grid.range()
      for i in xrange(len(edges)):
         if edges[i] == 1:
            edgeinds.append(i)
            #TODO: If grid doesn't start indices on 0, this should be grid-min.
            #in the actual drection.
            if indices[tmpi] == gr[1][i]:
               leftright.append(1)
            elif indices[tmpi] == gr[0][i]:
               leftright.append(-1)
            else:
                leftright.append(0)
            tmpi+=1

      # seek through stencil_coeff list, find all nodes which are ghost, apply
      # the neuman-condition.

      # New version - handle points with multiple ghost-indices (i.e., points 
      # that are ghosts in more than one dimension. - think left of and below a 
      # corner for instance). 
      # TODO: Handle such multiple-dimension ghosts for dimensions 
      # higher than 2!
      #
      # Algorithm: 
      #
      # 1. Find all "ghost-indices" in the key. Format: list where first entry is
      # the index of the first ghost-ind, second entry is index of the second etc.
      # Note: a len(ghostinds) == 0 indicate that key is not a ghost at all.
      #
      # 2. If len(ghostinds) == 1, proceed as before. But, if some of the other 
      # indices in the key are non-zero, offset calls to the condition-function 
      # (if it is a function) accordingly. Remove the key.
      #
      # 3. If len(ghostinds) == 2 and grid.nsd == 2, use the formula based on 
      # the taylor-expansion for corner-points. 
      #
      # 4. If len(ghostinds) >= 2 and grid.nsd >= 3, raise an exception, at least 
      # until the formula for higher dimension than 2 is developed.
      

      for key in self.stencil_coeff.keys():
          #print  "Now working on key ",key
          ghostinds = []
          for i in xrange(len(edgeinds)):
              edgei = edgeinds[i]
              #print "i :",i,", edgei: ",edgei
              if key[edgei] != 0:
                  if abs(leftright[i] + key[edgei]/abs(key[edgei])) == 2:
                      ghostinds.append((edgei,i))
                      #print "key: ",key," considered ghost (",edges,",",indices,")"
          #print len(ghostinds)," ghostinds in this key, ",ghostinds

          if len(ghostinds) == 0:
              #print "key: ",key," does not have ghosts, nothing to do"
              continue
          elif len(ghostinds) == 1:
              edgei = ghostinds[0][0]
              i = ghostinds[0][1]
              # This is a ghost-node, modify:
              #print "key: ",key," considered ghost in ",edgei," (",edges,",",indices,")"

              # flip sign in the key at the actual ghost-index (? Bad name)
              # Also, find distance from ghost to new contribution
              lkey=list(key)
              lkey[edgei] = -lkey[edgei]
              addkey=tuple(lkey)
              width = 2 * abs(lkey[edgei])

              # add contribution to self 
              contribCoeff = self.stencil_coeff[key]

              self.addNode(addkey,contribCoeff)


              #TODO: Remove this print
              #print "added node %s with value %s" % (addkey,contribCoeff)

              # TODO: offset where the condition function is called!
              # If some of the other values in key (apart from the one making it
              # a ghost) are non-zero, the condition function should be evaluated
              # in grid.division[i]*key[i] for this particular i.

              # If we have variable coefficients, we need to build a new
              # function wich includes a call to the condition function
              # In this case, "contribCoeff" may be a list of functions...
              # or it actually is. So we combine the contribCoeff with the
              # scalars, then we multiply with condition
              # (function-objects). Evaluation will be done later.
              # 
              # this one fail if varcoeff == 0, and condition happen to be
              # a function.


              if inspect.isfunction(condition):

                 # check other indices in key for non-zero values, and offset
                 # condition-call accordingly
                 offset_cond = []
                 useOffset=False
                 for j in xrange(len(key)):
                    if j == edgei:
                       offset_cond.append(0)
                    elif key[j] != 0:
                       offset_cond.append(key[j])
                       useOffset=True
                    else: 
                       offset_cond.append(0)

                 if not useOffset:
                     offset_cond = []
                 #print "Offset_cond: ",offset_cond

                 if not self.varcoeff:
                    # since condition is a function, stencil must be so too. 
                    # So, convert self to varcoeff-version. This will change the contribCoeff
                    # we got abobe, so lets re-set this variable.
                    self.convert_to_varcoeff()
                    contribCoeff = self.stencil_coeff[key]

                 if len(offset_cond) != 0:
                     
#                     condition_string = "lambda *x: condition("
#                     useargs = []
#                     for j in xrange(len(key)):
#                         print "j :",j
#                         if offset_cond[j] != 0:
#                             useargs.append('x[%s]+(%s)*%s' %(j,offset_cond[j],grid.division[j]))
#                             print "useargs: ",useargs
#                         else:
#                             useargs.append('x[%s]' % (j))
#                             print "useargs: ",useargs
#                     condition_string += ",".join(useargs)
#                     condition_string += ")"
#                     print "The new condition string is: ",condition_string
#                     condition_call = compile(condition_string,'string','eval')
#                     new_condition = eval(condition_call)

                     
                     tmp_condition = lambda X,oc=offset_cond,dx=grid.division: condition(*[x+o*d for (x,o,d) in zip(X,oc,dx)])
                     new_condition = lambda *x: tmp_condition(x)
                     self.addSource([_U.func_mulFuncFunc(_U.func_mulScalarFunclist(leftright[i]*width*grid.division[edgei],contribCoeff),new_condition)])
                 else:
                     self.addSource([_U.func_mulFuncFunc(_U.func_mulScalarFunclist(leftright[i]*width*grid.division[edgei],contribCoeff),condition)])
              else:
                 # Condition is numeric
                 if self.varcoeff:
                    # convert condition to function:
                    # No need to offset calls, as condition is a constant.
                    conditionfunc=constfunc(condition)
                    self.addSource([_U.func_mulFuncFunc(_U.func_mulScalarFunclist(leftright[i]*width*grid.division[edgei],contribCoeff),conditionfunc)])
                 else:
                    self.addSource(leftright[i]*contribCoeff*width*grid.division[edgei]*condition)
                    
              #TODO: Remove this debug print
              #print "added source %s, composed as this:" % (self.source)
              #print "contribCoeff: %e, width: %e, dx: %e, condition: %e" % (contribCoeff, width, grid.division[edgei],condition)

              # remove the original node.
              self.stencil_coeff.__delitem__(key)
              # No need to do anything more with this key:
              # TODO: If the key is "ghost in more than one direction" we may 
              # have just mapped a ghost into another. Say (-1,-1) => (1,-1), 
              # i.e., still a ghost... We can flag such things, and go through the
              # list again. But probably more care is needed, as the correct 
              # approximation to such nodes (as (-1,-1)) may need more threatment.
              #print "Key is removed, go on to next key"

          elif len(ghostinds) == 2 and grid.nsd == 2:
              #print "*** key with 2 ghostinds: ",key," ghostinds: ",ghostinds
              
              I = key[0]
              J = key[1]
              s_I = I/abs(I)
              s_J = J/abs(J)

              contribCoeff = self.stencil_coeff[key]
              if not self.varcoeff:
                  self.addNode((0,0),(-4./3)*contribCoeff)
                  self.addNode((-I,0),(4./3)*contribCoeff)
                  self.addNode((0,-J),(4./3)*contribCoeff)
                  self.addNode((-I,-J),(-1./3)*contribCoeff)
                  
                  #print "added lots of nodes...."

                  if inspect.isfunction(condition):
                      self.convert_to_varcoeff()
                      contribCoeff = self.stencil_coeff[key]
                      self.addSource(_U.func_mulFuncFunc(_U.func_mulScalarFunclist((8./3)*(s_I*grid.dx + s_J*grid.dy),contribCoeff),condition))

                      offset_cond = [0,-J]
                      tmp_condition_x = lambda X,oc=offset_cond,dx=grid.division: condition(*[x+o*d for (x,o,d) in zip(X,oc,dx)])
                      off_condition_x = lambda *x: tmp_condition_x(x)
                      offset_cond = [-I,0]
                      tmp_condition_y = lambda X,oc=offset_cond,dx=grid.division: condition(*[x+o*d for (x,o,d) in zip(X,oc,dx)])
                      off_condition_y = lambda *x: tmp_condition_y(x)
                      self.addSource(_U.func_mulFuncFunc(_U.func_mulScalarFunclist((2./3)*(-s_I*grid.dx),contribCoeff),off_condition_x))
                      self.addSource(_U.func_mulFuncFunc(_U.func_mulScalarFunclist((2./3)*(-s_J*grid.dy),contribCoeff),off_condition_y))
                      
                      #print "added lots of sources...."
                  else:
                      self.addSource(((8./3)*condition*(s_I*grid.dx + s_J*grid.dy)\
                              + (2./3)*-s_I*grid.dx*condition\
                              + (2./3)*-s_J*grid.dy*condition)*contribCoeff)
                      #print "added souce...."

              else:
                  if not inspect.isfunction(condition):
                      condition = constfunc(condition)
                  self.addNode((0,0),lambda *x: (-4./3)*contribCoeff(x))
                  self.addNode((-I,0),lambda *x: (4./3)*contribCoeff(x))
                  self.addNode((0,-J),lambda *x: (4./3)*contribCoeff(x))
                  self.addNode((-I,-J),lambda *x: (-1./3)*contribCoeff(x))

                  self.addSource(_U.func_mulFuncFunc(_U.func_mulScalarFunclist((8./3)*(s_I*grid.dx + s_J*grid.dy),contribCoeff),condition))

                  offset_cond = [0,-J]
                  tmp_condition_x = lambda X,oc=offset_cond,dx=grid.division: condition(*[x+o*d for (x,o,d) in zip(X,oc,dx)])
                  off_condition_x = lambda *x: tmp_condition_x(x)
                  offset_cond = [-I,0]
                  tmp_condition_y = lambda X,oc=offset_cond,dx=grid.division: condition(*[x+o*d for (x,o,d) in zip(X,oc,dx)])
                  off_condition_y = lambda *x: tmp_condition_y(x)
                  self.addSource(_U.func_mulFuncFunc(_U.func_mulScalarFunclist((2./3)*(-s_I*grid.dx),contribCoeff),off_condition_x))
                  self.addSource(_U.func_mulFuncFunc(_U.func_mulScalarFunclist((2./3)*(-s_J*grid.dy),contribCoeff),off_condition_y))

              self.stencil_coeff.__delitem__(key)
              #print "Key is removed, go on to next key"
              
          else:
              raise StencilError, "ghostinds >= 2 and/or nsd > 2, this is not implemented"


#
##
###
###  This is the old version of this code. Remove in due time!
###
##
#      for key in self.stencil_coeff.keys():
#         print "Now working on key ",key
#         #for i in xrange(len(edgeinds)):
#         i = 0
#         while i < len(edgeinds): 
#            edgei = edgeinds[i]
#            print "i :",i,", edgei: ",edgei
#            if key[edgei] != 0:
#               if abs(leftright[i] + key[edgei]/abs(key[edgei])) == 2:
#                  # This is a ghost-node, modify:
#                  print "key: ",key," considered ghost (",edges,",",indices,")"
#
#                  # flip sign in the key at the actual ghost-index (? Bad name)
#                  # Also, find distance from ghost to new contribution
#                  lkey=list(key)
#                  lkey[edgei] = -lkey[edgei]
#                  addkey=tuple(lkey)
#                  width = 2 * abs(lkey[edgei])
#
#                  # add contribution to self 
#                  contribCoeff = self.stencil_coeff[key]
#
#                  self.addNode(addkey,contribCoeff)
#
#
#                  #TODO: Remove this print
#                  print "added node %s with value %s" % (addkey,contribCoeff)
#
#                  # TODO: offset where the condition function is called!
#                  # If some of the other values in key (apart from the one making it
#                  # a ghost) are non-zero, the condition function should be evaluated
#                  # in grid.division[i]*key[i] for this particular i.
#
#                  # If we have variable coefficients, we need to build a new
#                  # function wich includes a call to the condition function
#                  # In this case, "contribCoeff" may be a list of functions...
#                  # or it actually is. So we combine the contribCoeff with the
#                  # scalars, then we multiply with condition
#                  # (function-objects). Evaluation will be done later.
#                  # 
#                  # this one fail if varcoeff == 0, and condition happen to be
#                  # a function.
#                  if inspect.isfunction(condition):
#                     if not self.varcoeff:
#                        # since condition is a function, stencil must be so too. 
#                        # So, convert self to varcoeff-version. This will change the contribCoeff
#                        # we got abobe, so lets re-set this variable.
#                        self.convert_to_varcoeff()
#                        contribCoeff = self.stencil_coeff[key]
#                     self.addSource([_U.func_mulFuncFunc(_U.func_mulScalarFunclist(leftright[i]*width*grid.division[edgei],contribCoeff),condition)])
#                  else:
#                     # Condition is numeric
#                     if self.varcoeff:
#                        # convert condition to function:
#                        conditionfunc=constfunc(condition)
#                        self.addSource([_U.func_mulFuncFunc(_U.func_mulScalarFunclist(leftright[i]*width*grid.division[edgei],contribCoeff),conditionfunc)])
#                     else:
#                        self.addSource(leftright[i]*contribCoeff*width*grid.division[edgei]*condition)
#                        
#                  #TODO: Remove this debug print
#                  #print "added source %s, composed as this:" % (self.source)
#                  #print "contribCoeff: %e, width: %e, dx: %e, condition: %e" % (contribCoeff, width, grid.division[edgei],condition)
#
#                  # remove the original node.
#                  self.stencil_coeff.__delitem__(key)
#                  # No need to do anything more with this key:
#                  # TODO: If the key is "ghost in more than one direction" we may 
#                  # have just mapped a ghost into another. Say (-1,-1) => (1,-1), 
#                  # i.e., still a ghost... We can flag such things, and go through the
#                  # list again. But probably more care is needed, as the correct 
#                  # approximation to such nodes (as (-1,-1)) may need more threatment.
#                  i = len(edgeinds)
#                  print "Key is removed, go on to next key"
#               else:
#                  i += 1
#            else:
#               i += 1
       
   def convert_to_varcoeff(self):
      """
      Convert the current stencil from a fixed coefficient version to a 
      variable coefficient one"""
      
      if not self.varcoeff:
         self.varcoeff=1
         val=self.source
         self.source = [constfunc(val)]


         for key in self.stencil_coeff:
            val = self.stencil_coeff[key]
            self.stencil_coeff[key] = [constfunc(val)]

         # We must also set the function pointers to the right instances
         self.addNode = self.addNode_vc
         self.addSource = self.addSource_vc
         self.call = self.call_vc
         self.effect = self.effect_vc
         self.effect_source = self.effect_source_vc

   def stencilWidth(self):
       """Compute the widht of this stencil
          Return: a tuple of tuples, one tuple for each direction, containg two 
          values, (l,r), where l is the maximum negative offset in that direction 
          and r is the maximum positive offset in that direction.
       """

       # when we create the stencilWidth, we have to use list, since tuples 
       # are imutable
       l = [0 for x in xrange(self.nsd)]
       r = [0 for x in xrange(self.nsd)]
       for k in self.stencil_coeff:
           for i in xrange(self.nsd):
               l[i] = min(l[i],k[i])
               r[i] = max(r[i],k[i])
       return (tuple(l),tuple(r))

# generator for constant function-objects (easy way of making closures).
# We might also do this with closures (maybe new consept in python 2.2/2.3:
# def make_constfu(value):
#    def constfu(*args): return value
#    return constfu
#
# For now, we leave the original construction with class.
class constfunc:
   def __init__(self,value):
      self.val = value

   def __call__(self, *args):
      return self.val

