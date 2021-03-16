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
# Grid; A simple grid class 
#

from numarray import array

from copy import copy, deepcopy
import Field as _f
import Utils as _u

# Import numarray/Numeric
from Utils import _n
# this should probably be removed!
array = _n.array

import operator
import math

GridError = "Error in Grid"

class Grid:
    """
      Some Info on this class:
      datamembers include:
      nsd:          number of space dimensions
      geometry:     a list with geometry data, length = nsd
      sizespec:     a list of length nsd with partition data
      spec:         a unified version of sizespec (always len(spec) == nsd)
      division:     the deltas 
      fields:       Pointers to Fields attached to this grid. A list.
      isParallel:   Boolean flag, true if we are running on a parallel machine

      in case we are parallel, several other datastructures are available:
     
      communicationPattern:  The communication requirements in each 
                             direction, aggregates upward on downward communication
           
      partitionStrategy:     How we partition the grid in the differnt spatial directions
      nodePartition:         How nodes are distributed in each direction
      myGridPosition:        Where "I" am in the grid lattice
      Pspec:                 Same as sizespec, but for the subprocess
      Pgeometry:             Same as geometry, but for the subprocess
      Gspec:                 Store the global sizespec when going parallel
      Ggeometry:             Store the global geometry when going parallel
      Ispec:                 The interior sizespec, excluding ghosts
      Igeometry:             The interior geometry, exluding ghosts
      InodeStart:            Index for first internal node in each direction
      neighbours:            Neighbours I need to send something to, hash where
                             key is the offset tuple, value is the real gridposition
                             of the neighbour.
      neighbourCommPattern:  What I need to send to my neighbours, hash where key
                             is the offset tuple, value is a list with the 
                             number of nodes to send.
      receivers:             Neighbours I need to receive from, hash as in 'neighbours'
      neighbourRecvPattern:  Similar to neighbourCommPattern
      myGlobalNodalBounds:   The global/real nodalbounds for this subprocess
      myGlobalNodalBoundsIncGhosts: Same as above, but include ghosts.
    
    """

   
    def __init__(self,div=None,domain=([0.0,1.0],[0.0,1.0]),division=None):

        if not div and division:
            div = division

        if isinstance(div,(int,float)):
            # div is just one number - convert to tuple:
            div = (div,)
            # in that case - check length of domain. If it is 
            # larger than 1, we assume that domain is given as a 
            # simple tuple/list, i.e., domain=[0,1] instead of the 
            # somehow more correct domain=([0,1],) in our context...
            if len(domain) > 1:
                domain=(domain,)

        # aliases for backward compatibility.
        sizespec = div
        geometry = domain
        
        self.nsd = len(geometry)
        # We will get items from the field, not the grid
        #self.getitem = self.getitemUnset
        #self.setitem = self.setitemUnset

        # make sure that geometry is given with floats, and attach:
        float_geo = []
        for i in xrange(self.nsd):
            geo_dir = list(geometry[i])
            for j in xrange(2):
                if isinstance(geo_dir[j],type(1)):
                    geo_dir[j] = float(geo_dir[j])
            float_geo.append(geo_dir)
        self.geometry = tuple(float_geo)
        
        if sizespec:
            self.setSize(sizespec)
        
        # initialize the fields list
        self.fields = []

        # for later use in case of parallel solver:
        self.isParallel = False

        # In case I am parallel, store rank for process, else -1
        if _u.parallelLoaded:
            self.myrank = _u.pypar.rank()
        else:
            self.myrank = -1

        

        # for consistency in names, the name of getCorners should be just corners?
        # Maybe find some other way to fix backward compat.
        self.corners = self.getCorners

        # want to store a basenode that take parallelism into account;
        self.basenode = _u.getZeroTuple(range(self.nsd))
        self.base = self.getBase()


    def getGlobalNodalBounds(self, withGhosts=False):
        """If this is a parallel grid, compute the global nodal bounds for this
           subgrid, and return that as ((start),(stop)), where both start and stop 
           are tuples of length nsd.
           
           Usually, this method is called internally on changes to the partitions-scheme.
           External accessors shouw preferable use the pre-computed, stored 
           self.myGlobalNodalBounds

           If the method is called with the flag withGhosts=True, the global
           nodal bounds returned include the ghost nodes.
        """


        start_bounds = []
        stop_bounds = []
        if withGhosts:
            # compute nodal bounds with ghosts
            if self.isParallel:
                for i in xrange(self.nsd):
                    start_bounds.append(0)
                    stop_bounds.append(-1) # node numbering start on 0 in each direction!
                    for j in xrange(self.myGridPosition[i]):
                        start_bounds[i] += self.nodePartition[i][j]
                    stop_bounds[i] += start_bounds[i] + self.nodePartition[i][self.myGridPosition[i]]
                # adjust for ghosts:
                for i in xrange(self.nsd):
                    # subtract the number in InodeStart, which describe how far into the
                    # domain the real nodes starts
                    start_bounds[i] -= self.InodeStart[i]
                    # next, we need to adjust stop_bounds  such that we agree with 
                    # sizespec
                    stop_bounds[i] -= (stop_bounds[i] - self.sizespec[i] - start_bounds[i])
        else:
            if self.isParallel:
                for i in xrange(self.nsd):
                    start_bounds.append(0)
                    stop_bounds.append(-1) # node numbering start on 0 in each direction!
                    for j in xrange(self.myGridPosition[i]):
                        start_bounds[i] += self.nodePartition[i][j]
                    stop_bounds[i] += start_bounds[i] + self.nodePartition[i][self.myGridPosition[i]]
        return tuple((start_bounds,stop_bounds))


    def partition(self, stencilList = None, partitionStrategy = None, stencilWidth = None):
        """
        This method signal that we want to utilize a parallel computer. If the pypar 
        environment is loaded, that is ok, and we can try to accomplish that, else we 
        just continue in a scalar fashion.
        Parameters: 
            stencilList: if not None, we assume that the user supply a stencillist
               which should be used when deciding on an optimal partitioning. 
            partitionStrategy: if stencillist is None, we check this argument for 
               a tuple. Starting from the first index in this tuple, we partition 
               each dimension in grid accordingly, until either all dimensions are 
               partitioned or there is no more values in the tuple. We also 
               check that the number of parts the grid is divided into agrees with 
               the number of processors available.
           stencilWidth: Since the width of the stencils in the stencilList is 
               crucial in the process of dividing the grid and defining the ghost
               nodes and the communication patterns, we also need to supply a 
               "virtual" stencilWidth pattern together with the partitionstrategy, in
               case a real stencillist is not supplied. Remark that this is a 
               required parameter when stencillist is None and partitionstrategy is 
               given.
        Returns:
            Various error codes may be returned, in the form of integers:
                1: No parallel environment available, no partition done.
        """

        ### TODO: We do not consider the where-parts in  Stencils in the list
        ### probably can the communication and storage of ghosts be optimized
        ### even further by considering this. From an algorithmic point of view
        ### it doesn't hurt to deal with the global stencilListWidth

        # define some variables to increase readability:
        left = 0
        right = 1

        # if we do not have a parallel environment loaded, it is no reason to 
        # continue:
        if not _u.parallelLoaded:
            # return with error-code 1
            return 1

        if self.nsd > 3:
            raise GridError, "Only 1, 2, and 3 space dimensions are supported by the parallell distribution algorithm"

        # we need to know the number or processors available to us:
        Psize = _u.pypar.size()
        # also handy to know who I am:
        Prank = _u.pypar.rank()

        # check that we have more than 1 cpu, if not stop partitioning.
        if not Psize > 1:
            return 2

        self.communicationPattern = [0 for x in xrange(self.nsd)]
        if stencilList:
            stencilWidth = stencilList.stencilListWidth()


            for i in xrange(self.nsd):
                self.communicationPattern[i] = stencilWidth[right][i] - stencilWidth[left][i]

            partitionStrategy = _u.createPartition(Psize,self.communicationPattern,self.sizespec)
            # save in self for later reference.
            self.partitionStrategy = partitionStrategy

        elif partitionStrategy and stencilWidth:
            # both partitionstrategy and stencilWidth are given as input
            self.partitionStrategy = partitionStrategy
            raise GridError, "Partitioning without stencilList is not implemented"

        else:
            raise GridError, "Either stencilList or both partitionStrategy and stencilWidth must be given for partitioning" 

        # compute the distribution of nodes on processors
        # remember that spec != nodes (there is one more node!)
        # Is that a clever or common choice? i.e., that a user specify the 
        # division into "segments" and not division into nodes (counting the ends)?

        nodePartition = _u.distributeNodes(_u.specToNodes(self.sizespec),partitionStrategy)

        self.nodePartition = nodePartition

        ### REMARK: From here I guess we go mostly parallel, i.e., everyone are not
        ###         doing the same thing anymore, at least not completely the same

        # it's about time to locate who I am in the system
        mgp = _u.rankToProc(Prank,partitionStrategy)
        self.myGridPosition = mgp

        # remember that nodePartition is given in terms of nodes, convert back 
        # to spec, that is, subtract one in each direction:
        Pspec = []
        self.Pspec = Pspec
        Pgeometry = []
        self.Pgeometry = Pgeometry
        for i in xrange(self.nsd):
            Pspec.append(nodePartition[i][mgp[i]]-1)
            tmpstart = self.geometry[i][left]
            for j in xrange(mgp[i]):
                tmpstart += nodePartition[i][j]*self.division[i]
            Pgeometry.append([tmpstart,tmpstart + Pspec[i]*self.division[i]])


        # we will now replace the original geometry and sizespec with parallel 
        # versjons. The original (global) versions are stored in G{geometry,spec}
        self.Ggeometry = self.geometry
        self.Gspec = self.sizespec
        self.geometry = self.Pgeometry
        self.sizespec = self.Pspec
        # also save "internal" datastructures, as we may extend the geometry/sizespec
        # with ghost nodes later on.
        self.Igeometry = deepcopy(self.geometry)
        self.Ispec = deepcopy(self.sizespec)
        # before we add ghosts, Internal nodes start on index 0 in each direction
        self.InodeStart = map(lambda x: 0,self.Ispec)

        # We need some concept of differntiate between physical (real) boundaries 
        # and internal "parallel" boundaries.
        physBound = []
        for i in xrange(self.nsd):
            tmpl = tmpr = False

            if mgp[i] == 0:
                # we are on the "left" side
                tmpl = True
            if mgp[i] == partitionStrategy[i]-1:
                # we are on the "right" side
                tmpr = True
            physBound.append([tmpl,tmpr])
            

        # We will now use the stencilWidth information to add ghost-nodes for 
        # those parts of the boundary which is not physical boundary
        for i in xrange(self.nsd):
            if not physBound[i][left]:
                # add the left (negative, hence add!) offset
                self.geometry[i][left] += self.division[i]*stencilWidth[left][i]
                # subtract the same (still negative) number of nodes for left offset
                self.sizespec[i] -= stencilWidth[left][i]
                # we also need to move the InodeStart the same number of nodes
                self.InodeStart[i] -= stencilWidth[left][i]
            if not physBound[i][right]:
                self.geometry[i][right] += self.division[i]*stencilWidth[right][i]
                self.sizespec[i] += stencilWidth[right][i]


        # Finally, we build a datastructure which know what I have to communicate,
        # and to whom I should communicate this.
        # communication is to be interpreted as slices to be sent, from the "real" 
        # nodes in the subgrid, that is, not ghost nodes, but counting cells from the 
        # virtual boundary and inwards.
        # Remark that unlike the communicationPattern, which aggregate the communication
        # in each direction, we here store the bare necessities of communication
        # in each direction, both upwards and downwards.

        self.neighbours = self.commNeighbours(stencilList)
        
        # for each neighbours we communicate with, decide what to communicate
        self.neighbourCommPattern = {}
        for n in self.neighbours:
            pattern = []
            for i in xrange(self.nsd):
                if n[i] == 0:
                    # communicate all Internal nodes in this direction
                    pattern.append(self.Ispec[i]+1)
                elif n[i] > 0:
                    # communicate upwards, i.e., downward part of the stencil:
                    pattern.append(stencilWidth[left][i])
                elif n[i] < 0:
                    # communicate downwards, i.e, upward part of the stencil:
                    pattern.append(stencilWidth[right][i])
            self.neighbourCommPattern[n] = pattern


        # we also need to store info about receiving data
        self.receivers = self.recvNeighbours(stencilList)
        self.neighbourRecvPattern = {}
        # for each neighbour I receive from, decide which ghost-nodes to fill
        # TODO: Check whether we always only fill into ghost, or could it 
        # happend that there should be some averageing between what I already have
        # and what I got (anyway, this is probably decisions which belongs in the 
        # Field).
        for n in self.receivers:
            pattern = []
            for i in xrange(self.nsd):
                if n[i] == 0:
                    # be prepared to reveive a full slice in this direction
                    pattern.append(self.Ispec[i]+1)
                elif n[i] > 0:
                    # receive from above:
                    pattern.append(stencilWidth[right][i])
                elif n[i] < 0:
                    # receive from below:
                    pattern.append(stencilWidth[left][i])
            self.neighbourRecvPattern[n] = pattern

        # Finish up by setting parallel to True:
        self.isParallel = True

        # and then compute my global bounds:
        self.myGlobalNodalBounds = self.getGlobalNodalBounds()
        self.myGlobalNodalBoundsIncGhosts = self.getGlobalNodalBounds(withGhosts=True)
        # update values for base
        self.basenode = self.myGlobalNodalBoundsIncGhosts[0]
        self.base = self.getBase()

        # If any Field is connected to the Grid from before we go into parallel, 
        # they should now be notified on the new status and change accordingly.

        for f in self.fields:
            f.doInitParallel()
        
        # The same must be done for all the where-parts in the stencilList in use 
        # here, as they now should reflect the real nodes in this subpart, not the
        # global nodes.

        if stencilList:
            stencilList.doInitParallel()

        # return
        return 0

    def commNeighbours(self, stencilList):
        """Based on the signature of the stencilList, decide which neighbours
           I need to communicate with
           This method assume that self.myGridPosition is assigned, as well as 
           self.partitionStrategy
        """

        mgp = self.myGridPosition
        ps = self.partitionStrategy

        neighbours = {}
        for s in stencilList.signature():
            # the center is also part of the stencil signature, but can not
            # result in any communication. Hence, make sure that the center 
            # is removed. The signature of the center is all zeros.
            if reduce(operator.add,map(lambda x: abs(x),s)) == 0:
                # this is the center
                continue
            s_l = list(s)
            for i in xrange(len(s_l)):
                if abs(s_l[i]) > 1:
                    s_l[i] = s_l[i]/abs(s_l[i])
                # turn around to find the neighbour who need this
                s_l[i] *= -1
            s_t = tuple(s_l)
            if not neighbours.has_key(s_t):
                # need to compute real pid of neighbour:
                real_n = _u.tupleAdd(mgp,s_t)
                # then we have to check that real_n is within ps
                flag = True
                for i in xrange(len(real_n)):
                    if real_n[i] >= ps[i] or real_n[i] < 0:
                        flag = False
                if flag:
                    neighbours[s_t] = real_n
        return neighbours


    def recvNeighbours(self, stencilList):
        """Based on the signature of the stencilList, decide which neighbours
           I will receive from.
           This method assume that self.myGridPosition is assigned, as well as
           self.partitionStrategy
        """

        mgp = self.myGridPosition
        ps = self.partitionStrategy

        receivers = {}
        for s in stencilList.signature():
            # I can't receive anything from myself, hence remove the center.
            if reduce(operator.add,map(lambda x: abs(x),s)) == 0:
                # this is the center
                continue
            s_l = list(s)
            for i in xrange(len(s_l)):
                if abs(s_l[i]) > 1:
                    # normalize
                    s_l[i] = s_l[i]/abs(s_l[i])
            # convert to tuple as only immutables can be used as index in a hash
            s_t = tuple(s_l)
            if not receivers.has_key(s_t):
                # need to compute real pid of neighbour:
                real_n = _u.tupleAdd(mgp,s_t)
                # check that real_n is within ps, to avoid communicating
                # of physical boundaries
                flag = True
                for i in xrange(len(real_n)):
                    if real_n[i] >= ps[i] or real_n[i] < 0:
                        flag = False
                if flag:
                    receivers[s_t] = real_n
        return receivers


    def setSize(self,sizespec):
        """Set the division for this grid. 
         The sizespec argument denotes the max index in each direction (also) 
         Initialize the Field datastructure. If
         fewer elements than nsd is given, the first are used for all
         directions. If length of sizespec is the same as nsd, the data will
         have sizes accordingly in each direction. If length of sizespec is
         larger than nsd, only the first nsd directions are used."""

        # sizespec given as one integer:
        if type(sizespec) == type(1):
            spec = _n.zeros(self.nsd) + sizespec
        # sizespec given as array of size 1:
        elif len(sizespec) < self.nsd:
            spec = _n.zeros(self.nsd) + sizespec[0]
        # sizespec given as indices set:
        elif len(sizespec) == self.nsd:
            spec = sizespec
        # sizespec given as too long indices set:
        else:
            spec = sizespec[:self.nsd]


        # Set the objects sizespec, and the division parameters (dx,dy, etc.)
        self.sizespec = spec 

        self.division = []
        divvars = ['dx','dy','dz']
        for i in xrange(self.nsd):
            g = self.geometry[i]
            self.division.append ((g[1] - g[0])/self.sizespec[i])
            if i < len(divvars):
                self.__dict__[divvars[i]] = self.division[i]

    def getCorners(self):
        """getCorners: return a list of corner-tuples"""

#        n_max=array(self.sizespec)
#        startup = _u.getZeroTuple(self.sizespec)
#        stoptup = _u.getOnesTuple(self.sizespec)
#        c_i = _u.tupleIterator(self.nsd,startup,stoptup)
#        return [tuple(n_max*array(i)) for i in c_i]
#
#        Reimplemented to use a proper iterator with reset and support for changing
#        bounds

        return _u.cornerTupleIterator(self.nsd,_u.getZeroTuple(self.sizespec),self.sizespec)

    def getBase(self):
        return map(lambda x: x[0],self.geometry)

    def indexToPoint(self,index):
        """Accept a index-tuple as argument"""
        
#        base = self.getBase()
#        gnb = _u.getZeroTuple(base)
#        if self.isParallel:
#            # index-adjuster according to global-nodalbounds:
#            gnb = self.myGlobalNodalBoundsIncGhosts[0]
#        for i in xrange(len(base)):
#            base[i] += (index[i]-gnb[i])*self.division[i]
#        return base

        #globindex = tuple(map(lambda x,y:x-y,index,self.basenode))
        globindex = tuple(map(operator.sub,index,self.basenode))
        #
        # This map should be replace with the list-comprehension:
        #
        # [b+x*y for b,x,y in zip(self.base,globindex,self.division)]
        #
        # as map may disappear in the future. The list-comprehension/zip is a
        # bit slower, so I may just leave the map for now...
        #
        return map(lambda b,x,y: b+x*y,self.base,globindex ,self.division)
        

    def indicesToPoints(self, i):
        """Accept a list of index-tuples i as argument, and convert to points"""

        globinds = map(lambda index: tuple(map(operator.sub,index,self.basenode)),i)

        return map(lambda globindex: map(lambda b,x,y: b+x*y, self.base,globindex, self.division),globinds)


    def allPoints(self):
        """Get an iterator for all points in this grid"""

        iter = _u.tupleIterator(self.nsd,_u.getZeroTuple(self.sizespec),self.sizespec)
        return iter

    def innerPoints(self, region=None):
        """Get an iterator for all inner points on this grid"""
        
        iter = _u.innerTupleIterator(self.nsd,_u.getZeroTuple(self.sizespec),self.sizespec)
        if region:
            # Restrict the iterator to the given region:
            bounds = self.regionToBounds(region)
            iter.setNewBounds(bounds,iOB=True)
        return iter

    def boundary(self, region=None, type='box',direction='in',radius=None,center=None):
        """Get an iterator for all boundary points for this grid"""
        
        iter = _u.boundaryTupleIterator(self.nsd,_u.getZeroTuple(self.sizespec),self.sizespec)
        if region:
            bounds = self.regionToBounds(region)
            iter.setNewBounds(bounds,iOB=True)

        if type == 'box':
            return iter
        elif type == 'circle':
            if not center:
                raise GridError, "center must be specified for circelIterator at boundary"
            if not radius:
                raise GridError, "radius must be specified for circelIterator at boundary"
            c_iter = _u.circleIterator(iter,self,center,radius,direction)
            return c_iter
        else:
            raise GridError, "Unknown boundary-type, %s" %(type)

    def regionToBounds(self,region):
        # Given a (coordinate-based) region, convert to a (index-based) bounds
        region_conv = [[r[0] for r in region],[r[1] for r in region]]
        d = self.division
        l = self.base
        start = [int(math.ceil((R - L)/D)) for (R,L,D) in zip(region_conv[0],l,d)]
        stop = [int(math.ceil((R - L)/D)) for (R,L,D) in zip(region_conv[1],l,d)]

#        start = []
#        stop = []
#        for i in xrange(len(region)):
#            d = self.division[i]
#            l = self.base[i]
#            start.append(int(math.ceil((region[i][0] - l)/d)))
#            stop.append(int(math.floor((region[i][1] - l)/d)))
        bounds= (start,stop)
        return bounds

    def range(self):
        """Return range for grid as (x0,y0,..),(xN,yM,...)
           A method which return the min and max corners in the grid.
        """
        range = []
        # all grids has its origin in (0,0,...)
        range.append(_u.getZeroTuple(self.sizespec))
        # then find the max-corner, which is also simple...
        range.append(self.sizespec)
        return range


    def attachField(self, f):
        """Attach a Field to the list of fields in this grid instance"""
        self.fields.append(f)
