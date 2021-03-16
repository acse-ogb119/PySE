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
# Field; A simple field class
#


import Grid as _g

import Utils as _u

# import Numeric/numarray!
from Utils import _n

import pyPDE.VizSurface
import pyPDE.VizCurve
from pyPDE import seq

from operator import mul, add
from math import sqrt
from copy import copy

from pyFDM.pyfdmcutils import closestint
import time

FieldError = "Error in Field"

class Field(object):
    """ 
         A simple Field class. 
         datamembers include:
         data:          A NumPy array, holds the actual data
         shapedata:     reference to 'data', shaped according to Grid.
         grid:          A field is data over a Grid. Here we have the Grid. The
                        Grid must be known on creation.

         datasize:      Size of the dataarray, given as one integer (the actual
                        datastorage is 1D)

         indexspec:     a list used to compute the global data index from a indices
                        set

         sizespec:      The max-indices in each direction, taken from the
                        underlaying grid.

         communicator:  Store slice-objects for data we want to communicate. A hash 
                        where the key is the rank of the receiver. Intended as internal
                        datastructure
         
         receiviator:   Store slice-objects for data we need to receive in order to 
                        update ghost nodes. A hash where the key is the rank of the 
                        sender. Intended as internal datastructure.
                        
         isParallel:    Boolean flag, true if we are running on a paralllel machine.

         myrank:        If running on a parallel machine, set to the mpi-rank of 
                        the current process, else -1.
    """

    def __init__(self, grid, data=0, comm=None, recv=None):
        """Initialize Field with a grid.
           Data can optionally be added. This should be a NumPy array of matching 
           size. No check is done on either size or content. 
           Normally, data should only be specified when making copies of a Field
           Other Arguments:
           comm: A parallel communicator, containing all slices we need to 
                 communicate with neighbours. A hash. 
           recv: A parallel "receiviator" (made-up word), containing all slices 
                 we should receive from neighbours. A hash.

                 These are an internal datastructures, but in case of a already 
                 parallel Field is copied or used as basis for a new Field, we can 
                 attach the communicator and receiviator to avoid recomputation.
        """
       
        # If we are parallel, store my rank as it is become useful later on.
        if _u.parallelLoaded:
            self.myrank = _u.pypar.rank()
        else:
            self.myrank = -1

        #TODO: Remove debug print
        #print self.myrank, "> Creating a new field ",id(self)

        self.grid = grid

        # If we go parallel later on, it is important that grid know me. So 
        # we create a relation from grid to field as well.
        self.grid.attachField(self)

        # use the simple indexToDataIndex mapping as long as we are scalar
        self.indexToDataIndex = self.localIndexToDataIndex

        self.initParallel = False
        self.isParallel = False
        if grid.isParallel:
            self.initParallel = True
            self.indexToDataIndex = self.globIndexToDataIndex

        #Get size-spec. from grid

        try:
            self.sizespec = self.grid.sizespec
        except:
            raise FieldError, "division of underlaying grid is not set"

        # default, we use origin as base-node for indexing.
        self.basenode = _u.getZeroTuple(self.grid.sizespec)
        # TODO: Remove debug print
        #print self.myrank,"> Set basenode in __init__: ",self.basenode,", me is ",id(self)

        # compute the datasize, and store it.
        #
        # this can be done much simpler (but it might not be necessary at all)
#        datasize = 1
#        for s in self.sizespec:
#            s += 1
#            datasize *= s
#        self.datasize = datasize

        # initialize the data
        self.initData(data)

#        #TODO: Remove debug print
#        if self.initParallel:
#            print self.myrank, "> initParallel is true, ",id(self)
#        else:
#            print self.myrank, "> initParallel is false, ",id(self)
#        #TODO: Remove debug print
#        if comm:
#            print self.myrank, "> comm is set, ",id(self)
#        else:
#            print self.myrank, "> comm is not set, ",id(self)

        if self.initParallel and not comm:
            # no need to check for recv, as comm and recv are computed the same 
            # place (at least for now), so it will be recomputed anyway.
            #TODO: Remove debug print
            #print self.myrank, "> do initParallel as requested, ",id(self)
            self.doInitParallel(startedParallel=True)
            self.isParallel = True
        elif comm and recv:
            # If a communicator and a receiviator are given, we assume that we 
            # started parallel, and that there is no other parallel initialization 
            # we have to do
            self.communicator = comm
            self.receiviator = recv
            self.isParallel = True
            # we have to reset the basenode though!
            self.basenode = self.grid.myGlobalNodalBoundsIncGhosts[0] 
            # TODO: Remove debug print
            #print self.myrank,"> Set basenode again in __init__, for parallel: ",self.basenode
            self.lns = self.localnodes_slice()


    def initData(self,data=0):
        self.datashape = map(lambda x: x+1,self.sizespec)
        self.datasize = reduce(mul, self.datashape)
        # initialize data:
        # If any data attached (preferable numarray, have to find a better way 
        # to test), assume that this have the right size and attach.
        # workaround here because data as numarray is not accepted in test.
        if not isinstance(data,int):
            self.data = data
        else:
            self.data = _n.zeros(self.datasize,_n.Float)
        self.shapedata = self.getShapedData()

        # indexspec 
        maph = MapHelper(1)
        self.indexspec = map(maph.multiply,[1] + map (lambda v: (v+1),self.sizespec[:-1]))

        # make sure that we have a communicator. 
        self.communicator = {}
        # probably need to make sure that I have one for receiving also
        self.receiviator = {}

        #create necessary data for mapping from shape-index to 1d-index:
        sh = self.datashape[1:]
        sh.append(1)
        sh.reverse()
        shapemul = []
        tmp = 1
        for i in xrange(len(sh)):
            tmp *= sh[i]
            shapemul.append(tmp)
        shapemul.reverse()
        self.shapemultiplicator = shapemul

    def getShapedData(self):
        return _n.reshape(self.data,self.datashape)

    def doInitParallel(self,startedParallel=False):
        """Do necessary initialization related to a parallel grid
           This method can be called in two different ways:
           1. Field is created from a non-parallel grid, grid is then parallized
              later on, and this method is called to change all attahced Fields in
              the grid.  In this case, startedParallel should be False
           2. Grid is already parallel when a Field is created over the Grid. In 
              this case this method is run directly from Field.__inti__, and 
              startedParallel = True
           
           Data in the field is preserved during partitioning.
        """


        # TODO: Check this code. It seems like most of this code is Grid-centric,
        # and may actually belong there and not here. This fact may be even more 
        # important if we switch to some other Grid/Mesh tool.
        # At leaset, a grid-based numbering is used, so I doubt this belongs here.
        #
        # but I am not su sure any more, as the slicing clearly belongs here
        # in the field.
        #
        
        mgp = self.grid.myGridPosition

        if not startedParallel:
            # Have to change datasize, sizespec and so forth (reattach...)
            self.sizespec = self.grid.sizespec # point to the new onw
            # re-initialize data
            # In case there is values in data, we should copy out the right slice of data
            # and reuse.
            # TODO: Check whether this is something we should do only on request (i.e.
            # copy data - if this is creation of a new Field, there is nothing to 
            # copy)

            globshape = tuple(_u.tupleAdd(self.grid.Gspec,_u.getOnesTuple(self.grid.Gspec)))
            shapedata = _n.reshape(self.data,globshape)
            myglobnodes = self.grid.myGlobalNodalBoundsIncGhosts
            # Endpoint is not included in slices, so I need to increase end with one.
            dataslices=[]
            # store the first element in myglobnodes as the basenode, 
            # i.e global indexing will be relative to this.
            self.basenode = tuple(myglobnodes[0])
            # TODO: Remove debug print
            #print self.myrank,"> Set basenode in doInitParallel: ",self.basenode,", me is ",id(self)
            for i in xrange(self.grid.nsd):
                start = myglobnodes[0][i]
                stop = myglobnodes[1][i]+1
                dataslices.append(slice(start,stop))
            # TODO: Check again whether we should have a .copy() here:
            savedata=shapedata[dataslices].flat
            self.initData(savedata)
        else:
            # we still need to set the basenode...
            self.basenode = tuple(self.grid.myGlobalNodalBoundsIncGhosts[0])
        
        # need to build the necessary communication "slices" 
        # use the grid.neighbourCommPattern to build these slices
        for n in self.grid.neighbourCommPattern:
            # compute the rank of the node who get this
            # n is just offset, not the actual position in the partition-grid
            # The actual position is stored in grid.neighbours:
            r = _u.procToRank(self.grid.neighbours[n],self.grid.partitionStrategy)
            # build the slice we want
            commNodes = self.grid.neighbourCommPattern[n]
            slices = []
            for i in xrange(self.grid.nsd):
                if commNodes[i] >= 0:
                    # if n[i] == 0, the interior nodes should be chosen, as 
                    # this is not the "ghost" direction
                    if n[i] == 0:
                        start = self.grid.InodeStart[i]
                        stop = start + commNodes[i]
                        slices.append(slice(start,stop))
                    else:
                        # start from bottom, slice from InodeStart[i] to the given number 
                        start = self.grid.InodeStart[i]
                        stop = start + commNodes[i]
                        slices.append(slice(start,stop))
                if commNodes[i] < 0:
                    # start from the back, and step this many nodes backwards
                    start = self.grid.InodeStart[i] + self.grid.Ispec[i]
                    # add commNodes to start, since we have negative number
                    stop = start + commNodes[i]
                    # create slice with backwards singlestepping
                    slices.append(slice(start,stop,-1))
            self.communicator[r] = tuple(slices)

        # need to build the necessary receiving "slices" - from whom I 
        # receive, and where should I put what i get.
        for n in self.grid.neighbourRecvPattern:
            # compute the rank of the process I receive from
            r = _u.procToRank(self.grid.receivers[n],self.grid.partitionStrategy)
            recvNodes = self.grid.neighbourRecvPattern[n]
            slices = []
            for i in xrange(self.grid.nsd):
                if recvNodes[i] >= 0:
                    # TODO: Figure out how to decide when we have a full slice
                    # and when we just receive according to the width of the 
                    # stencil in this direction.
                    # if n[i] == 0, we know that we are in the full-slice setting, 
                    # and we fill "internal" nodes in this direction - other parts
                    # of the slice will make sure that we fill ghosts.
                    if n[i] == 0:
                        start = self.grid.InodeStart[i]
                        stop = start + recvNodes[i]
                        slices.append(slice(start,stop))
                    else:
                        # we fill the "upper" ghosts
                        start = self.sizespec[i]
                        stop = start - recvNodes[i]
                        slices.append(slice(start,stop,-1))
                if recvNodes[i] < 0:
                    # receive from below, put into the first nodes:
                    start = 0
                    stop = -recvNodes[i]
                    slices.append(slice(start,stop))
            self.receiviator[r] = tuple(slices)

        # I am parllelized, so tell the world:
        self.isParallel = True
        self.lns = self.localnodes_slice()

        # Switch to the version of indexToDataIndex that maps the index from 
        # global to local
        self.indexToDataIndex = self.globIndexToDataIndex



    def updateField(self):
        """In case of a parallel Field, update ghost-nodes on each process
           to reflect the newest "real" values on neighbouring processes
        """
        
        if not self.isParallel:
            return None

        # update the Field, that is, communicate the slices 
        #
        # we assume that this method is called at the same time on 
        # all processors. Is that fair?

        # communicate based on rank: 
        #  - If myrank is 0, I just start sending to my peers
        #  - If myrank > 0 and I am suppose to receive from someone 
        #    with a lower rank than myself, I start listening for those
        #  - After finishing reception, I send what I need to send to those 
        #    with higher rank than me.
        #  - Next, I send to those with lower rank than me
        #  - Because we to all the sending at once, we may question to split this,
        #    but I need to look into how performance is affected if we just do all 
        #    the sending at once - maybe it is even better done that way? 
        #    TODO: chekc that!
        #  - Finally, I recevie what I should receive from those with higher
        #    rank than me.
        #
        #  Remark: as process with myrank 0 can not receive from anyone with 
        #  lower rank, we do not need to do any special arrangement for myrank
        #  0, it will go right to the send-section.
        #  
        #  If we can use non-blocking send/recv, this can be programmed a bit 
        #  more "parallel" - i.e., start all receivers non-blocking, then start
        #  sending, and at last, wait for finished reception. But pypar does not 
        #  suppoert non-blocking send/recv, so we leave that for now.

        # TODO: as of yet, pypar only support Numeric, so I have to copy data into 
        # a Numeric array, and receive in the same. Afterwards, I need to copy back 
        # into the original slices.
        import Numeric
        
        # pick r < myrank in receiviator:
        recvfirstkeys = filter(lambda x: x<self.myrank, self.receiviator.keys())
        for r in recvfirstkeys:
            # receive directly into a buffer - slice from shapedata.
            # TODO: fill directly into slice (hopefully possible...)
            #rdata = _u.pypar.receive(r, buffer=self.shapedata[self.receiviator[r]])

            # We let pypar create a receive-buffer, which will be of Numeric type
            rdata = _u.pypar.receive(r)
            self.shapedata[self.receiviator[r]] = rdata
            
        # pick r > myrank in communicator:
        sendfirstkeys = filter(lambda x: x>self.myrank, self.communicator.keys())
        for r in sendfirstkeys:
            # create a Numeric array for sending:
            tmpsend = Numeric.array(self.shapedata[self.communicator[r]])
            _u.pypar.send(tmpsend,r)

        # send to the rest, probably someone below is ready to listen now:
        sendrestkeys = filter(lambda x: x<self.myrank, self.communicator.keys())
        for r in sendrestkeys:
            # create a Numeric array for sending:
            tmpsend = Numeric.array(self.shapedata[self.communicator[r]])
            _u.pypar.send(tmpsend,r)

        # finally, we receive from those above:
        recvlastkeys = filter(lambda x: x>self.myrank, self.receiviator.keys())
        for r in recvlastkeys:
            # receive directly into a buffer - slice from shapedata.
            # TODO: fill directly into slice (hopefully possible...)
            #rdata = _u.pypar.receive(r, buffer=self.shapedata[self.receiviator[r]])

            # We let pypar create a receive-buffer, which will be of Numeric type
            rdata = _u.pypar.receive(r)
            self.shapedata[self.receiviator[r]] = rdata

    def __copy__(self):
        """Define how Field should respond to copy
           We want the reference to grid to be preserved, while we want new data
           to be allocated
        """

        # create a new field, with a copy of my data, reference to the same grid
        t = Field(self.grid, self.data.copy(), self.communicator.copy(), self.receiviator.copy())
        return t

    def copy(self):
        """convenience function, class implements __copy__"""
        return copy(self)

    def __iadd__(self, y):
        """Add the field y into this field"""
        if isinstance(y, Field):
            if id (y.grid) == id(self.grid):
                self.data += y.data
                return self
            else:
                raise FieldError, "Unable to iadd, grids does not match"
        else:
            raise FieldError, "Unable to iadd, operand is not a Field"

    def __add__(self, y):
        """Define how to add Fields together"""
        if isinstance(y, Field):
            if id(y.grid) == id(self.grid):
                # both are Fields, and have same grid.
                # use copy instead!
                #t = Field(self.grid, self.data.copy())
                t = self.copy()
                t.data += y.data
                return t
            else:
                raise FieldError, "Unable to add, grid does not match"
        else:
            raise FieldError, "Unable to add, operand is not a Field"

    def __mul__(self,val):
        """Implement multiplication with a scalar, when Field is the left operand
           If the right operand is another Field, a pointwise multiplication 
           is performed.
        """

        if isinstance(val, (int,float)):
            t = self.copy()
            t.data *= val
            return t
        elif isinstance(val, Field):
            t = self.copy()
            t.data *= val.data
            return t
        else:
            raise FieldError, "Unable to multiply with %s, unknown type %s" % (val,type(val).__name__)


    def __imul__(self, val):
        """Multiply the Field with a scalar, and modify the Field itself."""

        if isinstance(val, (int,float)):
            self.data *= val
            return self
        else:
            raise FieldError, "Only scalar supported in inplace multiplication"


    def __rmul__(self, val):
        """Define multiplication with a scalar, where self is the right 
           operand.
        """
        if isinstance(val, (int,float)):
            # use self.copy instead!
            #t = Field(self.grid, self.data.copy())
            t = self.copy()
            t.data *= val
            return t
        else:
            raise FieldError, "Unable to r-multiply with %s" % (val)

    def __sub__(self, y):
        """Define how to subtract Fields"""
        if isinstance(y, Field):
            if id(y.grid) == id(self.grid):
                # both are Fields and have same grid.
                # use self.copy instead!
                #t = Field(self.grid, self.data.copy())
                t = self.copy()
                t.data -= y.data
                return t
            else:
                raise FieldError, "Unable to mul, grid does not match"
        else:
            raise FieldError, "Unable to add, operand is not a Field"

    def norm(self, p=2):
        """Compute the norm of the field. 
           p=2 mean discrete L2-norm, other values of p is ok.
        """

        if p == 2:
            return sqrt(_n.dot(self.data,self.data))
        else:
            return reduce(add,self.data**p)**1./p
           
    def normL2(self, p=2):
        """Compute the norm of the field. 
           p=2 mean discrete L2-norm, other values of p is ok.
           weighted with grid divisions.
        """

        h = reduce(mul,self.grid.division)

        if p == 2:
            if self.isParallel:
                # pick right subset, i.e, exclude ghosts.
                import Numeric
                local_norm = Numeric.zeros(1,typecode='d')
                global_norm = Numeric.zeros(1,typecode='d')

                # compute local norm:
                localdata = self.shapedata[self.lns].flat
                local_norm[0] = _n.dot(localdata,localdata)

                # reduce and broadcast:
                _u.pypar.reduce(local_norm,_u.pypar.SUM,0,buffer=global_norm)
                _u.pypar.broadcast(global_norm,0)
                return sqrt(h*global_norm[0])
            else:
                return sqrt(h*_n.dot(self.data,self.data))
        else:
            if self.isParallel:
                raise FieldError, "Parallelversion of normL2 for p != 2 is not implemented"
            else:
                return (h*reduce(add,self.data**p))**1./p

    def supnorm(self):
        """ Compute the sup-norm, that is the largest absolute value 
        
            This function try to use scipy.amax/amin if available
        """

        try:
            import scipy
            datamax = scipy.amax(self.data)
            datamin = scipy.amin(self.data)
        except:
            datamax = max(self.data)
            datamin = min(self.data)

        return max(abs(datamax),abs(datamin))
           
    def inner(self, v):
        """Compute the inner product of self and v.
           Parameter: v: the right operand, must be another Field.
           Return:    integer: the innerproduct (self,v)

        """
        
        #TODO:
        #Remark: As pypar only support reduce of Numeric arrays, we 
        #        import Numeric here, and create a Numeric array of 
        #        length 1 for reduce of the innerproduct.
        #
        #        Also, there is no allreduce, so we have to broadcast the
        #        result after we are done compute it.

        if isinstance(v, Field):
            if self.isParallel:
                # need to pick the right subset from self.data.
                
                # build a slice list:
                sl = []
                for i in xrange(self.grid.nsd):
                    start = self.grid.InodeStart[i]
                    # remember: slice is up to, but not including, we need to add
                    # 1 'cause we use the spec
                    stop = start + self.grid.Ispec[i] + 1
                    tmpslice = slice(start,stop)
                    sl.append(tmpslice)

                # need to slice shapedata with this slice, and flatten it for 
                # use with dot:
                
                # for now, pypar can only handle Numeric correctly, so import 
                # Numeric, and create Numeric arrays of type double for local and 
                # global innerproduct:
                import Numeric
                local_ip = Numeric.zeros(1,typecode='d')
                global_ip = Numeric.zeros(1,typecode='d')

                # compute local innerproduct
                local_ip[0] = _n.dot(self.shapedata[sl].flat,v.shapedata[sl].flat)

                # next, we need to use reduce to compute the true innerproduct
                # We assume that all processers call inner at the same time (approx...)
                # We collect the sum on processor 0, using the operation SUM

                _u.pypar.reduce(local_ip,_u.pypar.SUM,0,buffer=global_ip)

                # finally, broadcast the value and return the scalar value:
                _u.pypar.broadcast(global_ip,0)

                return global_ip[0]
            else:
                return _n.dot(self.data,v.data)
        else:
            raise FieldError, "Unable to compute dot-product since right operand is not a Field"

    def __getitem__(self,index):
        """operator-method for field[x], where x is of length nsd.
            This method do not check index, it just assume that it correspond to
            a index in the grid, not in the 1D field indexing"""

#        # indexing is assumed to be relative to self.basenode, that is, we
#        # subtract the number in basenode[i] in each direction i.
#        globindex = 0
#        for i in range(len(index)):
#            globindex += (index[i]-self.basenode[i])*self.indexspec[i]
#
#        # TODO: Remove debug print (and computations)
#        maxindex = len(self.data)
#        if globindex >= maxindex or globindex < 0:
#            print "*************** start ****************"
#            print self.myrank,"> Get item: ",index,self.indexspec, self.data.shape
#            print self.myrank,"> Try to get: ",globindex
#            print self.myrank,"> my base: ",self.basenode
#            print self.myrank,"> my sizespec: ",self.sizespec
#            #print self.myrank,"> mgp: ",self.grid.myGridPosition
#            #print self.myrank,"> nodebounds: ",self.grid.myGlobalNodalBoundsIncGhosts
#            print "*************** end ****************"
#        if globindex >= maxindex or globindex < 0:
#            raise FieldError, "Try to get index out of range in field"
#        return self.data[globindex]

        # we assume that the shapedata is in shape, and can just get the 
        # value directly, after adjusting according to the basenode.
        #TODO: remove debug print
        #print self.myrank, "> request ",index," related to ",self.sizespec,self.basenode,self.indexspec
        globindex = tuple(map(lambda x,y:x-y,index,self.basenode))
        #print self.myrank, "> converted to ",globindex
        #print "this is what I return ",self.shapedata[globindex]
        return self.shapedata[globindex]

    def getValByPoint(self,*point):
        """Given a point, find an appropriate field value for the point. 
           The closest index will be used.
        """
#
#        theindex = [closestint((i-j)/d) for (i,j,d) in zip(point,self.grid.base, self.grid.division)]
#        try:
#            thevalue = self[theindex]
#        except:
#            print self.myrank,">Want to get value for point: ",point
#            print self.myrank,">grid: ",self.grid.division,self.grid.geometry
#            raise FieldError, "%s>index %s is out of range in this field. size: %s" % (self.myrank,theindex,self.datasize)

        # TODO:
        # this may fail when running in parallel, by adressing nodes outside the 
        # domain for grid. We must probably do some checking...
        return self[[closestint((i-j)/d) for (i,j,d) in zip(point,self.grid.base, self.grid.division)]]

        
    def getitem(self,index):
        """Get an element from the data. This should be
            optimized later
            
            This method accept both 1D (field) indexing and nsd-indexing as in the
            grid.
            The 1D indexing assume that indexing according to the data is used, 
            starting with the index 0, while the nsd-indexing assume that the given
            index is relative to the basenode.
        """ 

        #integer index
        if type(index) == type(1):
            return self.data[index]
        #index given as array of size 1:
        elif len(index) == 1:
            return self.data[index[0]]
        #index given as indices set:
        elif len(index) == self.grid.nsd:
            # use shapeddata:
            globindex = tuple(map(lambda x,y:x-y,index,self.basenode))
#            for i in range(len(index)):
#                globindex += (index[i]-self.basenode[i])*self.indexspec[i]
            return self.shapedata[globindex]
        else:
            raise GridError, 'unknown index given'

    def __setitem__(self,index,value):
        """operator-method for assignment to field[x], where x is of length nsd.
            This method do not check the index, it just assume that grid-indexing
            is used
        """

#        # indexing is assumed to be relative to self.basenode, that is, we
#        # subtract the number in basenode[i] in each direction i.
#        globindex = 0
#        for i in range(len(index)):
#            globindex += (index[i]-self.basenode[i])*self.indexspec[i]
#
#        # TODO: Remove debug print (and computations)
#        maxindex = len(self.data)
#        if globindex >= maxindex or globindex < 0:
#            print "*************** start ****************"
#            print self.myrank,"> Set item: ",index,self.indexspec, self.data.shape
#            print self.myrank,"> Try to set: ",globindex
#            print self.myrank,"> my base: ",self.basenode
#            print self.myrank,"> my sizespec: ",self.sizespec
#            #print self.myrank,"> mgp: ",self.grid.myGridPosition
#            #print self.myrank,"> nodebounds: ",self.grid.myGlobalNodalBoundsIncGhosts
#            print "*************** end ****************"
#
#        if globindex >= maxindex or globindex < 0:
#            raise FieldError, "Try to set index out of range in field"
        # use shapedata, adjust for basenode
        globindex = tuple(map(lambda x,y:x-y,index,self.basenode))
        self.shapedata[globindex] = value 

    def setitem(self,index,value):
        """Set an element in data with the given value. 
           This version should be optimized later
           
           This method accept both 1D (field) indexing and nsd-indexing as in 
           the grid.
           The 1D indexing assume that indexing according to the data is used, 
           starting with the index 0, while the nsd-indexing assume that the given
           index is relative to the basenode.
        """

        #integer index:
        if type(index) == type(1):
            self.data[index] = value
        #index given as array of size 1:
        elif len(index) == 1:
            self.data[index[0]] = value
        #index given as indices set:
        elif len(index) == self.grid.nsd:
#            globindex = 0
#            for i in range(len(index)):
#                globindex += (index[i]-self.basenode[i])*self.indexspec[i]
            # use shapedata, adjust for basenode
            globindex = tuple(map(lambda x,y:x-y,index,self.basenode))
            self.shapedata[globindex] = value 
        else: 
            raise GridError, 'unknown index given'

    def localIndexToDataIndex(self,index):
        # assume that index is already local (adjusted for basenode) 
        return sum([k*l for k,l in zip(index,self.shapemultiplicator)])
   
    def globIndexToDataIndex(self,index):
        # assume that index is global, not adjusted for basenode:
        globindex = tuple(map(lambda x,y:x-y, index, self.basenode))
        return sum([k*l for k,l in zip(globindex,self.shapemultiplicator)])

    def indexToPoint(self,index):
        """Accept a index-tuple as argument"""

        # TODO: seem to be duplicate of similar method in Grid - check, and remove
        #      one of them!
        base = self.grid.getBase()
        for i in xrange(len(base)):
            base[i] += index[i]*self.grid.division[i]

        return base

    def plot(self,**kwargs):
        """ Make a simple plot using pyPDE Viz-tools.
            Arguments:
              title='Title to display'
              movie='on'/'off'
            Check whether dimension is 1D or 2D, in that case, provide a simple 
            plot over the grid. Else issue a warning.

            If we run in a parallel context, we should at least make sure that 
            we do not plot ghosts.
            We have to run in non-batch mode to do plotting anyway, so it is 
            probably ok that all of us plot. 
            There should be a user-option to gather data on process 0 before 
            plotting
        """

        usemovie='off'
        usetitle='Field'
        for k in kwargs:
            if k == "movie":
                usemovie=kwargs[k]
            if k == "title":
                usetitle=kwargs[k]

        if self.isParallel:
            # do not want to plot ghosts, so create a slice list for interior nodes:
                sl = []
                for i in xrange(self.grid.nsd):
                    start = self.grid.InodeStart[i]
                    # remember: slice is up to, but not including, we need to add
                    # 1 'cause we use the spec
                    stop = start + self.grid.Ispec[i] + 1
                    tmpslice = slice(start,stop)
                    sl.append(tmpslice)

        if self.grid.nsd == 1:
            # plot with VizCurve
            #print "Make 1D plot with VizCurve"
            if self.isParallel:
                geo = self.grid.Igeometry
            else:
                geo=self.grid.geometry
            dx=self.grid.division
            x = seq(geo[0][0],geo[0][1],dx[0])
            # A better solution regarding kwargs, is to just pass the kwargs on
            #pyPDE.VizCurve.plot(x,self.data,title=usetitle,movie=usemovie)
            if self.isParallel:
                pyPDE.VizCurve.plot(x,self.shapedata[sl],**kwargs)
            else:
                pyPDE.VizCurve.plot(x,self.data,**kwargs)
        elif self.grid.nsd == 2:
            # plot with VizSurface
            # need x and y axis, and the array of values on this.
            #print "Make 2D plot with VizSurface"
            # create x and y vectors:
            if self.isParallel:
                geo = self.grid.Igeometry
            else:
                geo=self.grid.geometry
            dx=self.grid.division
            x = seq(geo[0][0],geo[0][1],dx[0])
            y = seq(geo[1][0],geo[1][1],dx[1])
            # shape data according to x and y
            xl = len(x)
            yl = len(y)
            if self.isParallel:
                z = self.shapedata[sl]
            else:
                z = _n.reshape(self.data,(xl,yl))
            # A better solution regarding kwargs, is to just pass the kwargs on
            #pyPDE.VizSurface.plot(x,y,z,color='black',line_width=0.1,color_map='def2',title=usetitle,movie=usemovie)
            pyPDE.VizSurface.plot(x,y,z,color='black',line_width=0.1,color_map='def2',**kwargs)
        elif self.grid.nsd == 3:
            if self.isParallel:
                geo = self.grid.Igeometry
            else:
                geo = self.grid.geometry
            x = seq(geo[0][0],geo[0][1],self.grid.dx)
            y = seq(geo[1][0],geo[1][1],self.grid.dy)
            z = seq(geo[2][0],geo[2][1],self.grid.dz)

            # check for a 'slice' value in kwargs:
            if 'slice' in kwargs:
                slice_value = kwargs['slice']
            else:
                # else in the middle of the second direction.
                slice_value = len(y)/2

            if self.isParallel:
                # skip ghosts.
                d = self.shapedata[sl]
            else:
                d = self.shapedata

            if 'allslice' in kwargs:
                if slice_value > 5:
                    for i in xrange(slice_value-3,slice_value+3):
                        f = d[:,i,:]
                        pyPDE.VizSurface.plot(x,z,f,color='black',line_width=0.1,color_map='def2',**kwargs)
                        time.sleep(0.4)
                else:
                    f = d[:,slice_value,:]
                    pyPDE.VizSurface.plot(x,z,f,color='black',line_width=0.1,color_map='def2',**kwargs)

            else:
                if 'face' in kwargs:
                    f = d[:,:,0]
                    z = y
                else:
                    f = d[:,slice_value,:]
                pyPDE.VizSurface.plot(x,z,f,color='black',line_width=0.1,color_map='def2',**kwargs)
                
        else:
            # Do nothing, issue a warning that plotting is not supported.
            print "Plotting in more dimensions is not implemented"

    def fill(self,function):
        """Fill the datagrid based on the function. The function need to take the
            same number of arguments as self.nsd, the argument to function will be
            a point."""
        
        if isinstance(function,(float,int)):
            function = _u.constfunc(function)

        #print "create data, fill"
        ## We may do this with tupleIterator, but when data is a 1D array, this
        ## is maybe not the most clever thing to do... ?
        ## Create a tupleIterator over the indices in the grid:
        #zerotupl = _u.getZeroTuple(self.sizespec)
        #allPoints = _u.tupleIterator(self.grid.nsd,zerotupl,self.sizespec)
        #for index in allPoints:
        #    self.setitem(index) = function(indexToPoint(index))


        # TODO: this is unefficient. Rewrite to use shape/reshape etc. from
        # numarray.
         
        # create a view of the data which resembles the grid layout
        #shapedata = _n.reshape(self.data,map(lambda x:x+1,self.sizespec))
        # fill with respect to shape

        runningind = map(lambda x: 0,range(self.grid.nsd))
        runningind[self.grid.nsd-1] -= 1
        
        for w in xrange(0,self.datasize):

            i=self.grid.nsd - 1
            while runningind[i] >= self.sizespec[i]:
                i -= 1
            # increment last possible:
            runningind[i] += 1
            # zero-out above:
            runningind[i+1:] = map(lambda x: 0, range(self.grid.nsd - i - 1))

            # TODO: Remove this debug print
            #print "w: %d, running: %s" % ( w,runningind )
            #print "Point: %s, fill value: %e" % (self.indexToPoint(runningind),function(*self.indexToPoint(runningind)))
            # unroll the arguments:
            self.data[w] = function(*self.indexToPoint(runningind))

        #print "data: ", self.data

    def fill_vec(self,function):
        """Fill the datagrid based on a function. The function 
           have to work vectorized, that is accept numarrays as input, 
           not individual points. The function should accept one argument, which
           will be a list with the coordinate vectors, shaped according to the
           grid"""

        # we will use variables from the grid
        nsd = self.grid.nsd
        geo=self.grid.geometry
        dx=self.grid.division
        # I will help with reshape each coordinate vector in the right direction
        I = _n.ones(nsd)
        # the coordinate vector, unshaped and shaped
        x = []
        xs = []
        for i in xrange(nsd):
            g = geo[i]
            x.append(seq(g[0],g[1],dx[i]))
            I[i] = len(x[i])
            xs.append(_n.reshape(x[i],I))
            I[i] = 1

        
        # create a view of the data which resembles the grid layout
        #shapedata = _n.reshape(self.data,map(lambda x:x+1,self.sizespec))
        # fill with respect to shape
        # cheavet: we create new storage her... Need to fill back in data later on.
        # hopefylly the function return rightly shaped data
        shapedata = function(*xs)
        self.data = shapedata.flat
        # since we create new data, we need to move the shapedata pointer:
        self.shapedata = self.getShapedData()

        #print "shapedata: ",shapedata
        #print "data: ", self.data

    def localnodes_slice(self):
        
        # if self.shapedata is sliced with this slice, only
        # localnodes, i.e, excluding ghosts, are used.
        
        # build a slice list:
        sl = []
        for i in xrange(self.grid.nsd):
            start = self.grid.InodeStart[i]
            # remember: slice is up to, but not including, we need to add
            # 1 'cause we use the spec
            stop = start + self.grid.Ispec[i] + 1
            tmpslice = slice(start,stop)
            sl.append(tmpslice)
        return sl



# A small helperclass. I need a map-multiplier with memory. I figured out this
# was the easiest way to do it.  
class MapHelper:
    def __init__(self,y=0):
        self.y = y

    def add(self,x):
        self.y += x
        return self.y

    def multiply(self,x):
        self.y *= x
        return self.y

    def reset(self,y=0):
        self.y = y

