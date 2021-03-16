# -*- coding: utf-8 -*-
#
# (c) Copyright 2003, 2004, 2005
#     Author: Hans Petter Langtanen,
#             Ola Skavhaug
#             Åsmund Ødegård
#     Simula Research Laboratory AS
#     
#     This file is part of PyPDE.
#     As PyPDE is unfinshed and not distributed, this is distributed as part
#     of PyFDM / PySE.
#
#     PyPDE is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     PyPDE is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with PyPDE; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



"""
Transparent use of Numeric and numarray.
The module also has some extensions to Numerical Python.

Important features of the unified interface:

RandomArray and LinearAlgebra are imported and available under these
names

Some frequently standard modules like sys, os, and operator are
imported (and available if you do from py4cs.numpytools import *)

sequence   : sequence(a,b,s) computes numbers from a up to and
             including b in steps of s and (default) type Float
seq        : same as sequence (short form)

isequence  : as sequence, but integer counters are computed
             (isequence is an alternative to range where the
             upper limit is included in the sequence)
iseq       : same as isequence (short form)

NumPyArray : the type used in isinstance(a,NumPyArray) for
             checking if a is a NumPy array

amin/amax  : compute maximum and minimum of array entries

exp_robust : special version of the exponential function that
             handles very small and very large arguments without
             raising exceptions

array_output_precision(n) : print arrays with n decimals

solve_tridiag_linear_system : returns the solution of a tridiagonal
             linear system

wrap2callable : tool for turning constants, discrete data, string
             formulas, function objects, or plain functions
             into an object that behaves as a function
"""

import os, sys, operator

py4csError = "Problem in pyPDE.numpytools"

use_numarray = False
use_numeric = False
try_import = True
# First, we check whether the user have a preference regarding numarray/Numeric:
# if the choosen package is already imported, there is no need to try the import
if os.environ.get('NUMPYARRAY','') == 'numarray':
    use_numarray = True
    if sys.modules.has_key('numarray'):
        try_import = False
elif os.environ.get('NUMPYARRAY','') == 'Numeric':
    use_numeric = True
    if sys.modules.has_key('Numeric'):
        try_import = False
elif sys.modules.has_key('numarray'):
    # user have choosen by importing something:
    use_numarray = True
    try_import = False
elif sys.modules.has_key('Numeric'):
    use_numeric = True
    try_import = False

if try_import:
    # the above test either told us that the user have a preference, but have
    # not imported anything relevant, or that user have no preference at all
    # (either by environment or by import)
    if use_numeric:
        try:
            from Numeric import *
        except:
            raise py4csError, "Unable to import Numeric as requested"
    elif use_numarray:
        try:
            from numarray import *
        except:
            raise py4csError, "Unable to import numarray as requested"
    else:
        # First try numarray, then Numeric, then raise exception:
        try: 
            from numarray import *
            use_numarray = True
        except:
            try: 
                from Numeric import *
                use_numeric = True
            except:
                raise py4csError, "Unable to import both numarray and Numeric"
                

# import Numeric or numarray?
#if os.environ.get('NUMPYARRAY','') == 'numarray' or \
#       sys.modules.has_key('numarray'):
if use_numarray:
    basic_NumPy = 'numarray'
    from numarray import *
    import numarray.random_array.RandomArray2 as RandomArray
    import numarray.linear_algebra.LinearAlgebra2 as LinearAlgebra

    def array_output_precision(no_of_decimals):
        """Set no of decimals in printout of arrays."""
        arrayprint.set_precision(no_of_decimals)

    def amax(a):
        """Compute the maximum of the entries in a."""
        try:
            return a.max()
        except AttributeError:
            # not a NumPy array
            if operator.isSequenceType(a):
                return max(a)
            elif operator.isNumberType(a):
                return a
            else:
                raise TypeError, 'amax of %s not supported' % type(a)        

    def amin(a):
        """Compute the minimum of the entries in a."""
        try:
            return a.min()
        except AttributeError:
            # not a NumPy array
            if operator.isSequenceType(a):
                return min(a)
            elif operator.isNumberType(a):
                return a
            else:
                raise TypeError, 'amin of %s not supported' % type(a)
        
    NumPyArray = NumArray
#else:
#elif os.environ.get('NUMPYARRAY','') == 'Numeric' or \
#         sys.modules.has_key('Numeric'):
elif use_numeric:
    basic_NumPy = 'Numeric'
    # doesn't hurt, and the user may have done import Numeric.
    from Numeric import *
    import RandomArray
    import LinearAlgebra
    # dangerous; they have different data attributes...
    from UserArray import UserArray as NumArray

    def array_output_precision(no_of_decimals):
        """Set no of decimals in printout of arrays."""
        sys.float_output_precision = no_of_decimals

    def amax(a):
        """Compute the maximum of the entries in a."""
        try:
            return max(a.flat)  # use Python's list min
        except AttributeError:
            # not a NumPy array
            if operator.isSequenceType(a):
                return max(a)
            elif operator.isNumberType(a):
                return a
            else:
                raise TypeError, 'amax of %s not supported' % type(a)

    def amin(a):
        """Compute the minimum of the entries in a."""
        try:
            return min(a.flat)
        except AttributeError:
            # not a NumPy array
            if operator.isSequenceType(a):
                return min(a)
            elif operator.isNumberType(a):
                return a
            else:
                raise TypeError, 'amin of %s not supported' % type(a)
            
    NumPyArray = ArrayType
else:
    raise py4csError, "No Numeric or numarray library available"

# with straight import one can say
# import Numeric as N or import numarray as N

def exp_robust(x):
    """
    Exponential function without over/under-flow.
    Works with arrays.
    """
    if x < -200:
        r = 0.0
    elif x > 700:
        #raise OverflowError, 'exp(%g) caused overflow' % x
        r = 1.0E+999  # inf
    else:
        r = exp(x)
    return r


def sequence(min=0.0, max=None, inc=1.0, type=Float,
             return_type='NumPyArray'):
    """
    Generate numbers from min to (and including) max,
    with increment of inc. Alternative to arange/arrayrange.
    """
    if max is None: # allow sequence(3) to be 0., 1., 2., 3.
        # take 1st arg as max, min as 0, and inc=1
        max = min; min = 0.0; inc = 1.0
    r = arrayrange(min, max + inc/2.0, inc, type)
    if return_type == 'NumPyArray':
        return r
    elif return_type == 'list':
        return r.tolist()
    elif return_type == 'tuple':
        return tuple(r.tolist())

seq = sequence # short form

def isequence(start=0, stop=None, inc=1):
    """
    Generate integers from start to (and including) stop,
    with increment of inc. Alternative to range/xrange.
    """
    if stop is None: # allow isequence(3) to be 0, 1, 2, 3
        # take 1st arg as stop, start as 0, and inc=1
        stop = start; start = 0; inc = 1
    return xrange(start, stop+inc, inc)

iseq = isequence


def arr(shape=None, element_type=Float, data=None, copy=True):
    """
    Compact and flexible interface for creating NumPy arrays.

    Input arguments:
    shape               length of each dimension
    data                list, tuple, or NumPy array
    copy                True: copy data, False: share memory with data
    element_type        Float, Int, Complex, Float32, etc.

    The array can be created as zeros (just shape specified), or
    a copy of or reference to (depending on copy=True,False resp.)
    a list, tuple, or NumPy array (provided as the data argument).

    The function calls the underlying NumPy functions zeros and array
    (see the NumPy manual for the functionality of these functions).

    >>> arr((3,4))
    array([[ 0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.]])
    >>> somelist=[[0,1],[5,5]]
    >>> a = arr(data=somelist)
    >>> b = a + 1
    >>> c = arr(data=b, copy=False)  # let c share data with b
    >>> b is c
    True
    >>> id(b) == id(c)
    True
    >>> arr(4, element_type=Int) + 4  # integer array
    array([4, 4, 4, 4])
    """
    if data is not None:

        if not operator.isSequenceType(data):
            raise TypeError, 'arr: data argument is not a sequence type'
        
        if isinstance(shape, (list,tuple)):
            # check that shape and data are compatible:
            if reduce(operator.mul, shape) != size(data):
                raise ValueError, \
                      'arr: shape=%s is not compatible with %d '\
                      'elements in the provided data' % (shape, size(data))
        elif isinstance(shape, int):
            if shape != size(data):
                raise ValueError, \
                      'arr: shape=%d is not compatible with %d '\
                      'elements in the provided data' % (shape, size(data))
        elif shape is None:
            if isinstance(data, (list,tuple)) and copy == False:
                print 'arr: cannot share data (data is list/tuple)'
            return array(data, element_type, copy=copy)
        else:
            raise TypeError, \
                  'shape is %s, must be list/tuple or int' % type(shape)
    else:
        # no data, no from_function, make just zeros

        if isinstance(shape, NumPyArray):
            raise TypeError, \
           'arr: shape (1st arg) is NumPy array, must be tuple or int'
        if shape is None:
            raise ValueError, \
            'arr: either shape, data, or from_function must be specified'

        return zeros(shape, element_type)

            
def solve_tridiag_linear_system(A, b):
    """
    Solve a tridiagonal linear system of the form
    A[0,1]*x[0] + A[0,2]*x[1]                                        = 0
    A[1,0]*x[0] + A[1,1]*x[1] + A[1,2]*x[2]                          = 0
    ...
    ...
             A[k,0]*x[k-1] + A[k,1]*x[k] + A[k,2]*x[k+1]             = 0
    ...
                 A[n-2,0]*x[n-3] + A[n-2,1]*x[n-2] + A[n-2,2]*x[n-1] = 0
    ...
                                   A[n-1,0]*x[n-2] + A[n-1,1]*x[n-1] = 0

    That is, the diagonal is stored in A[:,1], the subdiagonal
    is stored in A[1:,0], and the superdiagonal is stored in A[:n-2,2].

    The storage is not memory friendly in Python/C (diagonals in
    the columns of A, but when A is sent to F77 for high-performance
    computing, a copy is taken and the F77 routine works with the
    same algorithm and hence optimal Fortran storage.
    """
    n = len(b)
    x = zeros(n, Float)  # solution
    d = zeros(n, Float);  c = zeros(n, Float);  m = zeros(n, Float)

    d[0] = A[0,1]
    c[0] = b[0]

    for k in isequence(1, n-1, 1):
        m[k] = A[k,0]/d[k-1]
        d[k] = A[k,1] - m[k]*A[k-1,2]
        c[k] = b[k] - m[k]*c[k-1]
    x[n-1] = c[n-1]/d[n-1]

    # back substitution:
    for k in isequence(n-2, 0, -1):
        x[k] = (c[k] - A[k,2]*x[k+1])/d[k]
    return x




try:
    import Pmw
    class NumPy2BltVector(Pmw.Blt.Vector):
        """
        Copy a NumPy array to a BLT vector:
        # a: some NumPy array
        b = NumPy2BltVector(a)  # b is BLT vector
        g = Pmw.Blt.Graph(someframe)
        # send b to g for plotting
        """
        def __init__(self, array):
            Pmw.Blt.Vector.__init__(self, len(array))
            self.set(tuple(array))  # copy elements
except:
    class NumPy2BltVector:
        def __init__(self, array):
            raise ImportError, "Python is not compiled with BLT"

try:
    from py4cs.StringFunction import StringFunction
except:
    pass  # wrap2callable may not work

class WrapNo2Callable:
    """Turn a number (constant) into a callable function."""
    def __init__(self, constant):
        self.constant = constant

    def __call__(self, *args):
        """
        >>> w = WrapNo2Callable(4.4)
        >>> w(99)
        4.4000000000000004
        >>> # try vectorized computations:
        >>> x = seq(1, 4, 1)
        >>> y = seq(1, 2)
        >>> xv = x[:,NewAxis]; yv = y[NewAxis,:]
        >>> xv + yv
        array([[ 2.,  3.],
               [ 3.,  4.],
               [ 4.,  5.],
               [ 5.,  6.]])
        >>> w(xv, yv)
        array([[ 4.4,  4.4],
               [ 4.4,  4.4],
               [ 4.4,  4.4],
               [ 4.4,  4.4]])
        """
        if isinstance(args[0], NumPyArray):
            # vectorized version:
            r = args[0].copy()
            # to get right dimension of the return array,
            # compute with args in a simple formula (sum of args)
            for a in args[1:]:
                r = r + a  # in-place r+= won't work
                # (handles x,y,t - the last t just adds a constant)
            r[:] = self.constant
            return r
        else:
            # scalar version:
            return self.constant

class WrapDiscreteData2Callable:
    """
    Turn discrete data on a uniform grid into a callable function,
    i.e., equip the data with an interpolation function.

    >>> from py4cs.numpytools import *
    >>> x = seq(0,1,0.1)
    >>> y = 1+2*x
    >>> f = wrap2callable((x,y))
    >>> f(0.5)   # evaluate f(x)
    1.5
    >>> f(0.5, 0.1)  # discrete data with extra time prm: f(x,t)
    1.5
    """
    def __init__(self, data):
        self.data = data  # (x,y,f) data for an f(x,y) function
        from Scientific.Functions.Interpolation \
             import InterpolatingFunction # from ScientificPython
        self.interpolating_function = \
             InterpolatingFunction(self.data[:-1], self.data[-1])
        self.ndims = len(self.data[:-1])  # no of spatial dim.
        
    def __call__(self, *args):
        # allow more arguments (typically time) after spatial pos.:
        args = args[:self.ndims]
        # args can be tuple of scalars (point) or tuple of vectors
        if isinstance(args[0], (float, int)):
            return self.interpolating_function(*args)
        else:
            # args is tuple of vectors; Interpolation must work
            # with one point at a time:
            r = [self.interpolating_function(*a) for a in zip(*args)]
            return array(r)  # wrap in NumPy array

        
def wrap2callable(f, **kwargs):
    """
    Allow constants, string formulas, discrete data points,
    user-defined functions and (callable) classes to be wrapped
    in a new callable function. That is, all the mentioned data
    structures can be used as a function, usually of space and/or
    time.
    (kwargs is used for string formulas)

    >>> from py4cs.numpytools import *
    >>> f1 = wrap2callable(2.0)
    >>> f1(0.5)
    2.0
    >>> f2 = wrap2callable('1+2*x')
    >>> f2(0.5)
    2.0
    >>> f3 = wrap2callable('1+2*t', independent_variables='t')
    >>> f3(0.5)
    2.0
    >>> f4 = wrap2callable('a+b*t')
    >>> f4(0.5)
    Traceback (most recent call last):
    ...
    NameError: name 'a' is not defined
    >>> f4 = wrap2callable('a+b*t', independent_variables='t', \
                           a=1, b=2)
    >>> f4(0.5)
    2.0

    >>> x = seq(0,1,0.5); y=1+2*x
    >>> f5 = wrap2callable((x,y))
    >>> f5(0.5)
    2.0
    >>> def myfunc(x):  return 1+2*x
    >>> f6 = wrap2callable(myfunc)
    >>> f6(0.5)
    2.0
    >>> f7 = wrap2callable(lambda x: 1+2*x)
    >>> f7(0.5)
    2.0
    >>> class MyClass:
            def __call__(self, x):
                return 1+2*x

    >>> myclass = MyClass()
    >>> f8 = wrap2callable(myclass)
    >>> f8(0.5)
    2.0
    >>> # 3D functions:
    >>> f9 = wrap2callable('1+2*x+3*y+4*z', \
                           independent_variables=('x','y','z'))
    >>> f9(0.5,1/3.,0.25)
    4.0
    >>> # discrete 3D data:
    >>> y = seq(0,1,0.5); z = seq(-1,0.5,0.1)
    >>> xc = x[:,NewAxis,NewAxis]; yc = y[NewAxis,:,NewAxis]
    >>> zc = z[NewAxis,NewAxis,:]
    >>> def myfunc3(x,y,z):  return 1+2*x+3*y+4*z

    >>> values = myfunc3(xc,yc,zc)
    >>> f10 = wrap2callable((x,y,z,values))
    >>> f10(0.5,1/3.,0.25)
    4.0
    >>> 
    """
    if isinstance(f, str):
        return StringFunction(f, **kwargs)
    elif isinstance(f, (float, int)):
        return WrapNo2Callable(f)
    elif isinstance(f, (list,tuple)):
        return WrapDiscreteData2Callable(f)
    elif operator.isCallable(f):
        return f
    else:
        raise TypeError, 'f of type %s is not callable' % type(f)


# problem: setitem in ArrayGen does not support multiple indices
# relying on inherited __setitem__ works fine

def _doctest():
    import doctest, numpytools
    return doctest.testmod(numpytools)

if __name__ == '__main__':
    test_ArrayGen()
    #_doctest()  # does not work properly with wrap2callable
    
