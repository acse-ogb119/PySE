# -*- coding: utf-8 -*-
#
# (c) Copyright 2005
#     Author: Åsmund Ødegård
#             Ola Skavhaug
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
# ConjGrad: A generic implementation of the conjugate gradient method. 
#


ConjGradError = "Error in ConjGrad"

#import Numeric/numarray
from Utils import _n,inner
from math import sqrt

#def inner(u,v):
#    """Compute innerproduct of u and v.
#       It is not computed here, we just check type of operands and
#       call a suitable innerproduct method"""
#
#    if isinstance(u,_n.ArrayType):
#        # assume both are numarrays:
#        return _n.dot(u,v)
#    else:
#        # assume that left operand implements inner
#        return u.inner(v)

def conjgrad(A, x, b, tolerance=1.0E-05, relativeconv=False):
    """
    conjgrad(A,x,b): Solve Ax = b with the conjugate gradient method.

    Arguments: A: a matrix or discrete operator. Must support A*x to operate on x.

               x: anything that A can operate on, as well as support inner product 
                  and multiplication/addition/subtraction with scalar and something
                  of the same type

               b: right-hand side, probably the same type as x.

               tolerance: the convergence criterion

               relativeconv: boolean, control relative vs. absolute convergence 
                  criterion. That is, if relativeconv is True, ||r||_2/||r_init||_2 
                  is used as convergence criterion, else just ||r||_2 is used.
                  Actually ||r||_2^2 is used, since that save a few cycles.

    Returns:   x: the solution. 
    """

    #print "Solving with conjgrad"

    x = x.copy()
    r = b - A*x
    p = r.copy()
    r0 = inner(r,r)   
    if relativeconv:
        #rS = sqrt(r0)  # store rS in case relative convergence crit. is selected
        #actually, I think it is the norm of b i should relate to!
        #rS = sqrt(inner(b,b))
        tolerance *= sqrt(inner(b,b))
        #rS = 1.0
    #else:
        #rS = 1.0

    # later on, include iter in return (as tuple?) if user want that for monitoring 
    # convergence.
    iter = 0
    #print "Starting, ||r||=%e, tolerance=%e" % (sqrt(r0),tolerance)
    while sqrt(r0) > tolerance:
        w = A*p
        a = r0/inner(w,p)
        x += a*p
        #x.plot()
        r -= a*w
        r1 = inner(r,r)
        p = r + (r1/r0)*p
        r0 = r1
        iter += 1
        print "ConjGrad keep going..., iter: %d, ||e||=%e" % (iter,sqrt(r0))
    return x
