# -*- coding: iso-8859-1 -*-
#
# (c) Copyright 2005
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
# Jacobi: A generic implementation of the jacobi iteration
#

from Utils import inner
from math import sqrt

JacobiError = "Error in Jabobi"

def jacobi(A, x, b, tolerance=1.0E-05, relativeconv=False):
    """
    jacobi(A, x, b): Solve Ax = b with the jacobi iteration.

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
    # compute current residual:
    # x and b may be the same vector, so do not change any of them!
    #
    #print "Solve with Jacobi, tolerance: ",tolerance
    #if relativeconv: 
    #    print "Using relative convergence"
    r = b - A*x
    xn = x + r
    r0 = inner(r,r)
    if relativeconv:
        #rS = sqrt(r0)
        # scale against the norm of the right-hand side.
        tolerance *= sqrt(inner(b,b))
    #else:
    #    rS = 1.0
    #
    iter = 0
    #print "Starting, ||r||=%e, tolerance=%e" % (sqrt(r0),tolerance)
    while sqrt(r0) > tolerance:
        #print "Residual norm is %e" % (r.norm())
        xp = xn
        r = b - A*xp
        xn = xp + r
        #xn.plot()
        r0 = inner(r,r)
        iter += 1
        print " iter: %d, ||r||=%e" % (iter,sqrt(r0))
    return xn

