/**
     (c) Copyright 2005
     Author: Ola Skavhaug
     Simula Research Laboratory AS
     
     This file is part of DegMatSparse.

     DegMatSparse is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

     DegMatSparse is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DegMatSparse; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

**/

#ifndef FASTMATSPARSE_IS_INCLUDED
#define FASTMATSPARSE_IS_INCLUDED


/**
 * The main idea is that b = Ax, in such a way that
 * when b is initialized, say b_i = c_i
 * b[i] = a[i,:]*x where 
 *  i is a non-zero row
 *  x is a vector
 *  a is a Fastenerate compressed row storage matrix
 *
 *  If i = j is a zero-row in A, b_i = c_i.
 * */

#include <vector>
#include <valarray>
#include "MapMatSparse.h"
#include "isearch.h"

using namespace std;

class FastMatSparse {
    public:
        int n;
        int nnz;
        int *col;
        int *row;
        double *values;
        int offset;

        FastMatSparse() {};
        ~FastMatSparse();
        FastMatSparse(int n, int nnz);
        FastMatSparse(int n, int nnz, int *_col, int *_row, double *_values);
        FastMatSparse(MapMatSparse& mmat);

        void prod(double *xptr, double *bptr); // Compute b = A*x
        double operator()(int i, int j);
        int numberOfRows() {return n;};
        int numberOfNonZeros() {return nnz;};

};
#endif
