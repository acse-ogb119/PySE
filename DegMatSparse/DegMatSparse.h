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

#ifndef DEGMATSPARSE_IS_INCLUDED
#define DEGMATSPARSE_IS_INCLUDED


/**
 * The main idea is that b = Ax, in such a way that
 * when b is initialized, say b_i = c_i
 * b[i] = a[i,:]*x where 
 *  i is a non-zero row
 *  x is a vector
 *  a is a Degenerate compressed row storage matrix
 *
 *  If i = j is a zero-row in A, b_i = c_i.
 * */

#include <vector>
#include <valarray>
#include "MapMatSparse.h"
#include "isearch.h"

using namespace std;

class DegMatSparse {
    public:
        valarray<int> *col;
        valarray<int> *row;
        valarray<int> *m2i;
        valarray<double> *values;
        int offset;

        DegMatSparse();
        DegMatSparse(valarray<int>&  _col, valarray<int>&  _row, valarray<int>&  _m2i, valarray<double>&  _values);
        DegMatSparse(MapMatSparse& mmat);
        ~DegMatSparse(){};

        void prod(valarray<double>* x, valarray<double>* b); // Compute b = A*x
        void prod2(double *x, double *b); // Compute b = A*x
        double operator()(int i, int j) ;
        int numberOfRows() {return (*m2i)[(*m2i).size()-1]+1;};

};
#endif
