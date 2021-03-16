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

#include "FastMatSparse.h"
#include <iostream>


FastMatSparse:: ~FastMatSparse(){
    free(col);
    free(row);
    free(values);
}

FastMatSparse:: FastMatSparse(int _n, int _nnz){
    n = _n;
    nnz = _nnz;
    col = (int *)malloc(n*sizeof(int));
    row = (int *)malloc((nnz+1)*sizeof(int));
    values = (double *)malloc(nnz*sizeof(double));
    offset = 0;
}

FastMatSparse:: FastMatSparse(int _n, int _nnz, int*  _col, int*  _row, double*  _values){
    n = _n;
    nnz = _nnz;
    col = _col;
    row = _row;
    values = _values;
    offset = 0;
}

FastMatSparse:: FastMatSparse(MapMatSparse& mmat){
/** Construct a FastMatSparse from a MapMatSparse **/
    offset = 0;

    vector<int> rowtmp;

    map<pair<int, int>, double>::iterator iter = mmat.values.begin();
    map<pair<int, int>, double>::iterator end = mmat.values.end();

    nnz = mmat.values.size();
    col     = (int *)malloc(nnz*sizeof(int));
    values  = (double *)malloc(nnz*sizeof(double));

    int i, j;
    int i_prev = -1;
    int k = 0;
    while (iter != end) {
        i = (*iter).first.first;
        j = (*iter).first.second;
        if (i != i_prev){ 
            rowtmp.push_back(k);
            i_prev = i;
        }
        values[k] = (*iter).second;
        col[k] = j;
        iter++;
        k++;
    }
    n = rowtmp.size();
    rowtmp.push_back(k);
    row =  (int *)malloc((n+1)*sizeof(int));
    for (int i=0;i<n;i++) {
        row[i] = rowtmp[i];
    }
    row[n] = rowtmp[n];
}

double FastMatSparse:: operator()(int i, int j) {
    int r = interpolationSearch(col, j, row[i], row[i+1]-1);
    if (r == -1) return 0.0;
    return values[r];
}

void FastMatSparse:: prod(double *x, double *b) {
    int i;
    for (int i = 0; i<n; i++) {
        b[i] = 0.0;
        for ( int k = row[i]-offset; k < row[i+1]-offset; k++) {
            b[i] += values[k]*x[col[k]-offset];
        }
    }
}
