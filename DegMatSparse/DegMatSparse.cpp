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

#include "DegMatSparse.h"
#include <iostream>



DegMatSparse:: DegMatSparse(){
    col = new valarray<int>;
    row = new valarray<int>;
    m2i = new valarray<int>;
    values = new valarray<double>;
    offset = 0;
}

DegMatSparse:: DegMatSparse(valarray<int>&  _col, valarray<int>&  _row, valarray<int>&  _m2i, valarray<double>&  _values){
    col = &_col;
    row = &_row;
    m2i = &_m2i;
    values = &_values;
    offset = 0;
}

DegMatSparse:: DegMatSparse(MapMatSparse& mmat){
/** Construct a DegMatSparse from a MapMatSparse **/
    offset = 0;

    vector<int> m2itmp;
    vector<int> rowtmp;

    map<pair<int, int>, double>::iterator iter = mmat.values.begin();
    map<pair<int, int>, double>::iterator end = mmat.values.end();

    col     = new valarray<int>(mmat.values.size());
    values  = new valarray<double>(mmat.values.size());

    int i, j;
    int i_prev = -1;
    int k = 0;
    while (iter != end) {
        i = (*iter).first.first;
        j = (*iter).first.second;
        if (i != i_prev){ 
            m2itmp.push_back(i);
            rowtmp.push_back(k);
            i_prev = i;
        }
        (*values)[k] = (*iter).second;
        (*col)[k] = j;
        iter++;
        k++;
    }
    rowtmp.push_back(k);
    int nrows = m2itmp.size();
    row = new valarray<int>(nrows+1);
    m2i = new valarray<int>(nrows);
    for (int i=0;i<nrows;i++) {
        (*row)[i] = rowtmp[i];
        (*m2i)[i] = m2itmp[i];
    }
    (*row)[nrows] = rowtmp[nrows];
}

double DegMatSparse:: operator()(int i, int j) {
    int m = interpolationSearch(*m2i, i, 0, (*m2i).size()-1);    
    if (m == -1) return 0.0;
    int r = interpolationSearch(*col, j, (*row)[m], (*row)[m+1]-1);
    if (r == -1) return 0.0;
    return (*values)[r];
}

void DegMatSparse:: prod(valarray<double>* x, valarray<double>* b) {
    int i;
    int nnzr = (*m2i).size();
//    cout << "nnzr=" << nnzr << endl;
    for (int m = 0; m<nnzr; m++) {
        i = (*m2i)[m];
        (*b)[i] = 0.0;
        for ( int k = (*row)[m]-offset; k < (*row)[m+1]-offset; k++) {
            (*b)[i] += (*values)[k]*(*x)[(*col)[k]-offset];
        }
//        cout << b[i] << endl;
    }
}

void DegMatSparse:: prod2(double *x, double *b) {
    int i;
    int nnzr = (*m2i).size();
    for (int m = 0; m<nnzr; m++) {
        i = (*m2i)[m];
        b[i] = 0.0;
        for ( int k = (*row)[m]-offset; k < (*row)[m+1]-offset; k++) {
            b[i] += (*values)[k]*x[(*col)[k]-offset];
        }
    }
}
