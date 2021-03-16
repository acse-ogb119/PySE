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

#include "DMS.h"
#include <iostream>

//int interpolationSearch(int *vec, int key, int low, int high) {
//  // Interpolation search for the index of "key"
//    double low_diff, range_diff, count_diff;
//    while ( vec[high] >= key && key > vec[low] ) {
//        low_diff = (double)key - vec[low];
//        range_diff = (double)vec[high] - vec[low];
//        count_diff = (double)high - low;
//        int range = (int)( low_diff / range_diff * count_diff + low );
//        if ( key > vec[range] )
//            low = range + 1;
//        else if ( key < vec[range] )
//            high = range - 1;
//        else
//            low = range;
//    } 
//    if ( key == vec[low] ) 
//        return low;
//    else 
//        return -1;
//};
//


DMS:: DMS(){
    n=nnz=nnzr=offset=0;
}

DMS:: DMS(int _n, int _nnz, int _nnzr) {
    n = _n; nnz = _nnz; nnzr = _nnzr;
    col = (int *)malloc(nnz*sizeof(int));
    row = (int *)malloc((nnzr + 1)*sizeof(int));
    m2i = (int *)malloc(nnzr*sizeof(int));
    values = (double *)malloc(nnz*sizeof(double));
}


DMS:: DMS(int _n, int _nnz, int _nnzr, int *_col, int *_row, int *_m2i, double *_values){
    n    = _n;      nnz = _nnz;
    nnzr = _nnzr;   col = _col;
    row  = _row;    m2i = _m2i;
    values = _values;
    offset = 0;
}

DMS:: DMS(MapMatSparse& mmat){
/** Construct a DMS from a MapMatSparse **/
    offset = 0;

    vector<int> m2itmp;
    vector<int> rowtmp;

    map<pair<int, int>, double>::iterator iter = mmat.values.begin();
    map<pair<int, int>, double>::iterator end = mmat.values.end();

    nnz = mmat.values.size();
    col    = (int *)malloc(nnz*sizeof(int));
    values = (double *)malloc(nnz*sizeof(double));

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
        values[k] = (*iter).second;
        col[k] = j;
        iter++;
        k++;
    }
    rowtmp.push_back(k);

    nnzr = m2itmp.size();
    n = m2itmp[nnzr-1];
    row = (int *)malloc((nnzr+1)*sizeof(int));
    m2i = (int *)malloc((nnzr)*sizeof(int));
    for (int i=0;i<nnzr;i++) {
        row[i] = rowtmp[i];
        m2i[i] = m2itmp[i];
    }
    row[nnzr] = rowtmp[nnzr];
}

double DMS:: operator()(int i, int j) {
    int m = interpolationSearch(m2i, i, 0, nnzr-1);    
    if (m == -1) return 0.0;
    int r = interpolationSearch(col, j, row[m], row[m+1]-1);
    if (r == -1) return 0.0;
    return values[r];
}

//void DMS:: prod(double *x, double *b) {
//    int i;
//    int nnzr = (*m2i).size();
////    cout << "nnzr=" << nnzr << endl;
//    for (int m = 0; m<nnzr; m++) {
//        i = (*m2i)[m];
//        (*b)[i] = 0.0;
//        for ( int k = (*row)[m]-offset; k < (*row)[m+1]-offset; k++) {
//            (*b)[i] += (*values)[k]*(*x)[(*col)[k]-offset];
//        }
////        cout << b[i] << endl;
//    }
//}

void DMS:: prod(double *x, double *b) {
    int i;
    for (int m = 0; m<nnzr; m++) {
        i = m2i[m];
        b[i] = 0.0;
        for ( int k = row[m]-offset; k < row[m+1]-offset; k++) {
            b[i] += values[k]*x[col[k]-offset];
        }
    }
}
