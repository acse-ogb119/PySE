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
#include "MapMatSparse.h"
#include <iostream>


int main() {

// Make a sparse matrix

//    vector<int> col(3);
//    vector<int> row(1+1);
//    vector<int> m2i(1);
//    vector<double> values(3);
//
//    values[0] = 1.0; values[1] = 2.0; values[2] = 3.0;
//    col[0] = 0; col[1] = 1; col[2] = 2;
//
//    row[0] = 0; row[1] = 3; m2i[0] = 1;
//
//    DegMatSparse A(col, row, m2i, values);
//
// Make a map-matrix:

    MapMatSparse B;
    B(0,0) = 0.0;
    B(1,1) = 1.0;
    B(2,2) = 2.0;
    B(0,3) = 19;
    B(3,3) = 3.0;
    B(4,4) = 4.0;
    B(5,5) = 5.0;
    
    string filename, name;
    filename = "f.m";
    name = "X";


    B.save(filename, name);


// Convert!
    DegMatSparse C(B);

// Construct a right hand side and a result vector
    
    valarray<double>* b = new valarray<double>(6); *b = 0.0;
    valarray<double>* b2= new valarray<double>(6); *b2 = 0.0;
    valarray<double>* x = new valarray<double>(6); *x = 1.0;



    B.prod(x, b); // Compute b = A*x
    C.prod(x, b2); // Compute b = A*x

    for (int i=0;i<6;i++) {
        cout << "b["  << i << "] - b2[" << i << "] =" << (*b)[i] -(*b2)[i]  << endl;
//        cout << "b2[" << i << "]=" << (*b2)[i] << endl;
    }

    return 0;
}
