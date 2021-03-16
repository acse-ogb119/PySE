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

#include <fstream>
#include <iostream>
#include "MapMatSparse.h"
#include "DegMatSparse.h"
#include "FastMatSparse.h"
#include <math.h>
#include <sys/times.h>

using namespace std;

int main() {
    MapMatSparse m;
    int nloops = 10;

//    ifstream input ("../../data/A.m");
//    char buffer[256];
//    string content;
//    input.getline(buffer, 255);
//    while (!input.eof()){
//        content = buffer;
//        if (sscanf(content.c_str(), "A(%d,%d) =%f;", i, j, f) >= 3) {
////            if (fabs(f[0]) > 1e-10) {
//                m((*i)-1,(*j)-1) = f[0];
////            }
//        }
//        input.getline(buffer, 255);
//    }
//
    cout << "Loading data..." << endl ;
    m.load("../../data/A.m");

    cout << "Done" << endl;

//    MapMatSparse B;
//    B(0,0) = 0.0;
//    B(1,1) = 1.0;
//    B(2,2) = 2.0;
//    B(0,3) = 19;
//    B(3,3) = 3.0;
//    B(4,4) = 4.0;
//    B(5,5) = 5.0;
// 
    DegMatSparse d(m);
    FastMatSparse e(m);
    int n = d.numberOfRows();
    cout << "Rows = " << n << endl;
    valarray<double> *x = new valarray<double>(n);
    (*x) = 1.0;
    valarray<double> *b = new valarray<double>(n);  (*b)  = 0.0;
    valarray<double> *b2 = new valarray<double>(n); (*b2) = 0.0;

    cout << "Valarrays created" << endl;

    double *bptr = (double *)malloc(n*sizeof(double));
    double *xptr = (double *)malloc(n*sizeof(double));
    double *fbptr = (double *)malloc(n*sizeof(double));
    double *fxptr = (double *)malloc(n*sizeof(double));

    for (int i = 0; i<n; i++) {
        xptr[i] = 1.0;
        bptr[i] = 0.0;
        fxptr[i] = 1.0;
        fbptr[i] = 0.0;
    }


    cout << "C-arrays created" << endl;

//    ifstream input2 ("../../data/vec.m");
//    input2.getline(buffer, 255);
//    while (!input2.eof()){
//        content = buffer;
//        if (sscanf(content.c_str(), "X(%d) =%f;", i, f) >= 2) {
//            x[(*i)-1] = f[0];
//        }
//        input2.getline(buffer, 255);
//    }
//
    /* Timing structures */

    tms *t0 = new tms(); tms *t1 = new tms();
    int cps =  sysconf(_SC_CLK_TCK);
    double secmap, secdeg, secptr, secfast;

    cout << "Running tests, please wait"<<endl;

    /* Run the tests */

    times(t0); 
    for (int i=0; i<nloops;i++) d.prod2(xptr,bptr);
    times(t1);
    secptr = ((*t1).tms_utime - (*t0).tms_utime)*1.0/cps;


    times(t0); 
    for (int i=0; i<nloops;i++) m.prod(x,b);
    times(t1);
    secmap = ((*t1).tms_utime - (*t0).tms_utime)*1.0/cps;

    times(t0); 
    for (int i=0; i<nloops;i++) d.prod(x,b2);
    times(t1);
    secdeg = ((*t1).tms_utime - (*t0).tms_utime)*1.0/cps;

    times(t0); 
    for (int i=0; i<nloops;i++) e.prod(fxptr,fbptr);
    times(t1);
    secfast = ((*t1).tms_utime - (*t0).tms_utime)*1.0/cps;


    /* Write out the results */

    cout << "Time for "<< nloops << " MapMatSparse products: " << secmap << endl;
    cout << "Time for "<< nloops << " DegMatSparse products: " << secdeg << endl;
    cout << "Time for "<< nloops << " DegMatSparse ptr-products: " << secptr << endl;
    cout << "Time for "<< nloops << " FastMatSparse ptr-products: " << secfast << endl;
    cout << " secmap/secdeg = " << secmap/secdeg << endl;

    double normm = 0.0;
    double normd = 0.0;
    double normp = 0.0;
    double normf = 0.0;
    for (int k=0; k< n; k++){
        normm += (*b)[k]*(*b)[k];
        normd += (*b2)[k]*(*b2)[k];
        normp += bptr[k]*bptr[k];
        normf += fbptr[k]*fbptr[k];
    }
    cout << "MapMatSparse: b=ax norm:           " << normm <<endl;
    cout << "DegMatSparse valarray: b=ax norm:  " << normd <<endl;
    cout << "DegMatSparse carray: b=ax norm:    " << normp <<endl;
    cout << "FastMatSparse carray: b=ax norm:    " << normf <<endl;

    return 0;
}
