/*
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

*/


%include exception.i
%include valarray_typemaps.i

/*
The __setitem__ and __getitem__ special methods translates the function
arguments to tuples. This tuple must be parsed in order to access the
two-dimentional indexing in the external code. Hence, by forcing the indexing
variable names to idx0 and idx1, e.g.,
double SomeClass:: __getitem__(int idx0, int idx1)
we can make this work.
*/
%typemap(in) (int idx0, int idx1) {
    $1 = PyInt_AsLong(PyTuple_GetItem($input,0));
    $2 = PyInt_AsLong(PyTuple_GetItem($input,1));
}

/*
Our basic communication unit is Numeric arrays, and the C++-code accepts arrays
as double pointers. If the name of the double array is x, and the preceding argument 
is an integer, get the info from the PyArrayObject.
*/
%typemap(in) (int n, double *x) {
    if PyArray_Check($input) {
        PyArrayObject *xa = (PyArrayObject *)($input);
        if (xa->descr->type == 'd') {
            $1 = (*xa).dimensions[0];
            $2 = (double *)(*xa).data;
        } else {
            SWIG_exception(SWIG_TypeError, "Array of doubles expected");
        }
    } else {
        SWIG_exception(SWIG_TypeError, "Array expected");
    }
}
/*
This is really a work around a missing SWIG feature (multiple argument typechecks 
for overloaded methods). The swig extension  code of FastMatSparse uses DoubleVector 
as input argument for the overloaded methods. We need to access the vector in
the entire wrapper function, and therefore store it in a temporary variable dv.
*/

%typemap(in) DoubleVector * (DoubleVector dv) {
    if PyArray_Check($input) {
        PyArrayObject *xa = (PyArrayObject *)($input);
        if (xa->descr->type == 'd') {
            dv.n =  (*xa).dimensions[0];
            dv.data = (double *)(*xa).data;
            $1 = &dv;
        } else {
            SWIG_exception(SWIG_TypeError, "Array of doubles expected");
        }
    } else {
        SWIG_exception(SWIG_TypeError, "Array expected");
    }
}

/*

*/
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) DoubleVector * {
    $1 = PyArray_Check($input) ? 1 : 0;
}

%typemap(in) (double *) {
    if PyArray_Check($input) {
        PyArrayObject *xa = (PyArrayObject *)($input);
        if (xa->descr->type == 'd') {
            $1  = (double *)(*xa).data;
        } else {
            SWIG_exception(SWIG_TypeError, "Array of doubles expected");
        }
    } else {
        SWIG_exception(SWIG_TypeError, "Array expected");
    }
}

%typemap(in) (int *) {
    if PyArray_Check($input) {
        PyArrayObject *xa = (PyArrayObject *)($input);
        if (xa->descr->type == 'i') {
            $1  = (int *)(*xa).data;
        } else {
            SWIG_exception(SWIG_TypeError, "Array of integers expected");
        }
    } else {
        SWIG_exception(SWIG_TypeError, "Array expected");
    }
}

