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



%extend DMS {

%{
#include <iostream>
%}

PyObject* __mul__(DoubleVector *vec) {
    int dims[1];
    dims[0] = vec->n;
    PyArrayObject *arr = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_DOUBLE);
    self->prod(vec->data, (double *)arr->data);
    return PyArray_Return(arr);
}

void setOffset(int i) {
    (*self).offset = i;
}

}
