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


%newobject *::__mul__;

%extend FastMatSparse {

%{
#include <iostream>
%}

FastMatSparse* __mul__(double d) {
    int n = self->n;
    int nnz = self->nnz;
    FastMatSparse *nmat = new FastMatSparse();
    nmat->n = n;
    nmat->nnz = nnz;
    nmat->offset = self->offset;
    nmat->col = (int *)malloc(nnz*sizeof(int));
    nmat->row = (int *)malloc((n+1)*sizeof(int));
    nmat->values = (double *)malloc(nnz*sizeof(double));
    for (int i=0; i<nnz; i++) {
        if (i <= n) nmat->row[i] = self->row[i];
        nmat->col[i] = self->col[i];
        nmat->values[i] = d*self->values[i];
    }
    nmat->row[nnz] = self->row[nnz];
    return nmat;
}

FastMatSparse* __rmul__(double d) {
    return FastMatSparse___mul____SWIG_0(self, d);
}

PyObject* __mul__(DoubleVector *vec) {
    int dims[1];
    dims[0] = vec->n;
    PyArrayObject *arr = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_DOUBLE);
    self->prod(vec->data, (double *)arr->data);
    return PyArray_Return(arr);
}


double __getitem__(int i) {
    return self->values[i];
}

void setOffset(int i) {
    self->offset = i;
}

}
