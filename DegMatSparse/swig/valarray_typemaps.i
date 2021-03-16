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


%typemap(in) valarray<double>*{
    PyArrayObject *arr = (PyArrayObject *)($input);
    const double *__restrict__ data = (double *__restrict__)(*arr).data;
    ($1) = new valarray<double>(data, (*arr).dimensions[0]);
}

%typemap(in) valarray<int>&{
    PyArrayObject *arr = (PyArrayObject *)($input);
    ($1) =  new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
}

%typemap(in) valarray<int>*{
    PyArrayObject *arr = (PyArrayObject *)($input);
    ($1) =  new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
}


%typemap(out) valarray<double>* {
   int *dims = new int();
   double *data = new double();
   data = &(*$1)[0];
   int ndims = 1;
   dims[0] = (*$1).size();
   PyObject *obj = new PyObject();
   obj = (PyObject *)PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
   $result = obj;
}

%typemap(out) valarray<int>* {
   int *dims = new int();
   int *data = new int();
   data = &(*$1)[0];
   int ndims = 1;
   dims[0] = (*$1).size();
   PyObject *obj = new PyObject();
   obj = (PyObject *)PyArray_FromDimsAndData(ndims, dims,PyArray_INT, (char *)data);
   $result = obj;
}
