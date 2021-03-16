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


%extend DegMatSparse {

%{
#include <iostream>
%}

/*
PyObject* getValues() {
   int *dims = new int();
   double *data = new double();
   data = &(*(*self).values)[0];
   int ndims = 1;
   dims[0] = (*(*self).values).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
   return arr;
}

PyObject* getCol() {
   int *dims = new int();
   int *data = new int();
   data = &(*(*self).col)[0];
   int ndims = 1;
   dims[0] = (*(*self).col).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_INT, (char *)data);
   return arr;
}

PyObject* getRow() {
   int *dims = new int();
   int *data = new int();
   data = &(*(*self).row)[0];
   int ndims = 1;
   dims[0] = (*(*self).row).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_INT, (char *)data);
   return arr;
}




void setValues(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).values = new valarray<double>((double *)(*arr).data, (*arr).dimensions[0]);
}

void setCol(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).col = new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
}

void setRow(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   int size = (*arr).dimensions[0];
   (*self).row = new valarray<int>((int *)(*arr).data, size);
//   cout << "row.size() = " << (*(*self).row).size()<< endl;
//   cout << "size = " << size << endl;
}

void setm2i(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).m2i = new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
//   cout << "m2i.size() = " << (*(*self).m2i).size()<< endl;
}

PyObject* __mul__(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   valarray<double> *x = new valarray<double>((double *)(*arr).data, (*arr).dimensions[0]);
   valarray<double> *b = new valarray<double>((*x).size());
   (*self).prod(*x, *b);

   int *dims = new int();
   double *data = new double();
   data = &(*b)[0];
   int ndims = 1;
   dims[0] = (*b).size();
   PyObject * res = PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
   return res;
    
}

//PyObject* bmx(PyObject *obj1, PyObject *obj2) {
void  bmx(PyObject *obj1, PyObject *obj2) {
   PyArrayObject *xa = (PyArrayObject *)obj1;
   PyArrayObject *ba = (PyArrayObject *)obj2;
   valarray<double> *x = new valarray<double>((double *)(*xa).data, (*xa).dimensions[0]);
   valarray<double> *b = new valarray<double>((double *)(*ba).data, (*ba).dimensions[0]);
   (*self).prod(*x, *b);

   int *dims = new int();
   double *data = new double();
   data = &(*b)[0];
   int ndims = 1;
   dims[0] = (*b).size();
   obj2 = (PyObject *)PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
}


PyObject* getValues() {
   int *dims = new int();
   double *data = new double();
   data = &(*(*self).values)[0];
   int ndims = 1;
   dims[0] = (*(*self).values).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
   return arr;
}

PyObject* getCol() {
   int *dims = new int();
   int *data = new int();
   data = &(*(*self).col)[0];
   int ndims = 1;
   dims[0] = (*(*self).col).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_INT, (char *)data);
   return arr;
}

PyObject* getRow() {
   int *dims = new int();
   int *data = new int();
   data = &(*(*self).row)[0];
   int ndims = 1;
   dims[0] = (*(*self).row).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_INT, (char *)data);
   return arr;
}




void setValues(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).values = new valarray<double>((double *)(*arr).data, (*arr).dimensions[0]);
}

void setCol(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).col = new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
}

void setRow(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   int size = (*arr).dimensions[0];
   (*self).row = new valarray<int>((int *)(*arr).data, size);
//   cout << "row.size() = " << (*(*self).row).size()<< endl;
//   cout << "size = " << size << endl;
}

void setm2i(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).m2i = new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
//   cout << "m2i.size() = " << (*(*self).m2i).size()<< endl;
}

PyObject* __mul__(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   valarray<double> *x = new valarray<double>((double *)(*arr).data, (*arr).dimensions[0]);
   valarray<double> *b = new valarray<double>((*x).size());
   (*self).prod(*x, *b);

   int *dims = new int();
   double *data = new double();
   data = &(*b)[0];
   int ndims = 1;
   dims[0] = (*b).size();
   PyObject * res = PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
   return res;
    
}

//PyObject* bmx(PyObject *obj1, PyObject *obj2) {
void  bmx(PyObject *obj1, PyObject *obj2) {
   PyArrayObject *xa = (PyArrayObject *)obj1;
   PyArrayObject *ba = (PyArrayObject *)obj2;
   valarray<double> *x = new valarray<double>((double *)(*xa).data, (*xa).dimensions[0]);
   valarray<double> *b = new valarray<double>((double *)(*ba).data, (*ba).dimensions[0]);
   (*self).prod(*x, *b);

   int *dims = new int();
   double *data = new double();
   data = &(*b)[0];
   int ndims = 1;
   dims[0] = (*b).size();
   obj2 = (PyObject *)PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
}

PyObject* getValues() {
   int *dims = new int();
   double *data = new double();
   data = &(*(*self).values)[0];
   int ndims = 1;
   dims[0] = (*(*self).values).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
   return arr;
}

PyObject* getCol() {
   int *dims = new int();
   int *data = new int();
   data = &(*(*self).col)[0];
   int ndims = 1;
   dims[0] = (*(*self).col).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_INT, (char *)data);
   return arr;
}

PyObject* getRow() {
   int *dims = new int();
   int *data = new int();
   data = &(*(*self).row)[0];
   int ndims = 1;
   dims[0] = (*(*self).row).size();
   PyObject * arr = PyArray_FromDimsAndData(ndims, dims,PyArray_INT, (char *)data);
   return arr;
}




void setValues(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).values = new valarray<double>((double *)(*arr).data, (*arr).dimensions[0]);
}

void setCol(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).col = new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
}

void setRow(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   int size = (*arr).dimensions[0];
   (*self).row = new valarray<int>((int *)(*arr).data, size);
//   cout << "row.size() = " << (*(*self).row).size()<< endl;
//   cout << "size = " << size << endl;
}

void setm2i(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   (*self).m2i = new valarray<int>((int *)(*arr).data, (*arr).dimensions[0]);
//   cout << "m2i.size() = " << (*(*self).m2i).size()<< endl;
}
*/

PyObject* __mul__(PyObject *obj) {
   PyArrayObject *arr = (PyArrayObject *)obj;
   valarray<double> *x = new valarray<double>((double *)(*arr).data, (*arr).dimensions[0]);
   valarray<double> *b = new valarray<double>((*x).size());
   (*self).prod(x, b);

   int *dims = new int();
   double *data = new double();
   data = &(*b)[0];
   int ndims = 1;
   dims[0] = (*b).size();
   PyObject * res = PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
   return res;
    
}

/*
//PyObject* bmx(PyObject *obj1, PyObject *obj2) {
void  bmx(PyObject *obj1, PyObject *obj2) {
   PyArrayObject *xa = (PyArrayObject *)obj1;
   PyArrayObject *ba = (PyArrayObject *)obj2;
   valarray<double> *x = new valarray<double>((double *)(*xa).data, (*xa).dimensions[0]);
   valarray<double> *b = new valarray<double>((double *)(*ba).data, (*ba).dimensions[0]);
   (*self).prod(*x, *b);

   int *dims = new int();
   double *data = new double();
   data = &(*b)[0];
   int ndims = 1;
   dims[0] = (*b).size();
   obj2 = (PyObject *)PyArray_FromDimsAndData(ndims, dims,PyArray_DOUBLE, (char *)data);
}

*/

void setOffset(int i) {
    (*self).offset = i;
}

}
