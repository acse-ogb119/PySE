
#include "pyfdmcutils.h"
#include <math.h>

static PyObject* closestint (PyObject *self, PyObject *args) {

   int f_val,c_val,ret_val;
   double x;

   if (!PyArg_ParseTuple(args, "d", &x))
         return NULL;

   f_val = floor(x);
   c_val = ceil(x);

   ret_val =  (x - f_val > 0.5) ? c_val : f_val;

   return Py_BuildValue("i",ret_val);
}
