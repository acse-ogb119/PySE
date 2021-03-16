
#include <Python.h>


static PyObject* closestint(PyObject *self, PyObject* args);

static PyMethodDef UtilsMethods[] = {
   { (char *)"closestint", closestint, METH_VARARGS },
   { NULL, NULL }
};

void initpyfdmcutils () {
   Py_InitModule("pyfdmcutils", UtilsMethods);
};
