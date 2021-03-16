# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _MatSparse

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ DoubleVector instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["data"] = _MatSparse.DoubleVector_data_set
    __swig_getmethods__["data"] = _MatSparse.DoubleVector_data_get
    if _newclass:data = property(_MatSparse.DoubleVector_data_get, _MatSparse.DoubleVector_data_set)
    __swig_setmethods__["n"] = _MatSparse.DoubleVector_n_set
    __swig_getmethods__["n"] = _MatSparse.DoubleVector_n_get
    if _newclass:n = property(_MatSparse.DoubleVector_n_get, _MatSparse.DoubleVector_n_set)
    def __init__(self, *args):
        """__init__(self) -> DoubleVector"""
        _swig_setattr(self, DoubleVector, 'this', _MatSparse.new_DoubleVector(*args))
        _swig_setattr(self, DoubleVector, 'thisown', 1)
    def __del__(self, destroy=_MatSparse.delete_DoubleVector):
        """__del__(self)"""
        try:
            if self.thisown: destroy(self)
        except: pass


class DoubleVectorPtr(DoubleVector):
    def __init__(self, this):
        _swig_setattr(self, DoubleVector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DoubleVector, 'thisown', 0)
        _swig_setattr(self, DoubleVector,self.__class__,DoubleVector)
_MatSparse.DoubleVector_swigregister(DoubleVectorPtr)

class DegMatSparse(_object):
    """Proxy of C++ DegMatSparse class"""
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DegMatSparse, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DegMatSparse, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ DegMatSparse instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["col"] = _MatSparse.DegMatSparse_col_set
    __swig_getmethods__["col"] = _MatSparse.DegMatSparse_col_get
    if _newclass:col = property(_MatSparse.DegMatSparse_col_get, _MatSparse.DegMatSparse_col_set)
    __swig_setmethods__["row"] = _MatSparse.DegMatSparse_row_set
    __swig_getmethods__["row"] = _MatSparse.DegMatSparse_row_get
    if _newclass:row = property(_MatSparse.DegMatSparse_row_get, _MatSparse.DegMatSparse_row_set)
    __swig_setmethods__["m2i"] = _MatSparse.DegMatSparse_m2i_set
    __swig_getmethods__["m2i"] = _MatSparse.DegMatSparse_m2i_get
    if _newclass:m2i = property(_MatSparse.DegMatSparse_m2i_get, _MatSparse.DegMatSparse_m2i_set)
    __swig_setmethods__["values"] = _MatSparse.DegMatSparse_values_set
    __swig_getmethods__["values"] = _MatSparse.DegMatSparse_values_get
    if _newclass:values = property(_MatSparse.DegMatSparse_values_get, _MatSparse.DegMatSparse_values_set)
    __swig_setmethods__["offset"] = _MatSparse.DegMatSparse_offset_set
    __swig_getmethods__["offset"] = _MatSparse.DegMatSparse_offset_get
    if _newclass:offset = property(_MatSparse.DegMatSparse_offset_get, _MatSparse.DegMatSparse_offset_set)
    def __init__(self, *args):
        """
        __init__(self) -> DegMatSparse
        __init__(self, valarray<(int)> _col, valarray<(int)> _row, valarray<(int)> _m2i, 
            valarray<(double)> _values) -> DegMatSparse
        __init__(self, MapMatSparse mmat) -> DegMatSparse
        """
        _swig_setattr(self, DegMatSparse, 'this', _MatSparse.new_DegMatSparse(*args))
        _swig_setattr(self, DegMatSparse, 'thisown', 1)
    def __del__(self, destroy=_MatSparse.delete_DegMatSparse):
        """__del__(self)"""
        try:
            if self.thisown: destroy(self)
        except: pass

    def prod(*args):
        """prod(self, valarray<(double)> x, valarray<(double)> b)"""
        return _MatSparse.DegMatSparse_prod(*args)

    def prod2(*args):
        """prod2(self, double x, double b)"""
        return _MatSparse.DegMatSparse_prod2(*args)

    def __call__(*args):
        """__call__(self, int i, int j) -> double"""
        return _MatSparse.DegMatSparse___call__(*args)

    def numberOfRows(*args):
        """numberOfRows(self) -> int"""
        return _MatSparse.DegMatSparse_numberOfRows(*args)

    def __mul__(*args):
        """__mul__(self, PyObject obj) -> PyObject"""
        return _MatSparse.DegMatSparse___mul__(*args)

    def setOffset(*args):
        """setOffset(self, int i)"""
        return _MatSparse.DegMatSparse_setOffset(*args)


class DegMatSparsePtr(DegMatSparse):
    def __init__(self, this):
        _swig_setattr(self, DegMatSparse, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DegMatSparse, 'thisown', 0)
        _swig_setattr(self, DegMatSparse,self.__class__,DegMatSparse)
_MatSparse.DegMatSparse_swigregister(DegMatSparsePtr)

class DMS(_object):
    """Proxy of C++ DMS class"""
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DMS, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DMS, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ DMS instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["n"] = _MatSparse.DMS_n_set
    __swig_getmethods__["n"] = _MatSparse.DMS_n_get
    if _newclass:n = property(_MatSparse.DMS_n_get, _MatSparse.DMS_n_set)
    __swig_setmethods__["nnz"] = _MatSparse.DMS_nnz_set
    __swig_getmethods__["nnz"] = _MatSparse.DMS_nnz_get
    if _newclass:nnz = property(_MatSparse.DMS_nnz_get, _MatSparse.DMS_nnz_set)
    __swig_setmethods__["nnzr"] = _MatSparse.DMS_nnzr_set
    __swig_getmethods__["nnzr"] = _MatSparse.DMS_nnzr_get
    if _newclass:nnzr = property(_MatSparse.DMS_nnzr_get, _MatSparse.DMS_nnzr_set)
    __swig_setmethods__["col"] = _MatSparse.DMS_col_set
    __swig_getmethods__["col"] = _MatSparse.DMS_col_get
    if _newclass:col = property(_MatSparse.DMS_col_get, _MatSparse.DMS_col_set)
    __swig_setmethods__["row"] = _MatSparse.DMS_row_set
    __swig_getmethods__["row"] = _MatSparse.DMS_row_get
    if _newclass:row = property(_MatSparse.DMS_row_get, _MatSparse.DMS_row_set)
    __swig_setmethods__["m2i"] = _MatSparse.DMS_m2i_set
    __swig_getmethods__["m2i"] = _MatSparse.DMS_m2i_get
    if _newclass:m2i = property(_MatSparse.DMS_m2i_get, _MatSparse.DMS_m2i_set)
    __swig_setmethods__["values"] = _MatSparse.DMS_values_set
    __swig_getmethods__["values"] = _MatSparse.DMS_values_get
    if _newclass:values = property(_MatSparse.DMS_values_get, _MatSparse.DMS_values_set)
    __swig_setmethods__["offset"] = _MatSparse.DMS_offset_set
    __swig_getmethods__["offset"] = _MatSparse.DMS_offset_get
    if _newclass:offset = property(_MatSparse.DMS_offset_get, _MatSparse.DMS_offset_set)
    def __init__(self, *args):
        """
        __init__(self) -> DMS
        __init__(self, int _n, int _nnz, int __nnzr) -> DMS
        __init__(self, int _n, int _nnz, int _nnzr, int _col, int _row, int _m2i, 
            double _values) -> DMS
        __init__(self, MapMatSparse mmat) -> DMS
        """
        _swig_setattr(self, DMS, 'this', _MatSparse.new_DMS(*args))
        _swig_setattr(self, DMS, 'thisown', 1)
    def __del__(self, destroy=_MatSparse.delete_DMS):
        """__del__(self)"""
        try:
            if self.thisown: destroy(self)
        except: pass

    def prod(*args):
        """prod(self, double x, double b)"""
        return _MatSparse.DMS_prod(*args)

    def __call__(*args):
        """__call__(self, int i, int j) -> double"""
        return _MatSparse.DMS___call__(*args)

    def numberOfRows(*args):
        """numberOfRows(self) -> int"""
        return _MatSparse.DMS_numberOfRows(*args)

    def __mul__(*args):
        """__mul__(self, DoubleVector vec) -> PyObject"""
        return _MatSparse.DMS___mul__(*args)

    def setOffset(*args):
        """setOffset(self, int i)"""
        return _MatSparse.DMS_setOffset(*args)


class DMSPtr(DMS):
    def __init__(self, this):
        _swig_setattr(self, DMS, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DMS, 'thisown', 0)
        _swig_setattr(self, DMS,self.__class__,DMS)
_MatSparse.DMS_swigregister(DMSPtr)

class MapMatSparse(_object):
    """Proxy of C++ MapMatSparse class"""
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MapMatSparse, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MapMatSparse, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ MapMatSparse instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["values"] = _MatSparse.MapMatSparse_values_set
    __swig_getmethods__["values"] = _MatSparse.MapMatSparse_values_get
    if _newclass:values = property(_MatSparse.MapMatSparse_values_get, _MatSparse.MapMatSparse_values_set)
    def __init__(self, *args):
        """__init__(self) -> MapMatSparse"""
        _swig_setattr(self, MapMatSparse, 'this', _MatSparse.new_MapMatSparse(*args))
        _swig_setattr(self, MapMatSparse, 'thisown', 1)
    def __del__(self, destroy=_MatSparse.delete_MapMatSparse):
        """__del__(self)"""
        try:
            if self.thisown: destroy(self)
        except: pass

    def __call__(*args):
        """__call__(self, int i, int j) -> double"""
        return _MatSparse.MapMatSparse___call__(*args)

    def prod(*args):
        """prod(self, valarray<(double)> x, valarray<(double)> y)"""
        return _MatSparse.MapMatSparse_prod(*args)

    def save(*args):
        """save(self, string filename, string name)"""
        return _MatSparse.MapMatSparse_save(*args)

    def load(*args):
        """load(self, char filename)"""
        return _MatSparse.MapMatSparse_load(*args)

    def __setitem__(*args):
        """__setitem__(self, int idx0, double d)"""
        return _MatSparse.MapMatSparse___setitem__(*args)

    def __getitem__(*args):
        """__getitem__(self, int idx0) -> double"""
        return _MatSparse.MapMatSparse___getitem__(*args)


class MapMatSparsePtr(MapMatSparse):
    def __init__(self, this):
        _swig_setattr(self, MapMatSparse, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, MapMatSparse, 'thisown', 0)
        _swig_setattr(self, MapMatSparse,self.__class__,MapMatSparse)
_MatSparse.MapMatSparse_swigregister(MapMatSparsePtr)

class FastMatSparse(_object):
    """Proxy of C++ FastMatSparse class"""
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FastMatSparse, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FastMatSparse, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ FastMatSparse instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["n"] = _MatSparse.FastMatSparse_n_set
    __swig_getmethods__["n"] = _MatSparse.FastMatSparse_n_get
    if _newclass:n = property(_MatSparse.FastMatSparse_n_get, _MatSparse.FastMatSparse_n_set)
    __swig_setmethods__["nnz"] = _MatSparse.FastMatSparse_nnz_set
    __swig_getmethods__["nnz"] = _MatSparse.FastMatSparse_nnz_get
    if _newclass:nnz = property(_MatSparse.FastMatSparse_nnz_get, _MatSparse.FastMatSparse_nnz_set)
    __swig_setmethods__["col"] = _MatSparse.FastMatSparse_col_set
    __swig_getmethods__["col"] = _MatSparse.FastMatSparse_col_get
    if _newclass:col = property(_MatSparse.FastMatSparse_col_get, _MatSparse.FastMatSparse_col_set)
    __swig_setmethods__["row"] = _MatSparse.FastMatSparse_row_set
    __swig_getmethods__["row"] = _MatSparse.FastMatSparse_row_get
    if _newclass:row = property(_MatSparse.FastMatSparse_row_get, _MatSparse.FastMatSparse_row_set)
    __swig_setmethods__["values"] = _MatSparse.FastMatSparse_values_set
    __swig_getmethods__["values"] = _MatSparse.FastMatSparse_values_get
    if _newclass:values = property(_MatSparse.FastMatSparse_values_get, _MatSparse.FastMatSparse_values_set)
    __swig_setmethods__["offset"] = _MatSparse.FastMatSparse_offset_set
    __swig_getmethods__["offset"] = _MatSparse.FastMatSparse_offset_get
    if _newclass:offset = property(_MatSparse.FastMatSparse_offset_get, _MatSparse.FastMatSparse_offset_set)
    def __del__(self, destroy=_MatSparse.delete_FastMatSparse):
        """__del__(self)"""
        try:
            if self.thisown: destroy(self)
        except: pass

    def __init__(self, *args):
        """
        __init__(self) -> FastMatSparse
        __init__(self, int n, int nnz) -> FastMatSparse
        __init__(self, int n, int nnz, int _col, int _row, double _values) -> FastMatSparse
        __init__(self, MapMatSparse mmat) -> FastMatSparse
        """
        _swig_setattr(self, FastMatSparse, 'this', _MatSparse.new_FastMatSparse(*args))
        _swig_setattr(self, FastMatSparse, 'thisown', 1)
    def prod(*args):
        """prod(self, double xptr, double bptr)"""
        return _MatSparse.FastMatSparse_prod(*args)

    def __call__(*args):
        """__call__(self, int i, int j) -> double"""
        return _MatSparse.FastMatSparse___call__(*args)

    def numberOfRows(*args):
        """numberOfRows(self) -> int"""
        return _MatSparse.FastMatSparse_numberOfRows(*args)

    def numberOfNonZeros(*args):
        """numberOfNonZeros(self) -> int"""
        return _MatSparse.FastMatSparse_numberOfNonZeros(*args)

    def __rmul__(*args):
        """__rmul__(self, double d) -> FastMatSparse"""
        return _MatSparse.FastMatSparse___rmul__(*args)

    def __mul__(*args):
        """
        __mul__(self, double d) -> FastMatSparse
        __mul__(self, DoubleVector vec) -> PyObject
        """
        return _MatSparse.FastMatSparse___mul__(*args)

    def __getitem__(*args):
        """__getitem__(self, int i) -> double"""
        return _MatSparse.FastMatSparse___getitem__(*args)

    def setOffset(*args):
        """setOffset(self, int i)"""
        return _MatSparse.FastMatSparse_setOffset(*args)


class FastMatSparsePtr(FastMatSparse):
    def __init__(self, this):
        _swig_setattr(self, FastMatSparse, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FastMatSparse, 'thisown', 0)
        _swig_setattr(self, FastMatSparse,self.__class__,FastMatSparse)
_MatSparse.FastMatSparse_swigregister(FastMatSparsePtr)


