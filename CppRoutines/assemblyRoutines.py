# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _assemblyRoutines
else:
    import _assemblyRoutines

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class SwigPyIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _assemblyRoutines.delete_SwigPyIterator

    def value(self):
        return _assemblyRoutines.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _assemblyRoutines.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _assemblyRoutines.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _assemblyRoutines.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _assemblyRoutines.SwigPyIterator_equal(self, x)

    def copy(self):
        return _assemblyRoutines.SwigPyIterator_copy(self)

    def next(self):
        return _assemblyRoutines.SwigPyIterator_next(self)

    def __next__(self):
        return _assemblyRoutines.SwigPyIterator___next__(self)

    def previous(self):
        return _assemblyRoutines.SwigPyIterator_previous(self)

    def advance(self, n):
        return _assemblyRoutines.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _assemblyRoutines.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _assemblyRoutines.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _assemblyRoutines.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _assemblyRoutines.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _assemblyRoutines.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _assemblyRoutines.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self

# Register SwigPyIterator in _assemblyRoutines:
_assemblyRoutines.SwigPyIterator_swigregister(SwigPyIterator)

class VectDouble(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _assemblyRoutines.VectDouble_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _assemblyRoutines.VectDouble___nonzero__(self)

    def __bool__(self):
        return _assemblyRoutines.VectDouble___bool__(self)

    def __len__(self):
        return _assemblyRoutines.VectDouble___len__(self)

    def __getslice__(self, i, j):
        return _assemblyRoutines.VectDouble___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _assemblyRoutines.VectDouble___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _assemblyRoutines.VectDouble___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _assemblyRoutines.VectDouble___delitem__(self, *args)

    def __getitem__(self, *args):
        return _assemblyRoutines.VectDouble___getitem__(self, *args)

    def __setitem__(self, *args):
        return _assemblyRoutines.VectDouble___setitem__(self, *args)

    def pop(self):
        return _assemblyRoutines.VectDouble_pop(self)

    def append(self, x):
        return _assemblyRoutines.VectDouble_append(self, x)

    def empty(self):
        return _assemblyRoutines.VectDouble_empty(self)

    def size(self):
        return _assemblyRoutines.VectDouble_size(self)

    def swap(self, v):
        return _assemblyRoutines.VectDouble_swap(self, v)

    def begin(self):
        return _assemblyRoutines.VectDouble_begin(self)

    def end(self):
        return _assemblyRoutines.VectDouble_end(self)

    def rbegin(self):
        return _assemblyRoutines.VectDouble_rbegin(self)

    def rend(self):
        return _assemblyRoutines.VectDouble_rend(self)

    def clear(self):
        return _assemblyRoutines.VectDouble_clear(self)

    def get_allocator(self):
        return _assemblyRoutines.VectDouble_get_allocator(self)

    def pop_back(self):
        return _assemblyRoutines.VectDouble_pop_back(self)

    def erase(self, *args):
        return _assemblyRoutines.VectDouble_erase(self, *args)

    def __init__(self, *args):
        _assemblyRoutines.VectDouble_swiginit(self, _assemblyRoutines.new_VectDouble(*args))

    def push_back(self, x):
        return _assemblyRoutines.VectDouble_push_back(self, x)

    def front(self):
        return _assemblyRoutines.VectDouble_front(self)

    def back(self):
        return _assemblyRoutines.VectDouble_back(self)

    def assign(self, n, x):
        return _assemblyRoutines.VectDouble_assign(self, n, x)

    def resize(self, *args):
        return _assemblyRoutines.VectDouble_resize(self, *args)

    def insert(self, *args):
        return _assemblyRoutines.VectDouble_insert(self, *args)

    def reserve(self, n):
        return _assemblyRoutines.VectDouble_reserve(self, n)

    def capacity(self):
        return _assemblyRoutines.VectDouble_capacity(self)
    __swig_destroy__ = _assemblyRoutines.delete_VectDouble

# Register VectDouble in _assemblyRoutines:
_assemblyRoutines.VectDouble_swigregister(VectDouble)

class VectInt(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _assemblyRoutines.VectInt_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _assemblyRoutines.VectInt___nonzero__(self)

    def __bool__(self):
        return _assemblyRoutines.VectInt___bool__(self)

    def __len__(self):
        return _assemblyRoutines.VectInt___len__(self)

    def __getslice__(self, i, j):
        return _assemblyRoutines.VectInt___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _assemblyRoutines.VectInt___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _assemblyRoutines.VectInt___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _assemblyRoutines.VectInt___delitem__(self, *args)

    def __getitem__(self, *args):
        return _assemblyRoutines.VectInt___getitem__(self, *args)

    def __setitem__(self, *args):
        return _assemblyRoutines.VectInt___setitem__(self, *args)

    def pop(self):
        return _assemblyRoutines.VectInt_pop(self)

    def append(self, x):
        return _assemblyRoutines.VectInt_append(self, x)

    def empty(self):
        return _assemblyRoutines.VectInt_empty(self)

    def size(self):
        return _assemblyRoutines.VectInt_size(self)

    def swap(self, v):
        return _assemblyRoutines.VectInt_swap(self, v)

    def begin(self):
        return _assemblyRoutines.VectInt_begin(self)

    def end(self):
        return _assemblyRoutines.VectInt_end(self)

    def rbegin(self):
        return _assemblyRoutines.VectInt_rbegin(self)

    def rend(self):
        return _assemblyRoutines.VectInt_rend(self)

    def clear(self):
        return _assemblyRoutines.VectInt_clear(self)

    def get_allocator(self):
        return _assemblyRoutines.VectInt_get_allocator(self)

    def pop_back(self):
        return _assemblyRoutines.VectInt_pop_back(self)

    def erase(self, *args):
        return _assemblyRoutines.VectInt_erase(self, *args)

    def __init__(self, *args):
        _assemblyRoutines.VectInt_swiginit(self, _assemblyRoutines.new_VectInt(*args))

    def push_back(self, x):
        return _assemblyRoutines.VectInt_push_back(self, x)

    def front(self):
        return _assemblyRoutines.VectInt_front(self)

    def back(self):
        return _assemblyRoutines.VectInt_back(self)

    def assign(self, n, x):
        return _assemblyRoutines.VectInt_assign(self, n, x)

    def resize(self, *args):
        return _assemblyRoutines.VectInt_resize(self, *args)

    def insert(self, *args):
        return _assemblyRoutines.VectInt_insert(self, *args)

    def reserve(self, n):
        return _assemblyRoutines.VectInt_reserve(self, n)

    def capacity(self):
        return _assemblyRoutines.VectInt_capacity(self)
    __swig_destroy__ = _assemblyRoutines.delete_VectInt

# Register VectInt in _assemblyRoutines:
_assemblyRoutines.VectInt_swigregister(VectInt)


def GetBsplineFunctionsMatrixStructured(x, y, z, knotXi, knotEta, knotZeta, deg_xi, deg_eta, deg_zeta, valuesN, indexI, indexJ):
    return _assemblyRoutines.GetBsplineFunctionsMatrixStructured(x, y, z, knotXi, knotEta, knotZeta, deg_xi, deg_eta, deg_zeta, valuesN, indexI, indexJ)

def GetBsplineFunctionsAndDerivativesMatrixStructured(x, y, z, knotXi, knotEta, knotZeta, deg_xi, deg_eta, deg_zeta, valuesN, valuesdNdx, valuesdNdy, valuesdNdz, indexI, indexJ):
    return _assemblyRoutines.GetBsplineFunctionsAndDerivativesMatrixStructured(x, y, z, knotXi, knotEta, knotZeta, deg_xi, deg_eta, deg_zeta, valuesN, valuesdNdx, valuesdNdy, valuesdNdz, indexI, indexJ)

def FcmIntegrationTrilinearInterp(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lvlmax):
    return _assemblyRoutines.FcmIntegrationTrilinearInterp(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lvlmax)

def FcmTrilinearInterpStiffnessAndBfIntegral(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lvlmax, E, nu, indexI, indexJ, nnz_values, intBf):
    return _assemblyRoutines.FcmTrilinearInterpStiffnessAndBfIntegral(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lvlmax, E, nu, indexI, indexJ, nnz_values, intBf)

def FcmTrilinearInterpStiffnessAndBfIntegralParallel(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lvlmax, E, nu, indexI, indexJ, nnz_values, intBf):
    return _assemblyRoutines.FcmTrilinearInterpStiffnessAndBfIntegralParallel(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lvlmax, E, nu, indexI, indexJ, nnz_values, intBf)

def Laplacian_Structured(deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, indexI, indexJ, nnz_values):
    return _assemblyRoutines.Laplacian_Structured(deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, indexI, indexJ, nnz_values)

def Laplacian_Structured_Parallel(deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, indexI, indexJ, nnz_values):
    return _assemblyRoutines.Laplacian_Structured_Parallel(deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, indexI, indexJ, nnz_values)

def L2Projection(deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, deg_xi2, deg_eta2, deg_zeta2, knotXi2, knotEta2, knotZeta2, U2, indexI, indexJ, nnz_values, rhs):
    return _assemblyRoutines.L2Projection(deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, deg_xi2, deg_eta2, deg_zeta2, knotXi2, knotEta2, knotZeta2, U2, indexI, indexJ, nnz_values, rhs)

def VoxelIntegrationThreshold(fip, thrsh, nipex, nipey, nipez, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta):
    return _assemblyRoutines.VoxelIntegrationThreshold(fip, thrsh, nipex, nipey, nipez, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta)

def VoxelIntegrationMask(maskip, nipex, nipey, nipez, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta):
    return _assemblyRoutines.VoxelIntegrationMask(maskip, nipex, nipey, nipez, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta)

def DVC_LHS_thrsh(dfipdx, dfipdy, dfipdz, ipIndices, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, Nxi, Neta, Nzeta, indexI, indexJ, nnz_values):
    return _assemblyRoutines.DVC_LHS_thrsh(dfipdx, dfipdy, dfipdz, ipIndices, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, Nxi, Neta, Nzeta, indexI, indexJ, nnz_values)

def DVC_RHS_thrsh_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, ipIndices, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs):
    return _assemblyRoutines.DVC_RHS_thrsh_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, ipIndices, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs)

def DVC_RHS_ZN_thrsh_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, fmeane, fstde, ipIndices, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs):
    return _assemblyRoutines.DVC_RHS_ZN_thrsh_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, fmeane, fstde, ipIndices, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs)

def DVC_LHS_Structured(dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, Nxi, Neta, Nzeta, indexI, indexJ, nnz_values):
    return _assemblyRoutines.DVC_LHS_Structured(dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, Nxi, Neta, Nzeta, indexI, indexJ, nnz_values)

def DVC_LHS_Structured_Parallel(dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, Nxi, Neta, Nzeta, indexI, indexJ, nnz_values):
    return _assemblyRoutines.DVC_LHS_Structured_Parallel(dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, Nxi, Neta, Nzeta, indexI, indexJ, nnz_values)

def DVC_RHS_Structured_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs):
    return _assemblyRoutines.DVC_RHS_Structured_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs)

def DVC_RHS_Structured_TrilinearInterp_Parallel(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs):
    return _assemblyRoutines.DVC_RHS_Structured_TrilinearInterp_Parallel(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs)

def DVC_RHS_ZN_Structured_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, fmeane, fstde, dyne, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs, elementRes):
    return _assemblyRoutines.DVC_RHS_ZN_Structured_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, fmeane, fstde, dyne, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs, elementRes)

def GLR_Structured_TrilinearInterp(fip, g, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr):
    return _assemblyRoutines.GLR_Structured_TrilinearInterp(fip, g, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr)

def GLR_Structured_TrilinearInterp_Parallel(fip, g, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr):
    return _assemblyRoutines.GLR_Structured_TrilinearInterp_Parallel(fip, g, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr)

def GLR_ZN_Structured_TrilinearInterp(fip, fmeane, fstde, g, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr):
    return _assemblyRoutines.GLR_ZN_Structured_TrilinearInterp(fip, fmeane, fstde, g, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr)

def DVC_RHS_Structured_CBspline2(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs):
    return _assemblyRoutines.DVC_RHS_Structured_CBspline2(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs)

def GLR_Structured_CBspline2(fip, g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr):
    return _assemblyRoutines.GLR_Structured_CBspline2(fip, g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr)

def DVC_RHS_Structured_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs):
    return _assemblyRoutines.DVC_RHS_Structured_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, dfipdx, dfipdy, dfipdz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, NOELEM, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, rhs)

def GLR_Structured_CBspline3(fip, g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr):
    return _assemblyRoutines.GLR_Structured_CBspline3(fip, g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, Nxi, Neta, Nzeta, U, glr)

def DVC_LHS_FE(e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, dfipdx, dfipdy, dfipdz, indexI, indexJ, nnz_values):
    return _assemblyRoutines.DVC_LHS_FE(e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, dfipdx, dfipdy, dfipdz, indexI, indexJ, nnz_values)

def DVC_RHS_FE_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, U, rhs):
    return _assemblyRoutines.DVC_RHS_FE_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, dfipdx, dfipdy, dfipdz, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, U, rhs)

def GLR_FE_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, e, n, conn, N, w, U, glr):
    return _assemblyRoutines.GLR_FE_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, fip, e, n, conn, N, w, U, glr)

def Gophi_FEMesh_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, xi, eta, zeta, ie, e, n, conn, U, glr):
    return _assemblyRoutines.Gophi_FEMesh_TrilinearInterp(g, knotXiImage, knotEtaImage, knotZetaImage, xi, eta, zeta, ie, e, n, conn, U, glr)

def DVC_RHS_FE_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, dfipdx, dfipdy, dfipdz, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, U, rhs):
    return _assemblyRoutines.DVC_RHS_FE_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, dfipdx, dfipdy, dfipdz, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, U, rhs)

def GLR_FE_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, e, n, conn, N, w, U, glr):
    return _assemblyRoutines.GLR_FE_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, fip, e, n, conn, N, w, U, glr)

def Gophi_FEMesh_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, xi, eta, zeta, ie, e, n, conn, U, glr):
    return _assemblyRoutines.Gophi_FEMesh_CBspline3(g, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, xc, yc, zc, xi, eta, zeta, ie, e, n, conn, U, glr)

def DVC_RHS_FE_L2ProjLumped(lsc_g, knotXiImage_g, knotEtaImage_g, knotZetaImage_g, deg_xi_g, deg_eta_g, deg_zeta_g, fip, dfipdx, dfipdy, dfipdz, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, U, rhs):
    return _assemblyRoutines.DVC_RHS_FE_L2ProjLumped(lsc_g, knotXiImage_g, knotEtaImage_g, knotZetaImage_g, deg_xi_g, deg_eta_g, deg_zeta_g, fip, dfipdx, dfipdy, dfipdz, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, U, rhs)

def GLR_FE_L2ProjLumped(lsc_g, knotXiImage_g, knotEtaImage_g, knotZetaImage_g, deg_xi_g, deg_eta_g, deg_zeta_g, fip, e, n, conn, N, w, U, glr):
    return _assemblyRoutines.GLR_FE_L2ProjLumped(lsc_g, knotXiImage_g, knotEtaImage_g, knotZetaImage_g, deg_xi_g, deg_eta_g, deg_zeta_g, fip, e, n, conn, N, w, U, glr)

def Gophi_FEMesh_L2ProjLumped(lsc_g, knotXiImage_g, knotEtaImage_g, knotZetaImage_g, deg_xi_g, deg_eta_g, deg_zeta_g, xi, eta, zeta, ie, e, n, conn, U, glr):
    return _assemblyRoutines.Gophi_FEMesh_L2ProjLumped(lsc_g, knotXiImage_g, knotEtaImage_g, knotZetaImage_g, deg_xi_g, deg_eta_g, deg_zeta_g, xi, eta, zeta, ie, e, n, conn, U, glr)

def Laplacian_FE(e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, indexI, indexJ, nnz_values):
    return _assemblyRoutines.Laplacian_FE(e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, indexI, indexJ, nnz_values)

def Stiffness_FE(E, nu, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, indexI, indexJ, nnz_values):
    return _assemblyRoutines.Stiffness_FE(E, nu, e, n, conn, N, dNdxi, dNdeta, dNdzeta, w, indexI, indexJ, nnz_values)


