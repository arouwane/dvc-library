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
    from . import _imageRoutines
else:
    import _imageRoutines

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



def ComputeL2ProjLumpedCoefficients(image, xmin, ymin, zmin, zmax, dx, dy, dz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lsc):
    return _imageRoutines.ComputeL2ProjLumpedCoefficients(image, xmin, ymin, zmin, zmax, dx, dy, dz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lsc)

def ComputeL2ProjLumpedCoefficientsSumFact(image, xmin, ymin, zmin, zmax, dx, dy, dz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lsc):
    return _imageRoutines.ComputeL2ProjLumpedCoefficientsSumFact(image, xmin, ymin, zmin, zmax, dx, dy, dz, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, lsc)

def EvaluateTrilinearInterpolation(image, knotXi, knotEta, knotZeta, x, y, z, v):
    return _imageRoutines.EvaluateTrilinearInterpolation(image, knotXi, knotEta, knotZeta, x, y, z, v)

def EvaluateTrilinearInterpolationAndGradient(image, knotXi, knotEta, knotZeta, x, y, z, v, dvdx, dvdy, dvdz):
    return _imageRoutines.EvaluateTrilinearInterpolationAndGradient(image, knotXi, knotEta, knotZeta, x, y, z, v, dvdx, dvdy, dvdz)

def EvaluateTrilinearInterpolationAndGradientStructured(image, knotXi, knotEta, knotZeta, x, y, z, v, dvdx, dvdy, dvdz):
    return _imageRoutines.EvaluateTrilinearInterpolationAndGradientStructured(image, knotXi, knotEta, knotZeta, x, y, z, v, dvdx, dvdy, dvdz)

def EvaluateCardBspline2(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v):
    return _imageRoutines.EvaluateCardBspline2(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v)

def EvaluateCardBspline3(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v):
    return _imageRoutines.EvaluateCardBspline3(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v)

def EvaluateCardBsplineAndGradient2(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz):
    return _imageRoutines.EvaluateCardBsplineAndGradient2(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz)

def EvaluateCardBsplineAndGradient3(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz):
    return _imageRoutines.EvaluateCardBsplineAndGradient3(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz)

def EvaluateCardBsplineAndGradient2Structured(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz):
    return _imageRoutines.EvaluateCardBsplineAndGradient2Structured(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz)

def EvaluateCardBsplineAndGradient3Structured(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz):
    return _imageRoutines.EvaluateCardBsplineAndGradient3Structured(image, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, x, y, z, v, dvdx, dvdy, dvdz)

def GetMeanImageAndStdOnMesh_TrilinearInterp(f, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz):
    return _imageRoutines.GetMeanImageAndStdOnMesh_TrilinearInterp(f, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz)

def GetMeanImageAndStdOnMesh_CBspline2(f, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz):
    return _imageRoutines.GetMeanImageAndStdOnMesh_CBspline2(f, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz)

def GetMeanImageAndStdOnMesh_CBspline3(f, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz):
    return _imageRoutines.GetMeanImageAndStdOnMesh_CBspline3(f, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz)

def GetMeanImageAndStdOnMesh_thrsh_TrilinearInterp(f, thrsh, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, fip, dfipdx, dfipdy, dfipdz):
    return _imageRoutines.GetMeanImageAndStdOnMesh_thrsh_TrilinearInterp(f, thrsh, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, knotXi, knotEta, knotZeta, nipex, nipey, nipez, xg, yg, zg, fmean, fstd, fip, dfipdx, dfipdy, dfipdz)

def GetC8MeshFromVoxelsTrilinearInterp(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, knotXi, knotEta, knotZeta, elementMask):
    return _imageRoutines.GetC8MeshFromVoxelsTrilinearInterp(image, knotXiImage, knotEtaImage, knotZetaImage, thrsh, knotXi, knotEta, knotZeta, elementMask)

def GetMeanImageAndStdOnFE_Mesh_TrilinearInterp(f, knotXiImage, knotEtaImage, knotZetaImage, e, n, N, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz):
    return _imageRoutines.GetMeanImageAndStdOnFE_Mesh_TrilinearInterp(f, knotXiImage, knotEtaImage, knotZetaImage, e, n, N, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz)

def GetMeanImageAndStdOnFE_Mesh_CBspline3(f, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, e, n, N, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz):
    return _imageRoutines.GetMeanImageAndStdOnFE_Mesh_CBspline3(f, xmin, ymin, zmin, dx, dy, dz, xc, yc, zc, e, n, N, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz)

def GetMeanImageAndStdOnFE_Mesh_L2ProjLumped(cp, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, e, n, N, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz):
    return _imageRoutines.GetMeanImageAndStdOnFE_Mesh_L2ProjLumped(cp, knotXiImage, knotEtaImage, knotZetaImage, deg_xi, deg_eta, deg_zeta, e, n, N, fmean, fstd, dyne, fip, dfipdx, dfipdy, dfipdz)


