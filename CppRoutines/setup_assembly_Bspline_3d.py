from distutils.core import setup, Extension
import distutils.sysconfig
import numpy

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-O2", "")

for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-DNDEBUG", "")

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")
        
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


assembly_Bspline_3d_module = Extension('_assembly_Bspline_3d', 
                libraries=['CGAL','gmp'],
                library_dirs=['/usr/lib/x86_64-linux-gnu','/usr/lib/x86_64-linux-gnu'],
                sources=['assembly_Bspline_3d.cpp','tools.cpp','octreeDecomposition.cpp','meshUtils.cpp','memory.cpp','toolsRoutines.cpp',
                         'assembly_Bspline_3d.i'],
				extra_compile_args = ["-g","-frounding-math","-ffast-math","-std=c++11","-fopenmp"], 
				include_dirs = [numpy_include],
                swig_opts=['-c++'])

setup(name='assembly_Bspline_3d', ext_modules=[assembly_Bspline_3d_module], py_modules=["assembly_Bspline_3d"])         