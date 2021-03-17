from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
           "cc.pyx",                 # our Cython source
           language="# distutils: language=# distutils: language=c++",             # generate C++ code
      ))
#python3 setup.py build_ext --inplace
