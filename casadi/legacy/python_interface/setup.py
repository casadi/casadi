# -*- coding: utf-8 -*-
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("casadi_python", 
              ["casadi_python.pyx"],
              include_dirs=["../../"],
              libraries=["ipopt_interface","sundials_interface","casadi","sundials_idas","sundials_cvodes","sundials_nvecserial","ipopt","lapack","pthread","m"],
              language='c++')
]

setup(
  name = 'Casadi',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)