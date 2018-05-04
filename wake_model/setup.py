#!/usr/bin/env python

"""
setup.py file for SWIG module
"""

from distutils.core import setup, Extension
import numpy
#import os
#os.environ["CC"] = "g++-4.7"
#os.environ["CXX"] = "g++-4.7"

vwake_module = Extension('_vwake',
                           #extra_compile_args=["--fsanitize=leak"],
                           libraries = ['gsl',"gslcblas"],
                           include_dirs = ['/usr/local/lib/'],
                           sources=['vwake.i','vwake.cpp'],
                           swig_opts=["-c++"]
                           )

setup (name = 'vwake',
       version = '0.1',
       author      = "Dagan Pielstick",
       description = """C++ Code for VAWT Wake Model""",
       ext_modules = [vwake_module],
       include_dirs=[numpy.get_include()],
       py_modules = ["vwake"],
       )
