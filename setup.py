#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup, find_packages

setup(
    name='vawtwake',
    version='2.3.0',
    description='Parameterized wake model for a vertical-axis wind turbine',
    author='Eric B. Tingey',
    author_email='ebtingey@byu.edu',
    url='https://github.com/byuflowlab/vawt-wake-model',
    package_dir={'': 'wake_model'},
    py_modules=['vawtwake'],
    test_suite='test.test.py',
    license='MIT License',
    zip_safe=False
)

from numpy.distutils.core import setup, Extension
setup(
    name='vawtwake',
    version='2.3.0',
    package_dir={'': 'wake_model'},
    ext_modules=[Extension('_vawtwake', ['wake_model/VAWT_Wake_Model.f90'], extra_compile_args=['-O2'])],
)
