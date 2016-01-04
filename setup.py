from setuptools import setup
from numpy.distutils.core import setup, Extension


setup(
    name='VAWT_Wake_Model',
    version='1.0',
    description='Parameterized wake model for a vertical-axis wind turbine',
    author='Eric B. Tingey',
    author_email='ebtingey@byu.edu',
    package_dir={'': 'wake_model'},
    py_modules=['VAWT_Wake_Model'],
    # test_suite='test.test.py',
    license='MIT License'
)