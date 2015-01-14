'''
Created on 17.04.2013

@author: uwe
'''
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'L2A_AtmCorr.so',
    ext_modules = cythonize('L2A_AtmCorr.py'), # accepts a glob pattern
)

if __name__ == '__main__':
    pass