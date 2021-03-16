#!/usr/bin/env python
# -*- coding: utf-8 -*-

# author: Åsmund Ødegård
#         aasmund@simula.no
#         Simula Research Laboratory

from distutils.core import setup, Extension
from os import chdir,system
import sys
import re

# check for prefix in argv:
regexp = re.compile(r'--prefix=.*')

use_prefix=""
for a in sys.argv:
    if regexp.match(a):
        use_prefix=a


chdir("DegMatSparse/swig")


if 'build' in sys.argv:
    system("python setup.py build")
elif 'install' in sys.argv:
    # build first to get swig-options right
    system("python setup.py build")
    # pass on prefix if given
    system("python setup.py install %s" % (use_prefix))
elif 'clean' in sys.argv:
    system("python setup.py clean")


chdir("../..")

setup(name = 'pyFDM', 
      version = '0.3',
      description = 'A framework for solving PDEs with FDM using Python',
      author = 'Åsmund Ødegård',
      author_email = 'aasmund@simula.no',
      packages=['pyFDM','pyPDE'],
      ext_package='pyFDM',
      ext_modules =  [ Extension('pyfdmcutils', sources=['pyFDM/clib/pyfdmcutils.c']) ] )
