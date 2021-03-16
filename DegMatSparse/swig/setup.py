#!/usr/bin/env python
# -*- coding: utf-8 -*-

# author: Ola Skavhaug, Åsmund Ødegård
#
from distutils.core import setup,Extension
import distutils
from os import chdir,system
import os
import sys


# Maybe we should try to detect a C++ compiler?
# For now, if there is an CXX in env, use that one, else
# use g++

# make sure that required objects exist in parent dir:
chdir("..")

if 'clean' in sys.argv:
    system("make clean")
else:
    if os.environ.has_key('CXX'):
        system("make CXX=%s DegMatSparse.o MapMatSparse.o DMS.o FastMatSparse.o isearch.o" % (os.environ['CXX']))
    else:
        system("make DegMatSparse.o MapMatSparse.o DMS.o FastMatSparse.o isearch.o")


# cd back again.
chdir("swig")

# we need to use cpp as extension for swig-generated files.
# Due to a bug in distutils, the option must be set globally, not inside
# the Extension directive

swig_opt = '--swig-cpp'
if distutils.__version__ >= '2.4': swig_opt = '--swig-opts=-c++'

if sys.argv[1] == 'build':
    sys.argv[1] = 'build_ext'
    sys.argv.insert(2, swig_opt)

setup(name='MatSparse',
      version='0.5',
      py_modules=['MatSparse'],
      ext_modules=[Extension('_MatSparse',['MatSparse.i'], 
          swig_opts=['-c++'],
          include_dirs=['../'],
          extra_objects=['../DegMatSparse.o','../MapMatSparse.o','../DMS.o','../FastMatSparse.o','../isearch.o'])],
      author='Ola Skavhaug',
      author_email='skavhaug@simula.no',
      )
