#!/usr/bin/env python

from distutils.core import setup, Extension

setup(name = 'pyfdmcutils', 
      version = '0.3',
      description = 'A very tiny utility for PySE/pyFDM',
      author = '�smund �deg�rd',
      author_email = 'aasmund@simula.no',
      ext_modules =  [ Extension('pyfdmcutils', sources=['pyfdmcutils.c']) ] )
