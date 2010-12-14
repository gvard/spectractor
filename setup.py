#!/usr/bin/env python

"""Setup script for Spectractor"""

import os
from distutils.core import setup, Extension
from distutils.util import get_platform
from numpy import get_include

description = """Spectractor is a library for manipulating
astronomical spectra images and one-dimentional extracted spectra,
containing number of functions and classes in pure python.
"""

if get_platform()[:3] == 'win':
    spext = None
    #midworker_mod = []
else:
    spext = [Extension('_spectractor',
        define_macros = [('MAJOR_VERSION', '0'),
                         ('MINOR_VERSION', '1')],
        include_dirs = [get_include()],
        libraries = ['m'],
        extra_compile_args = ["-O", "-march=native"], # '-Wall', '-g'
        sources = [os.path.join('src', 'spectractor.c')]
        )]
    #midworker_mod = ["midworker"]

setup(name='spectractor',
    version='0.1.0',
    description=description,
    author='Dmitry Nasonov',
    author_email='gvardopo4ta@gmail.com',
    url='http://spectractor.sourceforge.net/',
    package_dir={'spectractor': ''},
    packages=['spectractor'],
    scripts=[os.path.join('scripts', 'liner.py'),
             os.path.join('scripts', 'preparer.py')],
    ext_modules=spext
    #py_modules=['spectractor', 'serve'] + midworker_mod,
    )
