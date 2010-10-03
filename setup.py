#!/usr/bin/env python

import os
from distutils.core import setup, Extension
from numpy import get_include

spext = Extension('_spectractor',
    define_macros = [('MAJOR_VERSION', '0'),
                     ('MINOR_VERSION', '1')],
    include_dirs = [get_include()],
    libraries = ['m'],
    extra_compile_args = ["-O", "-march=native"], # '-Wall', '-g'
    sources = [os.path.join('src', 'spectractor.c')])

description = """Spectractor is a library for manipulating one-dimentional
astronomical spectra, containing number of functions in pure python.
"""

setup(name='spectractor',
    version='0.1.0',
    description=description,
    author='Dmitry Nasonov',
    author_email='gvardopo4ta@gmail.com',
    url='http://spectractor.sourceforge.net/',
    py_modules=['spectractor'],
    ext_modules = [spext]
    )

