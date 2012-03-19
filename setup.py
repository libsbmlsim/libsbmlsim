#!/usr/bin/env python

"""
setup.py file for SWIG libsbmlsim
"""

from distutils.core import setup, Extension

libsbmlsim_module = Extension('_libsbmlsim', sources=['libsbmlsim_wrap.c'])
setup(name = 'libsbmlsim',
    version = '0.1',
    author = "Akito Tabira",
    description = """Takky""",
    ext_modules = [libsbmlsim_module],
    py_modules = ["libsbmlsim"],
    )
