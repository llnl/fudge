# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
# This file is required for the C extension
# The main configuration is in pyproject.toml, but setuptools still uses this for extensions

from setuptools import setup, Extension
import os
import glob

# Clean up old libraries if they exist
libs = glob.glob(os.path.join('build', 'lib*', 'crossSectionAdjustForHeatedTarget*'))
for lib in libs: 
    try:
        os.remove(lib)
    except:
        pass

# Gather source files
sources = glob.glob(os.path.join('Src', '*.c'))
sources += glob.glob(os.path.join('Python', 'Src', '*.c'))

# Define the extension
crossSectionAdjustForHeatedTarget = Extension(
    'crossSectionAdjustForHeatedTarget.crossSectionAdjustForHeatedTarget', 
    sources=sources, 
    include_dirs=['Src']
)

# This setup() function is called by setuptools from pyproject.toml
setup(
    ext_modules=[crossSectionAdjustForHeatedTarget],
)