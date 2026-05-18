# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
# This file provides compatibility with pyproject.toml
# Main configuration is in pyproject.toml

import glob
from setuptools import setup

# This setup() function is called by setuptools from pyproject.toml
setup(
    scripts=['legacy/bin/prepro.py'] + glob.glob('bin/*.py'),
)