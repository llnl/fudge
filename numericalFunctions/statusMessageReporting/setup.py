# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
# This file is required for the custom build process
# The main configuration is in pyproject.toml, but we still need this for the custom install command

import os
import sys
from setuptools import setup
from setuptools.command.install import install
import subprocess

    
def compile_with_make():
    cmd = "make -j"
    subprocess.check_call(cmd, shell=True)


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_with_make()
        super().run()


# This setup() function is called by setuptools from pyproject.toml
setup(
    cmdclass={'install': CustomInstall}
)