# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
This file is required for the custom build process.
The main configuration is in pyproject.toml, but we still need this for the custom install command.
"""

from setuptools import setup
from setuptools.command.install import install
import subprocess


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        cmd = f'make -j'
        subprocess.check_call(cmd, shell=True)
        super().run()


# This setup() function is called by setuptools from pyproject.toml
setup(
    data_files=[('bin', ['bin/merced'])],
    cmdclass={'install': CustomInstall}
)
