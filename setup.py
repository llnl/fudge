# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# This file is required for C extensions and custom installation
# The main configuration is in pyproject.toml, but we need this for extensions and customizations

import os, sys, glob, shutil
import setuptools
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
from setuptools import setup, Extension
from pathlib import Path
import subprocess
import numpy

cwd: str = Path(__file__).parent

class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        # copy C executables Merced/bin/merced and fudge/processing/deterministic/upscatter/bin/calcUpscatterKernel to Python environment bin folder
        workingFolder = os.getcwd()
        binFolder = os.path.join(sys.prefix, 'bin')
        os.chdir('Merced')
        subprocess.check_call('make -j', shell=True)
        executable = "bin/merced"
        if sys.platform.startswith('win'):
            executable = "bin/merced.exe"
        shutil.copy(executable, binFolder)
        os.chdir(workingFolder)

        os.chdir('fudge/processing/deterministic/upscatter')
        subprocess.check_call('make -j', shell=True)
        executable = "bin/calcUpscatterKernel"
        if sys.platform.startswith('win'):
            executable = "bin/calcUpscatterKernel.exe"
        shutil.copy(executable, binFolder)
        os.chdir(workingFolder)

        super().run()


# This setup() is called by setuptools from pyproject.toml
setup(
    scripts=glob.glob('bin/*.py'),
    ext_modules=[
        Extension('fudge.processing.resonances._getBreitWignerSums',
            sources=['fudge/processing/resonances/getBreitWignerSums.c'], include_dirs=[numpy.get_include()], ),
        Extension('fudge.processing.resonances._getScatteringMatrices',
            sources=['fudge/processing/resonances/getScatteringMatrices.c'], include_dirs=[numpy.get_include()], ),
        Extension('fudge.processing.resonances._getCoulombWavefunctions',
            sources=['fudge/processing/resonances/getCoulombWavefunctions.c', 'fudge/processing/resonances/coulfg2.c'], include_dirs=[numpy.get_include()], ),
    ],
    install_requires=[
        'numpy',
        f'crossSectionAdjustForHeatedTarget @ {(cwd / "crossSectionAdjustForHeatedTarget").as_uri()}',
        f'numericalFunctions @ {(cwd / "numericalFunctions").as_uri()}',
        f'pqu @ {(cwd / "pqu").as_uri()}',
        f'xData @ {(cwd / "xData").as_uri()}',
        f'PoPs @ {(cwd / "PoPs").as_uri()}',
        f'brownies @ {(cwd / "brownies").as_uri()}',
    ],
    cmdclass={'install': CustomInstall}
)
