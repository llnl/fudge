# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

import fudge as fudgeModule

"""
This module contains methods to find scripts or executables in the FUDGE bin folder.
The location of the bin folder differs if FUDGE was installed with make vs. pip install.
"""

def locateFileInBin( fileInBin, exceptionIfNotFound=False ):
    """
    Find the full path to a file located in the FUDGE bin folder or local directory.

    :param fileInBin:           The desired executable or script.
    :param exceptionIfNotFound: Option to raise a FileNotFoundError exception if the requested path is not found. 
    :return:                    Full path to the requested file path. 
    """

    filePath = os.path.join(os.path.abspath('./'), fileInBin)
    if not os.path.exists( filePath ):
        filePath = os.path.join(os.path.dirname(fudgeModule.__file__), 'bin', fileInBin)
        if not os.path.exists( filePath ):
            filePath = os.path.join(os.path.split(os.path.dirname(fudgeModule.__file__))[0], 'bin', fileInBin)
            if not os.path.exists(filePath):
                if exceptionIfNotFound:
                    raise FileNotFoundError(f'FUDGE bin folder file "{fileInBin}" not found.')
                else:
                    filePath = None

    return filePath

def locateMerced( exceptionIfNotFound=False ):
    """
    Find the full path to the merced executable.

    :param exceptionIfNotFound: Option to raise a FileNotFoundError exception if no "mercury" executable is found.
    :return:                    Full path of the "merced" executable. 
    """

    return locateFileInBin( 'merced', exceptionIfNotFound )

