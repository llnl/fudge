# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Defines the supported version numbers for the GNDS format (not the same as the FUDGE version number)."""

version_1_10 = '1.10'
version_2_0_LLNL_3 = '2.0.LLNL_3'
version_2_0_LLNL_4 = '2.0.LLNL_4'
version_2_0        = '2.0'
version_2_1        = '2.1'
version_2_2        = '2.2.rc1'

allowed = (version_1_10, version_2_0, version_2_1, version_2_2)
allowedPlus = allowed + (version_2_0_LLNL_3, version_2_0_LLNL_4)

default = version_2_0       # version written by default in translation scripts

MAJORVERSION, MINORVERSION = default.split( '.' )[:2]

def versionComponents(version):
    '''
    Returns the major, minor and beta components of the GNDS format version. Both major and minor are returned as integers, and
    beta is returns as a string. If beta is not present in *version*, it is returned as an empty string. For example, if version is 
    "1.10" then (1, 10, '') is returned; while, if version is "2.0.LLNL_4" then (2, 0, 'LLNL_4') is returned.

    :param version:         The string representing for version.

    :return:                The components of version as the tuple (major, minor, beta).
    '''

    if version.count('.') == 1:
        major, minor = version.split('.')
        beta = ''
    else:
        major, minor, beta = version.split('.')

    return int(major), int(minor), beta

def constructVersionFromComponents(major, minor, beta=''):
    '''
    Returns a verion string for the components major, minor and beta. Major and minor can be integers or strings. Beta must be a
    string. If beta is an empty string, then the return string only contains the major and minor components. For example, if 
    major, minor and beta are 2, 0 and '', respectively, then '2.0' is returned; while, if major, minor and beta are 2, 0 and
    'LLNL_4' , respectively, then '2.0.LLNL_4' is returned.

    :param major:       The major component of the version.
    :param minor:       The minor component of the version.
    :param beta:        The optional beta component of the version.

    :return:            String representation of the version components.
    '''

    version = '%s.%s' % (major, minor)
    if beta != '': version += '.%s' % beta

    return version
