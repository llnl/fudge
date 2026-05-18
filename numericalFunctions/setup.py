# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
# This file is required for the C extensions
# The main configuration is in pyproject.toml, but setuptools still uses this for extensions

import glob, os
from setuptools import setup, Extension

def setup_package():
    statusMessageReporting_c = glob.glob(os.path.join('statusMessageReporting', 'Src', '*.c'))

    ptwC_c = glob.glob(os.path.join('ptwC', 'Src', '*.c'))

    ptwX_c = glob.glob(os.path.join('ptwX', 'Src', '*.c'))
    ptwX_Py_c = glob.glob(os.path.join('ptwX', 'Python', 'Src', '*.c'))

    ptwXY_c = glob.glob(os.path.join('ptwXY', 'Src', '*.c'))
    ptwXY_Py_c = glob.glob(os.path.join('ptwXY', 'Python', 'Src', '*.c'))

    nf_Legendre_c = glob.glob(os.path.join('nf_Legendre', 'Src', '*.c'))
    nf_Legendre_Py_c = glob.glob(os.path.join('nf_Legendre', 'Python', 'Src', '*.c'))

    nf_specialFunctions_c = glob.glob(os.path.join('nf_specialFunctions', 'Src', 'nf_[egip]*.c'))
    nf_specialFunctions_Py_c = glob.glob(
        os.path.join('nf_specialFunctions', 'Python', 'Src', 'nf_specialFunctions_C.c'))

    nf_angularMomentumCoupling_c = glob.glob(
        os.path.join('nf_specialFunctions', 'Src', 'nf_angularMomentumCoupling.c'))
    nf_angularMomentumCoupling_Py_c = glob.glob(
        os.path.join('nf_specialFunctions', 'Python', 'Src', 'nf_angularMomentumCoupling_C.c'))

    nf_integration_c = glob.glob(os.path.join('nf_integration', 'Src', '*.c'))
    nf_integration_Py_c = glob.glob(os.path.join('nf_integration', 'Python', 'Src', '*.c'))

    statusMessageReporting_hDir = os.path.join('statusMessageReporting', 'Src')

    ptwC_hDir = os.path.join('ptwC', 'Src')
    ptwX_hDir = os.path.join('ptwX', 'Src')
    ptwXY_hDir = os.path.join('ptwXY', 'Src')
    ptwXY_Py_hDir = os.path.join('ptwXY', 'Python', 'Src')
    nf_Legendre_hDir = os.path.join('nf_Legendre', 'Src')
    nf_specialFunctions_hDir = os.path.join('nf_specialFunctions', 'Src')
    nf_specialFunctions_Py_hDir = os.path.join('nf_specialFunctions', 'Python', 'Src')
    nf_integration_hDir = os.path.join('nf_integration', 'Src')

    extra_compile_args = []

    listOfDoubles_C = Extension('numericalFunctions.listOfDoubles_C',
        extra_compile_args = extra_compile_args,
        sources = statusMessageReporting_c + ptwC_c + ptwX_c + ptwX_Py_c,
        include_dirs = [statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir])
    
    pointwiseXY_C = Extension('numericalFunctions.pointwiseXY_C',
        extra_compile_args = extra_compile_args,
        sources = statusMessageReporting_c + ptwC_c + ptwX_c + nf_Legendre_c + nf_integration_c + ptwXY_c + ptwXY_Py_c,
        include_dirs = [statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir, nf_Legendre_hDir, nf_integration_hDir, 
                ptwXY_hDir, ptwXY_Py_hDir])

    Legendre = Extension('numericalFunctions.Legendre',
        extra_compile_args = extra_compile_args,
        sources = statusMessageReporting_c + ptwC_c + ptwX_c + ptwXY_c + nf_integration_c + nf_Legendre_c + nf_Legendre_Py_c,
        include_dirs = [statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir, ptwXY_hDir, ptwXY_Py_hDir, 
                nf_integration_hDir, nf_Legendre_hDir])
    
    specialFunctions = Extension('numericalFunctions.specialFunctions',
        extra_compile_args = extra_compile_args,
        sources = statusMessageReporting_c + ptwC_c + nf_specialFunctions_c + nf_specialFunctions_Py_c,
        include_dirs = [statusMessageReporting_hDir, ptwC_hDir, nf_specialFunctions_hDir, nf_specialFunctions_Py_hDir])
    
    angularMomentumCoupling = Extension('numericalFunctions.angularMomentumCoupling',
        extra_compile_args = extra_compile_args,
        sources = statusMessageReporting_c + ptwC_c + nf_angularMomentumCoupling_c + nf_angularMomentumCoupling_Py_c,
        include_dirs = [statusMessageReporting_hDir, ptwC_hDir, nf_specialFunctions_hDir, nf_specialFunctions_Py_hDir])
    
    integration = Extension('numericalFunctions.integration',
        extra_compile_args = extra_compile_args,
        sources = statusMessageReporting_c + ptwC_c + ptwX_c + ptwXY_c + nf_Legendre_c + nf_integration_c + nf_integration_Py_c,
        include_dirs = [statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir, ptwXY_hDir, nf_Legendre_hDir, nf_integration_hDir])

    # This setup() is called by setuptools from pyproject.toml
    setup(
        ext_modules=[listOfDoubles_C, pointwiseXY_C, Legendre, specialFunctions, angularMomentumCoupling, integration],
    )

if __name__ == '__main__':
    setup_package()