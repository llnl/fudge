# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from math import *
C = 2.0
n = 10
Emax = 4.0
mumax = 1 - C*C/(Emax*Emax)
dmu = (1 + mumax)/n
for i in range(n+1):
    mu = -1 + i*dmu
    E = C/sqrt(1 - mu)
    print('(', mu, ',', E, ')%')
C = 3.0
mumax = 1 - C*C/(Emax*Emax)
dmu = (1 + mumax)/n
for i in range(n+1):
    mu = -1 + i*dmu
    E = C/sqrt(1 - mu)
    print('(', mu, ',', E, ')%')
