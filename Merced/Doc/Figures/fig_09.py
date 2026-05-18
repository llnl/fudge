# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# python to calculate the E_out energy bins in terms of E_in and mu
E0 = 2.5   # value of E_in at mu = 1
E_out = 0.2
m_e = 0.511  # electron rest mass
mu_blow = 1 - m_e/E_out  # where we divide by zero
if mu_blow < -1:
    mu_min = -1.0
else:
    mu_min = mu_blow
n = 15
dmu = (1 - mu_min)/n
for i in range(n):
    mu = mu_min + (i+1)*dmu
    E_in = E0/( 1 - (1-mu)*E_out/m_e )
    print('(', mu, ',', E_in, ')%')

E0 = 2.1   # value of E_in at mu = 1
E_out = 0.2
m_e = 0.511  # electron rest mass
mu_blow = 1 - m_e/E_out  # where we divide by zero
if mu_blow < -1:
    mu_min = -1.0
else:
    mu_min = mu_blow
n = 15
dmu = (1 - mu_min)/n
for i in range(n):
    mu = mu_min + (i+1)*dmu
    E_in = E0/( 1 - (1-mu)*E_out/m_e )
    print('(', mu, ',', E_in, ')%')
