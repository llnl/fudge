# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# the curve in Fig. 6

x0 = 0.3
x1 = 0.95
y0 = 0.8
y1 = 0.95
m = (y1 - y0)/(x1 - x0)
n = 8
dx = (x1 - x0)/n
a = 0.3
for i in range(n+1):
    x = x0 + i*dx
    y = y0 + m*(x - x0) + a*(x - x0)*(x1 - x)
    print('(', x, ', ', y, ')%')
