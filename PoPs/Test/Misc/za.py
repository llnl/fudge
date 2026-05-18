# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs import database as databaseModule

from PoPs import misc as miscModule

pops = databaseModule.read( 'p.xml' )

ids_ZAs = []
for particle in pops :
    ids_ZAs.append( [ miscModule.ZA( particle ), particle.id ] )

ids_ZAs.sort( )
for ZA, _id in ids_ZAs : print("%-8s %s" % (ZA, _id))

ids_ZAs = []
for particle in pops :
    ids_ZAs.append( [ miscModule.ZAInfo( particle ), particle.id ] )

ids_ZAs.sort( )
for ZA, _id in ids_ZAs : print("%-8s %s" % (ZA, _id))
