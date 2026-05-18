# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import unittest
from xml.etree import ElementTree as parser

from xData import table as tableModule

class TestTable(unittest.TestCase):

    def setUp(self):
        self.tt = tableModule.Table( columns=[
            tableModule.ColumnHeader( 0, 'energy', 'eV' ),
            tableModule.ColumnHeader( 1, 'spin', '' ),
            tableModule.ColumnHeader( 2, 'neutronWidth', 'eV' ),
            tableModule.ColumnHeader( 3, 'captureWidth', 'eV' ), ] )
        for dat in ([-1.7, tableModule.Blank(), 3.2, 0.089], [3.4, 0.5, 4.2, 0.072], [5.6, 1.5, 2.76, 0.064]):
            self.tt.addRow( dat )
        self.xmlstring = '\n'.join( self.tt.toXML_strList() )
        element = parser.fromstring( self.xmlstring )
        self.tt2 = tableModule.Table.parseNodeUsingClass(element, xPath = [],
                linkData = {'conversionTable' : {'index':int} })
        self.xmlstring2 = '\n'.join( self.tt2.toXML_strList() )


    def test_table(self):
        #print(self.xmlstring)
        self.assertEqual(self.xmlstring, self.xmlstring2)

    def test_getColumn(self):
        # extract a column, converting eV -> MeV
        self.assertListEqual(self.tt2.getColumn('captureWidth', 'MeV'), [8.899999999999999e-08, 7.2e-08, 6.4e-08])
        # THIS^ WILL NOT WORK WITH Blank() OBJECT TYPES, ONLY TESTS IF COLUMN IS EXTRACTED AND UNITS ARE CONVERTED

    def test_addrow(self):
        """ Test if the ValueError raises when adding row of unequal length"""
        with self.assertRaises(ValueError):
            self.tt.addRow([1,1,1])


# to run unit test
if __name__=="__main__":
    unittest.main()
