# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the FissionFragmentData class, used to store spontaneous fission product yields
and the delayed neutrons and gammas emitted by the decay of those fission products.

 CMM documentation comments:
  - for common stuff like convertUnits, define them once and refer to them from other classes?
  - similar question on toXML_strList
"""

from LUPY import ancestry as ancestryModule

from . import delayedNeutronData as delayedNeutronDataModule
from . import productYield as productYieldModule

class FissionFragmentData(ancestryModule.AncestryIO):
    """This class stores information about spontaneous fission fragments and their decays."""

    moniker = 'fissionFragmentData'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__(self)

        # FIXME doesn't conform with GNDS specs. Should be Suite of delayedNeutron.
        self.__delayedNeutronData = delayedNeutronDataModule.DelayedNeutronData( )
        self.__delayedNeutronData.setAncestor( self )

        # FIXME should also support fissionEnergyRelease section, even if no files currently use it.

        self.__productYields = productYieldModule.Suite( )
        self.__productYields.setAncestor( self )
        
    @property
    def delayedNeutronData( self ) :
        """This accesses the delayed neutron data within a **FissionFragmentData** instance."""

        return( self.__delayedNeutronData )

    @property
    def productYields( self ) :
        """This accesses the fission product yield data within a **FissionFragmentData** instance."""


        return( self.__productYields )

    def convertUnits( self, unitMap ) :
        """
        The convertUnits method recursively searches the **FissionFragmentData** and child classes
        for each unit in the **unitMap** dictionary and converts to the desired new unit.

        :param unitMap: dict storing {'oldUnit': 'newUnit'}. Units must be compatible, e.g. {'eV': 'MeV'} but not {'amu': 'MeV'}.
        """

        self.__delayedNeutronData.convertUnits( unitMap )
        self.__productYields.convertUnits( unitMap )

    def replicate( self, other ) :
        # FIXME why replicate instead of copy? Does this method get used? Recommend removing it if not.

        self.__delayedNeutronData = other.delayedNeutronData    # FIXME not making a copy?
        self.__delayedNeutronData.setAncestor( self )

        self.__productYields = other.productYields
        self.__productYields.setAncestor( self )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:    The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:    A keyword list.

        :return:          List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLStringList = [ '%s<%s> ' % (indent, self.moniker) ]
        XMLStringList += self.__delayedNeutronData.toXML_strList(indent2, **kwargs)
        XMLStringList += self.__productYields.toXML_strList(indent2, **kwargs)
        XMLStringList[-1] += '</%s>' % self.moniker

        if len(XMLStringList) == 1: XMLStringList = []

        return XMLStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Parse node from a file, returning a new **FissionFragmentData** instance unless errors were encountered during parsing.

        :param cls: FUDGE Python class to return.
        :param node: **FissionFragmentData** node to parse.
        :param xPath: List containing xPath to current node, useful mostly for debugging.
        :param linkData: dict that collects unresolved links.
        """

        xPath.append( node.tag )
        FFD = cls()
        for child in node:
            if child.tag == delayedNeutronDataModule.DelayedNeutronData.moniker:
                FFD.delayedNeutronData.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == productYieldModule.Suite.moniker:
                FFD.productYields.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Encountered unknown node '%s' in %s" % (child.tag, node.tag))

        xPath.pop()
        return FFD
