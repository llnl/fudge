# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Container for all resonance parameters (resolved, unresolved and/or scattering radius)
"""

from pqu import PQU as PQUModule
from LUPY import ancestry as ancestryModule


class Resonances(ancestryModule.AncestryIO):
    """ 
    This is the top-level class for storing resonance parameters.
    For light targets it may contain only a scattering radius. For heavier targets it typically
    contains a resolved and/or unresolved section.
    
    resonances also has a boolean flag 'reconstructCrossSection'. If False, either the cross section
    has already been reconstructed, or the parameters are given for information only and no reconstruction
    should be performed.
    """

    moniker = 'resonances'
    ancestryMembers = ('scatteringRadius', 'hardSphereRadius', 'resolved', 'unresolved')

    def __init__(self, scatteringRadius, hardSphereRadius=None, resolved=None, unresolved=None):

        ancestryModule.AncestryIO.__init__( self )

        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        self.resolved = resolved
        self.unresolved = unresolved

    def __str__(self):
        """ string representation """
        return self.toString(simpleString=False)

    @property
    def scatteringRadius(self):
        """Returns a reference to self's scatteringRadius."""

        return self.__scatteringRadius

    def getScatteringRadius(self):
        """Returns the scattering radius. Provided for consistency with getHardSphereRadius"""

        return self.scatteringRadius

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        from . import scatteringRadius as scatteringRadiusModule

        if not isinstance(value, (scatteringRadiusModule.ScatteringRadius,)):
            raise TypeError("Expected ScatteringRadius, got %s" % type(value))
        value.setAncestor(self)
        self.__scatteringRadius = value

    @property
    def hardSphereRadius(self):
        """Returns a reference to self's hardSphereRadius."""

        return self.__hardSphereRadius

    def getHardSphereRadius(self):
        """Returns hardSphereRadius, or reverts to scatteringRadius if no hard sphere radius is defined"""

        if self.__hardSphereRadius is not None:
            return self.__hardSphereRadius
        return self.__scatteringRadius

    @hardSphereRadius.setter
    def hardSphereRadius(self, value):
        """Can be set to None or to a HardSphereRadius instance."""

        if value is not None:
            from . import scatteringRadius as scatteringRadiusModule
            if not isinstance(value, scatteringRadiusModule.HardSphereRadius):
                raise TypeError("Expected HardSphereRadius, got %s" % type(value))
            value.setAncestor(self)
        self.__hardSphereRadius = value

    @property
    def resolved(self):
        """Returns a reference to self's resolved."""

        return self.__resolved

    @resolved.setter
    def resolved(self, value):
        """Can be set to None or to a Resolved instance."""

        if value is not None:
            from . import resolved as resolvedModule
            if not isinstance(value, resolvedModule.Resolved):
                raise TypeError("Expected Resolved, got %s" % type(value))
            value.setAncestor(self)
        self.__resolved = value

    @property
    def unresolved(self):
        """Returns a reference to self's unresolved."""

        return self.__unresolved

    @unresolved.setter
    def unresolved(self, value):
        """Can be set to None or to an Unresolved instance."""

        if value is not None:
            from . import unresolved as unresolvedModule
            if not isinstance(value, unresolvedModule.Unresolved):
                raise TypeError("Expected Unresolved, got %s" % type(value))
            value.setAncestor(self)
        self.__unresolved = value

    def convertUnits(self, unitMap):
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key
        (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        for child in Resonances.ancestryMembers:
            if getattr(self, child) is not None:
                getattr(self, child).convertUnits(unitMap)

    def check( self, info ):
        from fudge import warning
        warnings = []
        for child in Resonances.ancestryMembers:
            section = getattr(self, child)
            if section is not None:
                warningList = section.check(info)
                if warningList:
                    warnings.append(warning.Context(section.moniker, warningList))
        return warnings

    def containsStyle(self, styleLabel):
        # FIXME: should search scatteringRadius and hardSphereRadius too after converting them to suites
        for section in (self.resolved, self.unresolved):
            if section is not None and styleLabel in section:
                return True
        return False

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s>' % (indent, self.moniker)]
        for child in self.ancestryMembers:
            section = getattr(self, child)
            if section is not None:
                xmlString += section.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def domain(self, unitTo=None, asPQU=False):
        """ Return resonance region domain as a tuple of floats: (lowest edge, highest edge).

        options:
          unitTo: convert output to specified unit (given as a string).
          asPQU = True: return a tuple of PhysicalQuantityWithUncertainty instances instead of floats.
        """

        bounds = [(self.scatteringRadius.domainMin, self.scatteringRadius.domainMax)]
        if self.hardSphereRadius:
            bounds.append((self.hardSphereRadius.domainMin, self.hardSphereRadius.domainMax))
        if self.resolved:
            if self.resolved.multipleRegions:
                bounds += [(reg.domainMin, reg.domainMax) for reg in self.resolved.regions]
            else: bounds.append((self.resolved.domainMin, self.resolved.domainMax))
        if self.unresolved:
            bounds.append((self.unresolved.domainMin, self.unresolved.domainMax))

        for idx in range(len(bounds)-1):
            assert bounds[idx][1] == bounds[idx+1][0], "Resonance region boundaries don't match!"

        if asPQU:
            return (bounds[0][0], bounds[-1][1])
        elif unitTo:
            return (bounds[0][0].getValue(unitTo), bounds[-1][1].getValue(unitTo))
        else:
            return (bounds[0][0].value, bounds[-1][1].value)

    @property
    def reconstructCrossSection(self):
        if self.resolved is not None and not self.resolved.evaluated.useForSelfShieldingOnly: return True
        if self.unresolved is not None and not self.unresolved.evaluated.useForSelfShieldingOnly: return True
        return False

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        from .resolved import Resolved
        from .scatteringRadius import ScatteringRadius, HardSphereRadius
        from .unresolved import Unresolved

        xPath.append(element.tag)

        scatRadius, hsRadius, RRR, URR = None, None, None, None
        for child in element:
            if child.tag == ScatteringRadius.moniker:
                scatRadius = ScatteringRadius.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == HardSphereRadius.moniker:
                hsRadius = HardSphereRadius.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == Resolved.moniker:
                RRR = Resolved.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == Unresolved.moniker:
                URR = Unresolved.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise Exception("unknown element '%s' encountered in resonances!" % child.tag)

        res = cls(scatteringRadius=scatRadius, hardSphereRadius=hsRadius, resolved=RRR, unresolved=URR)
        xPath.pop()
        return res

    def toString(self, simpleString=False):
        """Returns a string representation of self. If simpleString is True, 
        the string contains only an overview without listing resonances"""
        s = 'Resonances:\n'
        if self.scatteringRadius:
            s += self.scatteringRadius.toString(simpleString=simpleString)
        if self.hardSphereRadius:
            s += self.hardSphereRadius.toString(simpleString=simpleString)
        if self.resolved:
            s += self.resolved.toString(simpleString=simpleString)
        if self.unresolved:
            s += self.unresolved.toString(simpleString=simpleString)
        return s
