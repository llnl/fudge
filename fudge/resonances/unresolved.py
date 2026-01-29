# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Containers for unresolved resonance parameters
"""

import fractions
import abc

from LUPY import ancestry as ancestryModule

from PoPs import database as PoPsDatabaseModule
from pqu import PQU as PQUModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule
from xData import constant as constantModule
from xData.Documentation import documentation as documentationModule

from fudge import suites as suitesModule, abstractClasses as abstractClassesModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from .common import EnergyIntervals, ResonanceReactions, floatOrint, getAttrs
from .scatteringRadius import ScatteringRadius, HardSphereRadius


class Unresolved(abstractClassesModule.Component):
    moniker = 'unresolved'
    _requiredAttributes = (
        ('domainMin', float),
        ('domainMax', float),
        ('domainUnit', str)
    )

    def __init__(self, domainMin, domainMax, domainUnit):
        abstractClassesModule.Component.__init__(self, allowedClasses=(TabulatedWidths, EnergyIntervals))
        self.__domainMin = domainMin
        self.__domainMax = domainMax
        self.__domainUnit = domainUnit

    @property
    def domainMin(self):

        return self.__domainMin

    @property
    def domainMax(self):

        return self.__domainMax

    @property
    def domainUnit(self):

        return self.__domainUnit

    def toString(self, simpleString=False):
        return "Unresolved resonances in %s form\n" % self.evaluated.moniker

    def check(self, info):
        """
        Check that widths and level spacings span full URR domain, and that the energy grid is 'dense enough'.
        TODO: level spacing upper limit should decrease with increasing target mass
        TODO: reaction-specific limits on widths
        """
        from fudge import warning
        warnings = []
        for L in self.evaluated.Ls:
            for J in L.Js:
                if J.levelSpacing.data.domainMin > self.__domainMin or J.levelSpacing.data.domainMax < self.__domainMax:
                    warnings.append(warning.URRdomainMismatch(L.L, J.J, 'level spacing', J.levelSpacing))
                if J.levelSpacing.data.rangeMin <= 0 or J.levelSpacing.data.rangeMax > 5e+5:
                    warnings.append(warning.URRunphysicalLevelSpacing(L.L, J.J, J.levelSpacing))

                for width in J.widths:
                    if width.data.domainMin > self.__domainMin or width.data.domainMax < self.__domainMax:
                        warnings.append(warning.URRdomainMismatch(L.L, J.J, f"{width.resonanceReaction} width", width))
                    if width.data.rangeMin <= 0:
                        if width.data.rangeMin == 0 and width.resonanceReaction == 'competitive':
                            continue  # competitive may include threshold reactions
                        warnings.append(warning.URRunphysicalWidth(L.L, J.J, width.resonanceReaction, width))
                    if width.data.rangeMax > 1e+4:
                        warnings.append(warning.URRunphysicalWidth(L.L, J.J, width.resonanceReaction, width))

                elist = J.levelSpacing.data.convertAxisToUnit(1, self.__domainUnit).domainGrid
                missingPoints = [i1 for i1 in range(1, len(elist)) if elist[i1] > 3 * elist[i1 - 1]]
                for idx in missingPoints:
                    warnings.append(warning.URRinsufficientEnergyGrid(
                        L.L, J.J, PQUModule.PQU(elist[idx - 1], self.__domainUnit),
                        PQUModule.PQU(elist[idx], self.__domainUnit), J))
        return warnings

    def convertUnits(self, unitMap):

        if self.__domainUnit in unitMap:
            newUnit = unitMap[self.__domainUnit]
            factor = PQUModule.PQU(1, self.__domainUnit).getValueAs(newUnit)
            self.__domainMin *= factor
            self.__domainMax *= factor
            self.__domainUnit = newUnit
        for form in self:
            form.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s domainMin="%s" domainMax="%s" domainUnit="%s">' %
                     (indent, self.moniker,
                      PQUModule.floatToShortestString(self.__domainMin, 12),
                      PQUModule.floatToShortestString(self.__domainMax, 12),
                      self.__domainUnit)]
        xmlString += self.evaluated.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker

        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)

        attrs = getAttrs(element, required=Unresolved._requiredAttributes)
        URR = Unresolved(**attrs)
        for child in element:
            formClass = {
                TabulatedWidths.moniker: TabulatedWidths,
                EnergyIntervals.moniker: EnergyIntervals,
            }.get(child.tag)
            if formClass is None: raise Exception("unknown Unresolved resonance form '%s'!" % child.tag)
            URR.add(formClass.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop()
        return URR


class TabulatedWidths(ancestryModule.AncestryIO):
    moniker = 'tabulatedWidths'
    ancestryMembers = ('PoPs', 'scatteringRadius', 'hardSphereRadius', 'resonanceReactions', 'Ls')
    keyName = 'label'

    _requiredAttributes = (
        ('label', str),
        ('approximation', str),
    )
    _optionalAttributes = (
        ('calculateChannelRadius', bool, False),
        ('useForSelfShieldingOnly', bool, False)
    )

    def __init__(self, label, approximation, resonanceReactions, Ls, scatteringRadius=None, hardSphereRadius=None,
                 PoPs=None, calculateChannelRadius=False, useForSelfShieldingOnly=False):
        """
        Contains unresolved resonance parameters

        :param label: unique label corresponding to a style (i.e. 'eval')
        :param approximation: currently only 'SingleLevelBreitWigner' is supported
        :param resonanceReactions: ResonanceReactions instance
        :param Ls: URR_Lsections instance
        :param scatteringRadius: optional scatteringRadius instance
        :param hardSphereRadius: optional hardSphereRadius instance (used for phase shift only)
        :param PoPs: optional PoPs Database
        :param calculateChannelRadius: boolean
        :param useForSelfShieldingOnly: boolean
        """

        ancestryModule.AncestryIO.__init__(self)

        self.label = label
        self.approximation = approximation
        self.resonanceReactions = resonanceReactions
        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        self.PoPs = PoPs
        self.Ls = Ls
        self.calculateChannelRadius = calculateChannelRadius
        self.useForSelfShieldingOnly = useForSelfShieldingOnly

        self.__documentation = documentationModule.Documentation()
        self.__documentation.setAncestor(self)

    def __len__(self):
        return len(self.Ls)

    @property
    def documentation(self):
        """Returns the documentation instance."""

        return (self.__documentation)

    @property
    def resonanceReactions(self):
        return self.__resonanceReactions

    @resonanceReactions.setter
    def resonanceReactions(self, value):
        if not isinstance(value, ResonanceReactions):
            raise TypeError("Must be a ResonanceReactions instance")
        value.setAncestor(self)
        self.__resonanceReactions = value

    @property
    def scatteringRadius(self):

        return self.__scatteringRadius

    def getScatteringRadius(self):
        """Return ScatteringRadius, looking up ancestry if necessary"""
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            from .resonances import Resonances
            return self.findClassInAncestry(Resonances).getScatteringRadius()

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """ Can be set to None or to a scatteringRadius instance. """
        if value is not None:
            if not isinstance(value, ScatteringRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__scatteringRadius = value

    @property
    def hardSphereRadius(self):
        return self.__hardSphereRadius

    def getHardSphereRadius(self):
        """
        Return HardSphereRadius, looking up ancestry if necessary.
        If no HardSphereRadius is defined, return ScatteringRadius instead.
        """
        if self.__hardSphereRadius is not None:
            return self.__hardSphereRadius
        else:
            from .resonances import Resonances
            return self.findClassInAncestry(Resonances).getHardSphereRadius()

    @hardSphereRadius.setter
    def hardSphereRadius(self, value):
        """ Can be set to None or to a scatteringRadius instance. """
        if value is not None:
            if not isinstance(value, HardSphereRadius):
                raise TypeError("Hard sphere radius can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__hardSphereRadius = value

    @property
    def Ls(self):
        return self.__Ls

    @Ls.setter
    def Ls(self, value):
        if not isinstance(value, Lsections):
            raise TypeError("Must be a URR_Lsections instance")
        value.setAncestor(self)
        self.__Ls = value

    @property
    def PoPs(self):
        if self.__PoPs is not None:
            return self.__PoPs
        else:
            return self.rootAncestor.PoPs

    @PoPs.setter
    def PoPs(self, value):
        """ Can be set to None or to a PoPs Database instance """

        self.__PoPs = value
        if value is not None:
            if not isinstance(value, PoPsDatabaseModule.Database):
                raise TypeError("PoPs can't be set to type '%s'" % type(value))
            self.__PoPs.setAncestor(self)

    def getLocalPoPs(self):
        """ Returns None unless a local PoPs database is defined. """
        return self.__PoPs

    def convertUnits(self, unitMap):
        for child in self.__PoPs, self.__scatteringRadius, self.__hardSphereRadius:
            if child is not None:
                child.convertUnits(unitMap)
        for resonanceReaction in self.resonanceReactions:
            resonanceReaction.convertUnits(unitMap)
        for lval in self.Ls:
            for jval in lval.Js:
                jval.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attrs = ' label="%s" approximation="%s"' % (self.label, self.approximation)
        if kwargs.get('formatVersion', GNDS_formatVersionModule.default) >= GNDS_formatVersionModule.version_2_2:
            # always write it out
            attrs += ' calculateChannelRadius="%s"' % str(self.calculateChannelRadius).lower()
        if self.useForSelfShieldingOnly: attrs += ' useForSelfShieldingOnly="true"'

        xml = ['%s<%s%s>' % (indent, self.moniker, attrs)]

        xml += self.__documentation.toXML_strList(indent=indent2, **kwargs)

        if self.__PoPs:
            xml += self.__PoPs.toXML_strList(indent2, **kwargs)

        xml += self.resonanceReactions.toXML_strList(indent2, **kwargs)
        if self.__scatteringRadius:
            xml += self.__scatteringRadius.toXML_strList(indent2, **kwargs)
        if self.__hardSphereRadius:
            if self.__hardSphereRadius.isEnergyDependent() and kwargs.get(
                    'formatVersion', GNDS_formatVersionModule.default) < GNDS_formatVersionModule.version_2_2:
                # pretend it's a scattering radius for backwards compatibility
                scatteringRadius = ScatteringRadius(self.__hardSphereRadius.evaluated)
                xml += scatteringRadius.toXML_strList(indent2, **kwargs)
            else:
                xml += self.__hardSphereRadius.toXML_strList(indent2, **kwargs)
        xml += self.Ls.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker

        return xml

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        resonanceReactions = ResonanceReactions()
        child = node.find(ResonanceReactions.moniker)
        if child is not None: resonanceReactions.parseNode(child, xPath, linkData, **kwargs)

        radius = node.find(ScatteringRadius.moniker)
        if radius is not None:
            radius = ScatteringRadius.parseNodeUsingClass(radius, xPath, linkData, **kwargs)
        hsRadius = node.find(HardSphereRadius.moniker)
        if hsRadius is not None:
            hsRadius = HardSphereRadius.parseNodeUsingClass(hsRadius, xPath, linkData, **kwargs)
        pops = node.find(PoPsDatabaseModule.Database.moniker)
        if pops is not None:
            pops = PoPsDatabaseModule.Database.parseNodeUsingClass(pops, xPath, linkData, **kwargs)
        Ls = Lsections()
        Ls.parseNode(node.find(Ls.moniker), xPath, linkData, **kwargs)

        result = cls(node.get('label'), node.get('approximation'),
                     resonanceReactions, Ls=Ls, scatteringRadius=radius, hardSphereRadius=hsRadius, PoPs=pops,
                     calculateChannelRadius=node.get('calculateChannelRadius') == 'true',
                     useForSelfShieldingOnly=node.get('useForSelfShieldingOnly') == 'true'
                     )
        if kwargs.get('formatVersion', GNDS_formatVersionModule.default) < GNDS_formatVersionModule.version_2_2:
            # calculateChannelRadius wasn't stored prior to GNDS-2.2, but is usually True. Exception: Re185
            result.calculateChannelRadius = True

        documentation = node.find(documentationModule.Documentation.moniker)
        if documentation is not None:
            result.documentation.parseNode(documentation, xPath, linkData, **kwargs)

        xPath.pop()
        return result


class InterpolationTable(ancestryModule.AncestryIO, metaclass=abc.ABCMeta):
    """
    Base class for unresolved level spacing and widths
    Data are stored as XYs1d or Regions1d.
    """

    allowedSubclasses = (XYs1dModule.XYs1d, regionsModule.Regions1d, constantModule.Constant1d)

    def __init__(self, data):

        ancestryModule.AncestryIO.__init__(self)
        self.data = data

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, value):
        if not isinstance(value, self.allowedSubclasses):
            raise TypeError("levelSpacing must be an XYs1d, Regions1d or Constant1d instance")
        value.setAncestor(self)
        self.__data = value

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        child, = element[:]
        data = None
        for subclass in cls.allowedSubclasses:
            if child.tag == subclass.moniker:
                data = subclass.parseNodeUsingClass(child, xPath, linkData, **kwargs)
                break
        result = cls(data)
        xPath.pop()
        return result


class Lsections(suitesModule.Suite):
    """
    tabulated widths contains a list of 'L' sections, each of which contains a list of 'J' sections
    """

    moniker = "Ls"

    def __init__(self):
        suitesModule.Suite.__init__(self, (Lsection,))


class Lsection(ancestryModule.AncestryIO):
    """ unresolved average widths, grouped by L. Contains list of J-values: """

    moniker = 'L'
    keyName = 'label'

    def __init__(self, label, L, Js):
        ancestryModule.AncestryIO.__init__(self)
        self.label = label
        self.L = L
        self.Js = Js

    @property
    def value(self): return self.L

    @property
    def Js(self):
        return self.__Js

    @Js.setter
    def Js(self, value):
        if not isinstance(value, (Jsections)):
            raise TypeError("Must be URR_Jsections instance")
        value.setAncestor(self)
        self.__Js = value

    def toXML_strList(self, indent='', **kwargs):
        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s label="%s" value="%d">' % (indent, self.moniker, self.label, self.value)]
        xml += self.Js.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append('%s[@label="%s"]' % (element.tag, element.get('label')))
        Js = Jsections()
        Js.parseNode(element.find(Js.moniker), xPath, linkData, **kwargs)
        result = cls(element.get('label'), int(element.get('value')), Js)
        xPath.pop()
        return result


class Jsections(suitesModule.Suite):
    moniker = "Js"

    def __init__(self):
        suitesModule.Suite.__init__(self, (Jsection,))


class Jsection(ancestryModule.AncestryIO):
    """
    Unresolved average parameters for a specific L/J.
    Contains interpolation tables for the level spacing and widths, plus degrees of freedom for each open
    channel  (degrees of freedom are typically 1 or 2 indicating how many channel spins contribute)
    """

    moniker = 'J'
    keyName = 'label'

    def __init__(self, label, J, levelSpacing, widths):

        ancestryModule.AncestryIO.__init__(self)
        self.label = label
        self.J = J
        self.levelSpacing = levelSpacing
        self.widths = widths

    @property
    def value(self):
        return str(self.J)

    @property
    def levelSpacing(self):
        return self.__levelSpacing

    @levelSpacing.setter
    def levelSpacing(self, spacing):
        if not isinstance(spacing, LevelSpacing):
            raise TypeError("Expected levelSpacing instance, got %s instead" % type(spacing))
        spacing.setAncestor(self)
        self.__levelSpacing = spacing

    @property
    def widths(self):
        return self.__widths

    @widths.setter
    def widths(self, widths):
        if not isinstance(widths, Widths):
            raise TypeError("Expected URR_widths instance, got %s instead" % type(widths))
        widths.setAncestor(self)
        self.__widths = widths

    def convertUnits(self, unitMap):

        self.levelSpacing.data.convertUnits(unitMap)
        for width in self.widths:
            width.data.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s label="%s" value="%s">' % (indent, self.moniker, self.label, self.J)]
        xml += self.levelSpacing.toXML_strList(indent2, **kwargs)
        xml += self.widths.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append('%s[@label="%s"]' % (element.tag, element.get('label')))
        levelSpacing = LevelSpacing.parseNodeUsingClass(element.find(LevelSpacing.moniker), xPath, linkData, **kwargs)
        widths = Widths()
        widths.parseNode(element.find(Widths.moniker), xPath, linkData, **kwargs)
        Jsec = Jsection(element.get('label'), fractions.Fraction(element.get('value')), levelSpacing, widths)

        xPath.pop()
        return Jsec


class LevelSpacing(InterpolationTable):
    """
    Contains the average level spacing, stored in XYs1d or Regions1d.
    """

    moniker = "levelSpacing"


class Widths(suitesModule.Suite):
    """
    Contains all average channel widths for an L/J combination.
    """

    moniker = 'widths'

    def __init__(self):
        suitesModule.Suite.__init__(self, allowedClasses=(Width,))


class Width(InterpolationTable):
    """
    Stores the average width and number of degrees of freedom for a single channel

    Degrees of freedom are typically 0, 1 or 2 depending on how many channel spins contribute, but it may be
    non-integer for a reaction that has a threshold somewhere inside the unresolved region.

    Average width is stored as XYs1d or Regions1d.
    """

    moniker = 'width'
    keyName = 'label'

    def __init__(self, label, resonanceReaction, degreesOfFreedom, data):
        InterpolationTable.__init__(self, data)
        self.label = label
        self.resonanceReaction = resonanceReaction
        self.degreesOfFreedom = degreesOfFreedom

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s label="%s" resonanceReaction="%s" degreesOfFreedom="%g">' %
               (indent, self.moniker, self.label, self.resonanceReaction, self.degreesOfFreedom)]
        xml += self.data.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker

        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)
        child, = element[:]
        data = None
        for subclass in cls.allowedSubclasses:
            if child.tag == subclass.moniker:
                data = subclass.parseNodeUsingClass(child, xPath, linkData, **kwargs)
                break
        result = cls(element.get('label'), element.get('resonanceReaction'),
                     floatOrint(element.get('degreesOfFreedom', 0)), data)
        xPath.pop()
        return result
