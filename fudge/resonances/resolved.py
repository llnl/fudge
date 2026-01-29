# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Containers for resolved resonance parameters
"""

from pqu import PQU as PQUModule
from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule

from xData.Documentation import documentation as documentationModule

from PoPs import database as PoPsDatabaseModule

from fudge import abstractClasses as abstractClassesModule, suites as suitesModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge.resonances.resonances import Resonances
from fudge.resonances import externalRMatrix as externalRMatrixModule
from fudge.resonances.scatteringRadius import ScatteringRadius, HardSphereRadius
from fudge.resonances.common import Parity, Spin, ResonanceParameters, getAttrs, ResonanceReactions, EnergyIntervals

class Resolved(abstractClassesModule.Component):
    """ class for resolved resonances """

    moniker = 'resolved'
    _requiredAttributes = (
        ('domainMin', float),
        ('domainMax', float),
        ('domainUnit', str)
    )

    def __init__(self, domainMin, domainMax, domainUnit):

        abstractClassesModule.Component.__init__(self, allowedClasses=(EnergyIntervals, BreitWigner, RMatrix))
        self.__domainMin = domainMin
        self.__domainMax = domainMax
        self.__domainUnit = domainUnit

    @property
    def domainMin(self):

        return self.__domainMin

    @domainMin.setter
    def domainMin(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("domainMin must be numeric, got %s" % type(value))
        self.__domainMin = value

    @property
    def domainMax(self):

        return self.__domainMax

    @domainMax.setter
    def domainMax(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("domainMax must be numeric, got %s" % type(value))
        self.__domainMax = value

    @property
    def domainUnit(self):

        return self.__domainUnit

    def toString(self, simpleString=False):
        if isinstance(self.evaluated, EnergyIntervals):
            return "Resolved region with DEPRECATED multiple regions\n"
        else:
            return "Resolved resonances in %s form\n" % self.evaluated.moniker

    def check(self, info):
        from fudge import warning
        warnings = []
        if isinstance(self.evaluated, EnergyIntervals):
            warnings.append(warning.RRmultipleRegions())
        warningList = self.evaluated.check(info)
        if warningList:
            warnings.append(warning.Context(self.evaluated.moniker, warningList))
        return warnings

    def convertUnits(self, unitMap):

        if self.__domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            factor = PQUModule.PQU(1, self.domainUnit).getValueAs(newUnit)
            self.__domainMin *= factor
            self.__domainMax *= factor
            self.__domainUnit = newUnit
        for form in self:
            form.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s domainMin="%s" domainMax="%s" domainUnit="%s">' % (indent, self.moniker,
                     PQUModule.floatToShortestString(self.__domainMin, 12),
                     PQUModule.floatToShortestString(self.__domainMax, 12), self.__domainUnit)]
        for form in self: xmlString += form.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker

        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)

        RRR = cls(**getAttrs(element, required=cls._requiredAttributes))
        for child in element:
            formClass = {
                BreitWigner.moniker: BreitWigner,
                RMatrix.moniker: RMatrix,
                EnergyIntervals.moniker: EnergyIntervals,
            }.get(child.tag)
            if formClass is None: raise Exception("unknown resolved resonance form '%s'!" % child.tag)
            RRR.add(formClass.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop()
        return RRR


class BreitWigner(ancestryModule.AncestryIO):

    moniker = 'BreitWigner'
    ancestryMembers = ('PoPs', 'scatteringRadius', 'hardSphereRadius', 'resonanceParameters')
    keyName = 'label'

    _requiredAttributes = (
        ('label', str),
        ('approximation', str),
    )
    _optionalAttributes = (
        ('calculateChannelRadius', bool, False),
        ('useForSelfShieldingOnly', bool, False),
        ('computeAngularDistribution', bool, False),    # not supported by GNDS-2.0
    )

    class Approximation(enumsModule.Enum):
        """Defines the allowed approximations for a BreitWigner instance."""

        SingleLevel = enumsModule.auto()
        MultiLevel = enumsModule.auto()

        @classmethod
        def checkEnumOrString(cls, value):
            """Some GNDS-2.0 files were published with incorrect 'approximation' strings. Support reading those."""
            value = {
                'singleLevel': cls.SingleLevel,
                'multiLevel': cls.MultiLevel
            }.get(value, value)
            return super().checkEnumOrString(value)

    def __init__(self, label, approximation, resonanceParameters=None, scatteringRadius=None,
                 hardSphereRadius=None, PoPs=None, **kwargs):
        """
        Container for resonance parameters using Single-Level or Multi-Level Breit-Wigner approximation

        :param label: corresponds to a style (i.e. 'eval')
        :param approximation: 'SingleLevel' or 'MultiLevel'
        :param resonanceParameters: optional ResonanceParameters instance
        :param scatteringRadius: optional ScatteringRadius instance
        :param hardSphereRadius: optional HardSphereRadius instance
        :param PoPs: optional PoPs.Database instance
        :param kwargs: see _optionalAttributes
        """

        ancestryModule.AncestryIO.__init__(self)

        for attr, Type, default in self._optionalAttributes:
            val = kwargs.get(attr, default)
            assert type(val) is Type
            setattr(self, attr, val)

        self.label = label

        self.approximation = BreitWigner.Approximation.checkEnumOrString(approximation)

        self.resonanceParameters = resonanceParameters
        if self.resonanceParameters is not None:
            self.resonanceParameters.setAncestor(self)
        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        self.PoPs = PoPs

        self.__documentation = documentationModule.Documentation()
        self.__documentation.setAncestor(self)

    @property
    def documentation(self):
        """Returns the documentation instance."""

        return self.__documentation

    def check(self, info):
        return _resonance_checker(self, info, [self.scatteringRadius, self.resonanceParameters])

    def convertUnits(self, unitMap):

        for child in self.ancestryMembers:
            if getattr(self, child) is not None:
                getattr(self, child).convertUnits(unitMap)

    def __getitem__(self, idx):
        return self.resonanceParameters.table[idx]

    def __len__(self):
        return len(self.resonanceParameters.table)

    def addResonance(self, resonance):
        """ insert a new resonance in the resonance parameter table """
        # resonance = (energy, J, l, ... )
        self.resonanceParameters.table.addRow(resonance)

    @property
    def scatteringRadius(self):

        return self.__scatteringRadius

    def getScatteringRadius(self):
        """Return ScatteringRadius, looking up ancestry if necessary."""
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            return self.findClassInAncestry(Resonances).getScatteringRadius()

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """Can be set to None or to a scatteringRadius instance."""

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
            return self.findClassInAncestry(Resonances).getHardSphereRadius()

    @hardSphereRadius.setter
    def hardSphereRadius(self, value):
        """Can be set to None or to a hardSphereRadius instance."""

        if value is not None:
            if not isinstance(value, HardSphereRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__hardSphereRadius = value

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

    def toXML_strList(self, indent='', **kwargs):

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        xmlString = '%s<%s label="%s" approximation="%s"' % (indent, self.moniker, self.label, self.approximation)
        for attr, Type, default in self._optionalAttributes:
            if getattr(self, attr) != default:
                attrVal = getattr(self, attr)
                if Type is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (attr, attrVal)
        xmlString = [xmlString + '>']

        xmlString += self.__documentation.toXML_strList(indent=indent2, **kwargs)

        if self.__PoPs is not None:
            xmlString += self.__PoPs.toXML_strList(indent2, **kwargs)
        if self.__scatteringRadius is not None:
            xmlString += self.__scatteringRadius.toXML_strList(indent2, **kwargs)
        if self.__hardSphereRadius is not None:
            if self.__hardSphereRadius.isEnergyDependent() and kwargs.get(
                    'formatVersion', GNDS_formatVersionModule.default) < GNDS_formatVersionModule.version_2_2:
                # pretend it's a scattering radius for backwards compatibility
                scatteringRadius = ScatteringRadius(self.__hardSphereRadius.evaluated)
                xmlString += scatteringRadius.toXML_strList(indent2, **kwargs)
            else:
                xmlString += self.__hardSphereRadius.toXML_strList(indent2, **kwargs)
        if self.resonanceParameters:
            xmlString.extend(self.resonanceParameters.toXML_strList(indent2, **kwargs))
        xmlString[-1] += '</%s>' % self.moniker

        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        pops = element.find(PoPsDatabaseModule.Database.moniker)
        if pops is not None:
            pops = PoPsDatabaseModule.Database.parseNodeUsingClass(pops, xPath, linkData, **kwargs)

        radius = element.find(ScatteringRadius.moniker)
        if radius is not None:
            radius = ScatteringRadius.parseNodeUsingClass(radius, xPath, linkData, **kwargs)
        hsradius = element.find(HardSphereRadius.moniker)
        if hsradius is not None:
            hsradius = HardSphereRadius.parseNodeUsingClass(hsradius, xPath, linkData, **kwargs)
        parameters = ResonanceParameters.parseNodeUsingClass(
                element.find(ResonanceParameters.moniker), xPath, linkData, **kwargs)
        attrs = getAttrs(element, required=cls._requiredAttributes, optional=cls._optionalAttributes)
        label = attrs.pop("label")
        approximation = attrs.pop("approximation")
        resonanceData = cls(label, approximation, parameters, scatteringRadius=radius,
                            hardSphereRadius=hsradius, PoPs=pops, **attrs)

        documentation = element.find(documentationModule.Documentation.moniker)
        if documentation is not None:
            resonanceData.documentation.parseNode(documentation, xPath, linkData, **kwargs)

        xPath.pop()
        return resonanceData


class BoundaryCondition(enumsModule.Enum):
    """
    Defines allowed values for the 'boundaryCondition' attribute.
    """

    EliminateShiftFunction = enumsModule.auto()
    NegativeOrbitalMomentum = enumsModule.auto()
    Brune = enumsModule.auto()
    Given = enumsModule.auto()


class RMatrix(ancestryModule.AncestryIO):
    """
    RMatrix is a general container for resolved resonance parameters.
    It can handle standard Lane & Thomas R-Matrix, but can also use various approximations
    (Single- and Multi-level Breit Wigner plus Reich-Moore).

    Internally, resonances are sorted into spin groups, each with a conserved total angular momentum and parity.
    """

    moniker = 'RMatrix'
    ancestryMembers = ('resonanceReactions', 'spinGroups', 'PoPs', 'documentation')
    keyName = 'label'

    _requiredAttributes = (
        ('label', str),
        ('approximation', str),
    )
    _optionalAttributes = (
        ('boundaryCondition', BoundaryCondition, BoundaryCondition.EliminateShiftFunction),
        ('boundaryConditionValue', float, None),
        ('calculateChannelRadius', bool, False),
        ('calculatePenetrability', bool, True),
        ('useForSelfShieldingOnly', bool, False),
        ('supportsAngularReconstruction', bool, False),
        ('LValuesNeededForAngularConvergence', int, None),  # not supported by GNDS-2.1
        ('relativisticKinematics', bool, False),  # not supported by GNDS-2.0
        ('reducedWidthAmplitudes', bool, False),  # not supported by GNDS-2.0
        ('calculateShift', bool, False),  # not supported by GNDS-2.0
    )
    _writeIfDefault = ('boundaryCondition',)

    class Approximation(enumsModule.Enum):
        """Defines the allowed approximations for an RMatrix instance."""

        ReichMoore = 'ReichMoore'
        RMatrix = 'FullRMatrix'

        @classmethod
        def checkEnumOrString(cls, value):
            """Some GNDS-2.0 files were published with incorrect 'approximation' strings. Support reading those."""
            value = {
                'Reich_Moore': cls.ReichMoore,
                'Full R-Matrix': cls.RMatrix
            }.get(value, value)
            return super().checkEnumOrString(value)

    def __init__(self, label, approximation, resonanceReactions, spinGroups, PoPs=None, **kwargs):
        """
        Container for R-Matrix resonance parameters

        :param label: unique label corresponding to a style (i.e., 'eval')
        :param approximation: one of 'Reich_Moore' or 'Full R-Matrix'
        :param resonanceReactions: resonanceReactions instance
        :param spinGroups: SpinGroups instance
        :param PoPs: optional PoPs Database instance
        :param kwargs: see RMatrix._optionalAttributes for supported options
        """

        ancestryModule.AncestryIO.__init__(self)
        self.label = label
        self.approximation = RMatrix.Approximation.checkEnumOrString(approximation)
        self.resonanceReactions = resonanceReactions
        self.spinGroups = spinGroups
        self.PoPs = PoPs

        boundaryCondition = kwargs.get('boundaryCondition', BoundaryCondition.EliminateShiftFunction)
        kwargs['boundaryCondition'] = BoundaryCondition.checkEnumOrString(boundaryCondition)
        for (attr, Type, default) in self._optionalAttributes:
            if attr in kwargs:
                value = kwargs[attr]
                if value != default:
                    assert type(value) is Type, "Incorrect type provided for attribute %s" % (attr)
                setattr(self, attr, value)
            else:
                setattr(self, attr, default)

        self.__documentation = documentationModule.Documentation()
        self.__documentation.setAncestor(self)

    def __getitem__(self, idx):
        return self.spinGroups[idx]

    def __len__(self):
        return len(self.spinGroups)

    @property
    def documentation(self):
        """Returns the documentation instance."""

        return self.__documentation

    @property
    def resonanceReactions(self):
        return self.__resonanceReactions

    @resonanceReactions.setter
    def resonanceReactions(self, value):
        if not isinstance(value, ResonanceReactions):
            raise TypeError("Must be a resonanceReactions instance")
        value.setAncestor(self)
        self.__resonanceReactions = value

    @property
    def spinGroups(self):
        return self.__spinGroups

    @spinGroups.setter
    def spinGroups(self, value):
        if not isinstance(value, SpinGroups):
            raise TypeError("Must be a SpinGroups instance")
        value.setAncestor(self)
        self.__spinGroups = value

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

    def check(self, info):
        return _resonance_checker(self, info, [self.resonanceReactions] + list(self.spinGroups))

    def convertUnits(self, unitMap):

        if self.__PoPs is not None:
            self.__PoPs.convertUnits(unitMap)
        for reac in self.resonanceReactions:
            reac.convertUnits(unitMap)
        for sg in self.spinGroups:
            sg.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + '  '
        xmlString = ['%s<%s label="%s" approximation="%s"' % (indent, self.moniker, self.label, self.approximation)]
        for attr, Type, default in self._optionalAttributes:
            if attr == "LValuesNeededForAngularConvergence" and kwargs.get(
                    'formatVersion', GNDS_formatVersionModule.default) < GNDS_formatVersionModule.version_2_2:
                continue
            if attr in self._writeIfDefault or getattr(self, attr) != default:
                attrVal = getattr(self, attr)
                if Type is bool: attrVal = str(attrVal).lower()
                xmlString[0] += ' %s="%s"' % (attr, attrVal)
        xmlString[0] += '>'

        xmlString += self.__documentation.toXML_strList(indent=indent2, **kwargs)

        if self.__PoPs is not None:
            xmlString += self.__PoPs.toXML_strList(indent=indent2, **kwargs)
        xmlString += self.resonanceReactions.toXML_strList(indent=indent2, **kwargs)
        xmlString += self.spinGroups.toXML_strList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        pops = element.find(PoPsDatabaseModule.Database.moniker)
        if pops is not None:
            pops = PoPsDatabaseModule.Database.parseNodeUsingClass(pops, xPath, linkData, **kwargs)

        RRs = ResonanceReactions()
        RRs.parseNode(element.find(ResonanceReactions.moniker), xPath, linkData, **kwargs)

        SGs = SpinGroups()
        SGs.parseNode(element.find(SpinGroups.moniker), xPath, linkData, **kwargs)

        renames = ()
        formatVersion = kwargs.get( 'formatVersion', GNDS_formatVersionModule.default )
        if formatVersion in [GNDS_formatVersionModule.version_1_10, GNDS_formatVersionModule.version_2_0_LLNL_3, GNDS_formatVersionModule.version_2_0_LLNL_4]:
            renames = (('boundaryCondition',
                       {'S': BoundaryCondition.EliminateShiftFunction,
                        '-L': BoundaryCondition.NegativeOrbitalMomentum}),)

        attrs = getAttrs(element, required=RMatrix._requiredAttributes, optional=RMatrix._optionalAttributes,
                         attributeRenames=renames)

        tmp = cls(attrs.pop('label'), attrs.pop('approximation'), RRs, SGs, PoPs=pops, **attrs)

        documentation = element.find(documentationModule.Documentation.moniker)
        if documentation is not None:
            tmp.documentation.parseNode(documentation, xPath, linkData, **kwargs)

        xPath.pop()
        return tmp


class Channels(suitesModule.Suite):
    """
    Stores a list of channels (used to override global definitions for a single SpinGroup)
    """

    moniker = 'channels'

    def __init__(self):
        suitesModule.Suite.__init__(self, [Channel])


class Channel(ancestryModule.AncestryIO):
    """
    Defines an open channel for a single R-Matrix spin group.
    May be used to override the scattering radius or hard-sphere radius, etc. for this channel.
    """

    moniker = 'channel'
    keyName = 'label'
    ancestryMembers = ('externalRMatrix', 'scatteringRadius', 'hardSphereRadius')
    _requiredAttributes = (
        ('label', str),
        ('resonanceReaction', str),
        ('L', int),
        ('channelSpin', Spin),
        ('columnIndex', int),
    )
    _optionalAttributes = (
        ('boundaryConditionValue', float, None),
    )

    def __init__(self, label, resonanceReaction, L, channelSpin, columnIndex, *,
                 scatteringRadius=None, hardSphereRadius=None, externalRMatrix=None,
                 boundaryConditionValue=None):

        ancestryModule.AncestryIO.__init__(self)
        self.label = label
        self.resonanceReaction = resonanceReaction
        self.L = L
        self.channelSpin = channelSpin
        self.columnIndex = columnIndex
        self.boundaryConditionValue = boundaryConditionValue

        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        self.externalRMatrix = externalRMatrix

    @property
    def scatteringRadius(self):

        return self.__scatteringRadius

    def getScatteringRadius(self):
        """Return ScatteringRadius, looking up ancestry if necessary."""
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            return self.findClassInAncestry(RMatrix).resonanceReactions[self.resonanceReaction].getScatteringRadius()

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """Can be set to None or to a scatteringRadius instance."""

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
        if self.hardSphereRadius is not None:
            return self.hardSphereRadius
        else:
            resonanceReaction = self.findClassInAncestry(RMatrix).resonanceReactions[self.resonanceReaction]
            if resonanceReaction.hardSphereRadius is not None:
                return resonanceReaction.hardSphereRadius

            resonancesHSR = self.findClassInAncestry(Resonances).hardSphereRadius
            if resonancesHSR is not None:
                return resonancesHSR

            # no hardSphereRadius found
            return self.getScatteringRadius()

    @hardSphereRadius.setter
    def hardSphereRadius(self, value):
        """Can be set to None or to a HardSphereRadius instance."""

        if value is not None:
            if not isinstance(value, HardSphereRadius):
                raise TypeError("Hard sphere radius can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__hardSphereRadius = value

    @property
    def externalRMatrix(self):

        return self.__externalRMatrix

    @externalRMatrix.setter
    def externalRMatrix(self, value):
        """Can be set to None or to an ExternalRMatrix instance."""

        if value is not None:
            if not isinstance(value, externalRMatrixModule.ExternalRMatrix):
                raise TypeError("External RMatrix can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__externalRMatrix = value

    def convertUnits(self, unitMap):
        for child in ('scatteringRadius', 'hardSphereRadius', 'externalRMatrix'):
            if getattr(self, child) is not None:
                getattr(self, child).convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + '  '
        attrs = ""
        if self.boundaryConditionValue is not None:
            attrs += ' boundaryConditionValue="%s"' % self.boundaryConditionValue
        xmlString = ['%s<%s label="%s" resonanceReaction="%s" L="%d" channelSpin="%s" columnIndex="%d"%s>' %
                     (indent, self.moniker, self.label, self.resonanceReaction, self.L, self.channelSpin,
                      self.columnIndex, attrs)]
        if self.externalRMatrix is not None:
            xmlString += self.externalRMatrix.toXML_strList(indent=indent2, **kwargs)
        if self.__scatteringRadius is not None:
            xmlString += self.__scatteringRadius.toXML_strList(indent=indent2, **kwargs)
        if self.__hardSphereRadius is not None:
            xmlString += self.__hardSphereRadius.toXML_strList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)
        attrs = getAttrs(element, required=Channel._requiredAttributes, optional=Channel._optionalAttributes)
        for child in element:
            childClass = {'externalRMatrix': externalRMatrixModule.ExternalRMatrix,
                          'scatteringRadius': ScatteringRadius,
                          'hardSphereRadius': HardSphereRadius}.get(child.tag)
            attrs[child.tag] = childClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        tmp = cls(**attrs)

        xPath.pop()
        return tmp


class SpinGroups(suitesModule.Suite):
    """
    Contains a list of spinGroup nodes
    """

    moniker = 'spinGroups'

    def __init__(self):
        suitesModule.Suite.__init__(self, [SpinGroup])


class SpinGroup(ancestryModule.AncestryIO):
    """
    Single group with same Jpi (conserved). Each spin group contains an AP (scattering radius),
    along with 1 or more resonance widths.
    """

    moniker = 'spinGroup'
    keyName = 'label'
    ancestryMembers = ('channels', 'resonanceParameters')

    _requiredAttributes = (
        ('label', str),
        ('spin', Spin),
        ('parity', Parity)
    )

    def __init__(self, label, spin, parity, channels, resonanceParameters):

        ancestryModule.AncestryIO.__init__(self)
        self.label = label
        self.spin = spin
        self.parity = parity
        self.channels = channels
        self.resonanceParameters = resonanceParameters

    def __getitem__(self, idx):
        return self.resonanceParameters[idx]

    def __len__(self):
        return len(self.resonanceParameters.table)

    def __lt__(self, other):
        """ for sorting spin groups by Jpi. group J values together """
        return (self.spin, self.parity) < (other.spin, other.parity)

    @property
    def channels(self):
        return self.__channels

    @channels.setter
    def channels(self, value):
        if not isinstance(value, Channels):
            raise TypeError("Must be a channels instance")
        value.setAncestor(self)
        self.__channels = value

    @property
    def resonanceParameters(self):
        return self.__resonanceParameters

    @resonanceParameters.setter
    def resonanceParameters(self, value):
        if not isinstance(value, ResonanceParameters):
            raise TypeError("Must be a resonanceParameters instance")
        value.setAncestor(self)
        self.__resonanceParameters = value

    def check(self, info):
        # FIXME dummy method
        from fudge import warning
        warnings = []
        for thing in []:
            if thing is None: continue
            warningList = thing.check(info)
            if warningList:
                warnings.append(warning.Context(thing.moniker, warningList))
        return warnings

    def convertUnits(self, unitMap):
        for chan in self.channels:
            chan.convertUnits(unitMap)
        self.resonanceParameters.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s label="%s" spin="%s" parity="%s">' % (indent, self.moniker, self.label, self.spin, self.parity)]
        xml += self.channels.toXML_strList(indent=indent2, **kwargs)
        if self.resonanceParameters.table.columns:  # need not contain any data
            xml += self.resonanceParameters.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker

        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        parameters = ResonanceParameters.parseNodeUsingClass(element.find(ResonanceParameters.moniker), xPath, linkData, **kwargs)
        chs = Channels()
        chs.parseNode(element.find(Channels.moniker), xPath, linkData, **kwargs)
        attrs = getAttrs(element, required=SpinGroup._requiredAttributes)
        SG = cls(channels=chs, resonanceParameters=parameters, **attrs)
        xPath.pop()

        return SG


def _resonance_checker(self, info, things):
    """
    Common resolved resonance checks

    :param self: a reference to the BreitWigner or RMatrix class
    :param info: the dictionary of checking options, shared by all .check() member functions.
                 For the potential scattering convergence test, the 'potentialScatteringCoverganceCriteria'
                 value is the limit on the ratio of the Lth cross section term to the 0th one.  It is printed
                 as a percent but should be entered as a fraction.
    :param things: list of things to check the obvious way (with a .check() function)
    :return:
    """
    from fudge import warning
    from collections import OrderedDict
    warnings = []

    # Check the member data
    for thing in things:
        if thing is None: continue
        warningList = thing.check(info)
        if warningList:
            warnings.append(warning.Context(thing.moniker, warningList))

    # Setup for the physics checks
    import fudge.processing.resonances.reconstructResonances as rrReconstructModule
    rrReconstructor = rrReconstructModule.getResonanceReconstructionClass(self)(self, verbose=False)

    # Check for mismatched spins
    warnings += rrReconstructor.setResonanceParametersByChannel(warnOnly=True)  # multipleSScheme='NJOY', 'ENDF' or None

    # Check for missing channels via sum rule of statistical weights
    sumgJ = {}
    for cc in rrReconstructor.channels:
        if isinstance(cc, OrderedDict):
            iterateThroughThis = cc  # for SLBW
        else:
            iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            if (c.reaction, c.l) not in sumgJ:
                sumgJ[(c.reaction, c.l)] = 0.0
            sumgJ[(c.reaction, c.l)] += c.gfact
    keys = sorted(sumgJ.keys())
    for rxn, L in keys:
        if abs(abs(sumgJ[(rxn, L)]) - abs(2.0 * L + 1.0)) > 1e-6:
            warnings.append(warning.BadSpinStatisticalWeights(L, sumgJ[(rxn, L)], 2. * L + 1, rxn))

    # setup for checking allowed angular momentum:
    def getSList(Ia, Ib):
        """Get possible spins"""
        smin = abs(Ia - Ib)
        smax = Ia + Ib
        nS = int(smax - smin) + 1
        return [iS + smin for iS in range(nS)]

    LMax = 0
    allSs = []
    rxnList = []

    for cc in rrReconstructor.channels:
        if isinstance(cc, OrderedDict):
            iterateThroughThis = cc  # for SLBW
        else:
            iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            if c.eliminated:
                continue
            if c.channelClass is rrReconstructModule.GAMMACHANNEL and isinstance(self, BreitWigner):
                continue
            LMax = max(c.l, LMax)
            if c.reaction not in rxnList:
                rxnList.append(c.reaction)
            if c.channelClass not in [rrReconstructModule.FISSIONCHANNEL, rrReconstructModule.COMPETITIVECHANNEL]:
                try:
                    spinList = getSList(*rrReconstructor.getParticleSpins(c.reaction))
                    for s in spinList:
                        if s not in allSs:
                            allSs.append(s)
                except:
                    warnings.append(warning.UnknownSpinParity(c.reaction))

    # determine the min & max J allowed
    Jmin = min(allSs)
    Jmax = LMax + max(allSs)
    nJ = int(Jmax - Jmin)

    # Check the allowed angular momenta
    for L in range(0, LMax + 1):
        for iJ in range(0, nJ + 1):
            J = iJ + Jmin
            for rxn in rxnList:
                # FIXME ignore capture for spin parity check unless it's R-Matrix and not an eliminated channel
                if 'ission' not in rxn:
                    try:
                        spinList = getSList(*rrReconstructor.getParticleSpins(rxn))
                    except:
                        warnings.append(warning.UnknownSpinParity(rxn))
                        continue

                    for S in spinList:
                        if J not in rrReconstructModule.getAllowedTotalSpins(L, S, useFactor2Trick=False):
                            continue
                        gotIt = False
                        for cc in rrReconstructor.channels:
                            if isinstance(cc, OrderedDict):
                                iterateThroughThis = cc  # for SLBW
                            else:
                                iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
                            for c in iterateThroughThis:
                                if rxn == c.reaction and \
                                        rrReconstructModule.spins_equal(c.l, L) and \
                                        rrReconstructModule.spins_equal(c.J, J) and \
                                        rrReconstructModule.spins_equal(c.s, S):
                                    gotIt = True
                            if gotIt:
                                break
                        if not (gotIt or 'apture' in rxn):
                            theWarning = warning.MissingResonanceChannel(L, S, J, rxn)
                            if str(theWarning) not in [str(w) for w in warnings]:
                                warnings.append(theWarning)

    # Check for convergence in L
    import numpy
    almostXS = {}
    for cc in rrReconstructor.channels:
        if isinstance(cc, OrderedDict):
            iterateThroughThis = cc  # for SLBW
        else:
            iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            egrid = numpy.array([max(rrReconstructor.lowerBound, rrReconstructor.lowerBound + c.Xi),
                                 rrReconstructor.upperBound])
            phis = rrReconstructor.phiByChannel(c, egrid)
            almostXS.setdefault(c.l, []).extend(list(pow(numpy.sin(phis) / rrReconstructor.k(egrid), 2.0)))
    fom = max(almostXS[max(almostXS.keys())]) / max(almostXS[min(almostXS.keys())])
    fomTarget = info.get('potentialScatteringCoverganceCriteria', 0.001)
    if fom > fomTarget:
        warnings.append(warning.PotentialScatteringNotConverged(c.l, rrReconstructor.upperBound, fom, fomTarget))

    return warnings
