# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import math

from pqu import PQU as PQUModule

from numericalFunctions import pointwiseXY as pointwiseXYClass

from PoPs.chemicalElements import chemicalElement as chemicalElementPoPsModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import base as xDataBaseModule
from xData import axes as axesModule
from xData import link as linkModule
from xData import values as valuesModule
from xData import Ys1d as Ys1dModule
from xData import XYs1d as XYs1dModule
from xData import uncertainties as uncertaintiesModule
from xData import regions as regionsModule
from xData import gridded as griddedModule

from fudge import enums as enumsModule
from fudge.processing import group as groupModule

from fudge import abstractClasses as abstractClassesModule
from fudge import styles as stylesModule

from . import URR_probabilityTables as URR_probabilityTablesModule
from fudge.covariances import enums as covarianceEnumsModule

lowerEps = 1e-8
upperEps = 1e-8


def defaultAxes(energyUnit):

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis(Component.moniker, 0, 'b')
    axes[1] = axesModule.Axis('energy_in', 1, energyUnit)
    return axes


#
# crossSection forms.
#
class BaseCrossSectionForm(abstractClassesModule.Form, abc.ABC):

    keyName = 'label'

    @staticmethod
    def toLinearXYsClass():

        return XYs1d

class Ys1d(BaseCrossSectionForm, Ys1dModule.Ys1d):

    def __init__(self, **kwargs):

        BaseCrossSectionForm.__init__(self)
        Ys1dModule.Ys1d.__init__(self, **kwargs)

    def toPointwise_withLinearXYs(self, **kwargs):

        if 'cls' not in kwargs: kwargs['cls'] = XYs1d
        return Ys1dModule.Ys1d.toPointwise_withLinearXYs(self, **kwargs)

class XYs1d(BaseCrossSectionForm, XYs1dModule.XYs1d):

    mutableYUnit = False

    def __init__(self, **kwargs):

        BaseCrossSectionForm.__init__(self)
        XYs1dModule.XYs1d.__init__(self, **kwargs)

    def changeInterpolation(self, interpolation, accuracy, lowerEps=0, upperEps=0, cls=None):

        def func(x, parameters):

            self, T, subParameters = parameters
            x1, y1, x2, y2, f1, f2 = subParameters
            if (None in (x, x1, x2)) or (x < x1) or (x > x2):
                index = self.lowerIndexBoundingX(x)
                x1, y1 = self[index]
                if index + 1 == len(self): return y1
                x2, y2 = self[index+1]
                f1 = 1 / math.sqrt(x1 - T)
                f2 = 1 / math.sqrt(x2 - T)
                parameters[2] = [x1, y1, x2, y2, f1, f2]

            if x == x1: return y1
            if x == x2: return y2
            f = 1 / math.sqrt(x - T)
            alpha = (f - f1) / (f2 - f1)
            return math.pow(x1 * y1, 1 - alpha) * math.pow(x2 * y2, alpha) / x

        if cls is None: cls = XYs1d
        if self.interpolation == xDataEnumsModule.Interpolation.chargedParticle:
            if interpolation != xDataEnumsModule.Interpolation.linlin:
                raise TypeError('Only "%s" interpolation for conversion from %s' %
                                (xDataEnumsModule.Interpolation.linlin, self.interpolation))
            temp = self.cloneToInterpolation(interpolation)
            parameters = [temp, 0, [None, None, None, None, None, None]]
            return temp.applyFunction(func, parameters, accuracy=accuracy, biSectionMax=-1, checkForRoots=False)

        return XYs1dModule.XYs1d.changeInterpolation(self, interpolation, accuracy,
                                                     lowerEps=lowerEps, upperEps=upperEps, cls=cls)

    def effectiveThreshold(self):
        """
        Some reactions have a cross section with multiple zero values at their beginning. In such a case, the effective threshold is the point
        before the first non-zero cross section value. Otherwise, it is the energy of the first point.
        """

        for index, [x, y] in enumerate(self):
            if y != 0: break

        if index > 0: index -= 1
        return self[index][0]

    def evaluate(self, xValue):

        if self.interpolation == xDataEnumsModule.Interpolation.chargedParticle:
            if xValue > self.domainMax:
                raise Exception('Evaluation of "%s" interpolation not supported above domain' % self.interpolation)
            elif xValue < self.domainMin:
                index = 0
            else:
                index = self.lowerIndexBoundingX(xValue)

            T = 0.0
            x1, y1 = self[index]
            x2, y2 = self[index+1]
            f1 = 1.0 / math.sqrt(x1 - T)
            f2 = 1.0 / math.sqrt(x2 - T)
            fValue = 1.0 / math.sqrt(xValue - T)
            alpha = (fValue - f1) / (f2 - f1)
            return math.pow(x1 * y1, 1.0 - alpha) * math.pow(x2 * y2, alpha) / xValue

        return XYs1dModule.XYs1d.evaluate(self, xValue)

    def heat(self, currentTemperature, newTemperature, massRatio, EMin, lowerlimit=None, upperlimit=None,
             interpolationAccuracy=0.001, heatAllPoints=False, doNotThin=True, heatBelowThreshold=True,
             heatAllEDomain=True, setThresholdToZero=False):
        """
        Returns a linear version of the cross-section heated to 'newTemperature'. If the current temperature of the
        cross-section, given by 'currentTemperature', is greater than 'newTemperature' a raise is executed.
        If lowerlimit is None, it is set to 'oneOverV' except when the reaction is determined to be a threshold reaction,
        then it is set to 'threshold'. Any cross-section with domainMin greater than 2.5e-4 eV is determined to be a
        threshold reaction. If upperlimit is None it is set to 'constant'.
        If heatBelowThreshold is False, then EMin is set to the larger of EMin and self's domainMin.

        The unit of 'currentTemperature' and 'newTemperature' must be the same as the self's domain unit.

        For more information on EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints,
        doNotThin and heatAllEDomain see the module crossSectionAdjustForHeatedTarget.
        """

        import crossSectionAdjustForHeatedTarget as heat

        # FIXME what are we checking for here? That C extension was built? Python2 vs 3 differences?
        import types
        if not hasattr(heat, 'crossSectionAdjustForHeatedTarget') and hasattr(heat, 'heat'):
            from crossSectionAdjustForHeatedTarget import heat

        if not isinstance(heat.crossSectionAdjustForHeatedTarget, types.BuiltinFunctionType):
            from crossSectionAdjustForHeatedTarget import crossSectionAdjustForHeatedTarget as heat

        assert isinstance(heat.crossSectionAdjustForHeatedTarget, types.BuiltinFunctionType)

        dT = newTemperature - currentTemperature
        if abs(dT) <= 1e-2 * newTemperature: dT = 0.
        if dT < 0:
            raise Exception('Current temperature "%s" (in energy units) higher than desired temperature "%s"' %
                            (currentTemperature, newTemperature))

        heated = unheated = self
        if not unheated.nf_pointwiseXY.isInterpolationLinear():
            unheated = self.toPointwise_withLinearXYs(accuracy=1e-5, upperEps=upperEps)
        if (len(unheated) > 0) and (dT > 0.):
            if massRatio == 0.0: raise Exception("Heating with gamma as projectile not supported.")
            if not heatBelowThreshold: EMin = max(unheated.domainMin, EMin)
            if lowerlimit is None:
                lowerlimit = "oneOverV"
                if len(self.axes) > 0:
                    domainMin = PQUModule.PQU(1e-5, 'eV').getValueAs(self.domainUnit)
                    if domainMin < 0.04 * unheated.domainMin:
                        lowerlimit = "threshold"
            if upperlimit is None:
                upperlimit = 'constant'

            heated = heat.crossSectionAdjustForHeatedTarget(
                    massRatio, dT, EMin, unheated, lowerlimit=lowerlimit, upperlimit=upperlimit,
                    interpolationAccuracy=interpolationAccuracy, heatAllPoints=heatAllPoints, doNotThin=doNotThin,
                    heatAllEDomain=heatAllEDomain)
        if (unheated[0][1] == 0.0) and setThresholdToZero: heated[0][1] = 0.0
        heated = XYs1d(data=heated, axes=unheated.axes)
        # more thinning likely needed after crossSectionAdjustForHeatedTarget
        heated = heated.thin(interpolationAccuracy / 2)
        return heated

    def processMultiGroup(self, style, tempInfo, indent):

        from fudge.processing import miscellaneous as miscellaneousModule

        def styleFilter(style):

            if isinstance(style, stylesModule.GriddedCrossSection): return False
            return True

        if tempInfo['verbosity'] > 2: print('%sMulti-grouping XYs1d cross section' % indent)

        crossSectionGrouped = miscellaneousModule.groupOneFunctionAndFlux(
                style, tempInfo, self, styleFilter=styleFilter)
        return groupModule.toMultiGroup1d(Gridded1d, style, tempInfo, self.axes, crossSectionGrouped,
                                          zeroPerTNSL=tempInfo['zeroPerTNSL'])

    def toPointwise_withLinearXYs(self, **kwargs):

        ptw = XYs1dModule.XYs1d.toPointwise_withLinearXYs(self, cls=XYs1d, hh=True, **kwargs)
        ptw.label = self.label
        return ptw

    def toXML_strList(self, indent="", **kwargs):

        def dataToString(self, dummy, indent='', **kwargs):

            maxXFormat = "%.17e"

            valuesPerLine = int(kwargs.get('valuesPerLine', 100))
            if (valuesPerLine % 2) == 1: valuesPerLine += 1
            xFormat = kwargs.get("xFormat", "%14.8e")

            XMLList = []
            Line = []
            i1 = 1
            for value in self:
                if (i1 % 2) != 0:
                    xString = xFormat % value
                    if i1 > 1:
                        if xString == priorXString:
                            if xFormat == maxXFormat:
                                raise ValueError('x-value strings identical: "%s" and "%s"' % (priorXString, xString))
                            kwargs['xFormat'] = maxXFormat
                            return dataToString(self, dummy, indent=indent, **kwargs)
                    Line.append(xString)
                    priorXString = xString
                else:
                    Line.append(xFormat % value)
                if (i1 % valuesPerLine) == 0:
                    XMLList.append(indent + ' '.join(Line))
                    Line = []
                i1 += 1
            if len(Line) > 0: XMLList.append(indent + ' '.join(Line))
            return XMLList

        kwargs['dataToString'] = dataToString
        kwargs['dataToStringParent'] = None
        return XYs1dModule.XYs1d.toXML_strList(self, indent, **kwargs)


class Regions1d(BaseCrossSectionForm, regionsModule.Regions1d):
    """
    This class stores a cross-section in two or more regions, which may have different interpolations.
    Each region must contain at least two points. Each pair of adjacent regions must overlap at exactly one point.
    """

    def __init__(self, **kwargs):

        BaseCrossSectionForm.__init__(self)
        regionsModule.Regions1d.__init__(self, **kwargs)

    def effectiveThreshold(self):
        """
        Some reactions have a cross section with multiple zero values at their beginning. In such a case, the effective threshold is the point
        before the first non-zero cross section value. Otherwise, it is the energy of the first point.
        """

        for region in self:
            if isinstance(region, XYs1d):
                if region.rangeMax != 0.0:
                    return region.effectiveThreshold()
            else:
                raise TypeError('Region of type "%s" not supported.' % type(region))

        return self[-1].domainMax

    def processMultiGroup(self, style, tempInfo, indent):

        linear = self.toPointwise_withLinearXYs(accuracy=1e-5, upperEps=1e-8)
        return linear.processMultiGroup(style, tempInfo, indent)

    def toPointwise_withLinearXYs(self, **kwargs):

        xys = regionsModule.Regions1d.toPointwise_withLinearXYs(self, **kwargs)
        return XYs1d(data=xys, axes=xys.axes)

    def ysMappedToXs(self, cls, grid, label=None):

        start = 0
        allYs = []
        for i1, region in enumerate(self):
            offset, Ys = pointwiseXYClass.ysMappedToXs(region, grid.values)
            if i1 == 0:
                allYs = Ys
                start = offset
            else:
                allYs += Ys[1:]

        allYs += (len(grid.values) - (start + len(allYs))) * [0.0]

        ys1d = cls(Ys=valuesModule.Values(allYs, start=start), axes=self.axes.copy(), label=label)
        ys1d.axes[1] = linkModule.Link2(grid.moniker, instance=grid)
        return ys1d

    @staticmethod
    def allowedSubElements():

        return (XYs1d,)


class Gridded1d(BaseCrossSectionForm, griddedModule.Gridded1d):

    def __init__(self, axes, array, **kwargs):

        BaseCrossSectionForm.__init__(self)
        griddedModule.Gridded1d.__init__(self, axes, array, **kwargs)


class ResonanceLink(linkModule.Link):

    moniker = 'resonances'


class BackgroundTerm(ancestryModule.AncestryIO, abc.ABC):

    def __init__(self, data):
        ancestryModule.AncestryIO.__init__(self)
        self.data = data

    @property
    @abc.abstractmethod
    def moniker(self): pass

    @property
    def data(self): return self.__data

    @data.setter
    def data(self, data):
        if not isinstance(data, (XYs1d, Regions1d)):
            raise TypeError(f"Background cross section must be XYs1d or regions1d, not {type(data)}")
        self.__data = data
        self.__data.setAncestor(self)

    @property
    def XYs1d(self):
        if isinstance(self.__data, XYs1d): return self.__data
        raise KeyError("No XYs1d data found in the %s background section" % self.moniker)

    @property
    def regions1d(self):
        if isinstance(self.__data, Regions1d): return self.__data
        raise KeyError("No regions1d data found in the %s background section" % self.moniker)

    def convertUnits(self, unitMap):
        self.data.convertUnits(unitMap)

    def fixDomains(self, energyMin, energyMax, ancestorReference):
        """
        Calls the **fixDomains** for the **data** member.
        """

        counter = self.data.fixDomains(energyMin, energyMax, xDataEnumsModule.FixDomain.upper)

        return counter

    def toXML_strList( self, indent="", **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        child = element[0]
        if child.tag == XYs1d.moniker:
            data = XYs1d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        elif child.tag == Regions1d.moniker:
            data = Regions1d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else:
            raise ValueError("Encountered unexpected element '%s' in %s" % (child.tag, element.tag))

        background_ = cls(data)

        xPath.pop()

        return background_


class ResolvedRegion(BackgroundTerm):

    moniker = 'resolvedRegion'


class UnresolvedRegion(BackgroundTerm):

    moniker = 'unresolvedRegion'


class FastRegion(BackgroundTerm):

    moniker = 'fastRegion'


class Background(ancestryModule.AncestryIO):

    moniker = 'background'
    ancestryMembers = ('resolvedRegion', 'unresolvedRegion', 'fastRegion')

    def __init__(self, resolvedRegion=None, unresolvedRegion=None, fastRegion=None):

        ancestryModule.AncestryIO.__init__(self)
        self.resolvedRegion = resolvedRegion
        self.unresolvedRegion = unresolvedRegion
        self.fastRegion = fastRegion

    def __iter__(self):
        for term in self.ancestryMembers[:-1]:
            if getattr(self, term) is not None: yield getattr(self, term)

    @property
    def resolvedRegion(self):
        return self.__resolvedRegion

    @resolvedRegion.setter
    def resolvedRegion(self, data):
        if isinstance(data, ResolvedRegion):
            data.setAncestor(self)
        elif data is not None:
            raise TypeError("Expected resolvedRegion instance, got %s instead" % type(data))
        self.__resolvedRegion = data

    @property
    def unresolvedRegion(self):
        return self.__unresolvedRegion

    @unresolvedRegion.setter
    def unresolvedRegion(self, data):
        if isinstance(data, UnresolvedRegion):
            data.setAncestor(self)
        elif data is not None:
            raise TypeError("Expected unresolvedRegion instance, got %s instead" % type(data))
        self.__unresolvedRegion = data

    @property
    def fastRegion(self):
        return self.__fastRegion

    @fastRegion.setter
    def fastRegion(self, data):
        if isinstance(data, FastRegion):
            data.setAncestor(self)
        elif data is not None:
            raise TypeError("Expected fastRegion instance, got %s instead" % type(data))
        self.__fastRegion = data

    @property
    def domainMin(self):
        for term in (self.resolvedRegion, self.unresolvedRegion, self.fastRegion):
            if term is not None:
                return term.data.domainMin

    @property
    def domainMax(self):
        for term in (self.fastRegion, self.unresolvedRegion, self.resolvedRegion):
            if term is not None:
                domainMax = term.data.domainMax
                if term.data.domainUnit != self.domainUnit:
                    domainMax = PQUModule.PQU(domainMax, term.data.domainUnit).getValueAs(self.domainUnit)
                return domainMax

    @property
    def domainUnit(self):
        for term in (self.resolvedRegion, self.unresolvedRegion, self.fastRegion):
            if term is not None:
                return term.data.domainUnit

    def convertUnits(self, unitMap):
        for term in self.ancestryMembers :
            term_ = getattr(self, term)
            if term_ is not None:
                term_.convertUnits(unitMap)

    def fixDomains(self, energyMin, energyMax):
        """
        Calls the **fixDomains** for the *fastRegion* member.
        """

        counter = 0
        if self.__resolvedRegion is not None:
            counter += self.__resolvedRegion.fixDomains(energyMin, energyMax, self.__resolvedRegion)
        if self.__unresolvedRegion is not None:
            counter += self.__unresolvedRegion.fixDomains(energyMin, energyMax, self.__unresolvedRegion)
            if len(self.__unresolvedRegion.data) == 0:
                self.__unresolvedRegion = None
        if self.__fastRegion is not None:
            counter += self.__fastRegion.fixDomains(energyMin, energyMax, self.__fastRegion)
            if len(self.__fastRegion.data) == 0:
                self.__fastRegion = None

        return counter

    def toRegions(self):
        """
        Merge all background terms into a single Regions1d
        :return:
        """

        regions = Regions1d()
        axes = None
        for term in (self.resolvedRegion, self.unresolvedRegion, self.fastRegion):
            if term is None: continue

            if axes is None:
                axes = term.data.axes
            elif term.data.axes != axes:
                raise NotImplementedError("Inconsistent axes/units inside background cross section")

            if isinstance(term.data, XYs1d):
                data = term.data.copy()
                data.label = term.moniker
                regions.append(data)
            elif isinstance(term.data, Regions1d):
                for index, region in enumerate(term.data):
                    data = region.copy()
                    data.label = '%s: %s' % (term.moniker, index)
                    regions.append(data)
        regions.axes = axes
        return regions

    def toPointwise_withLinearXYs(self, **kwargs):

        regions = self.toRegions()
        xys = regionsModule.Regions1d.toPointwise_withLinearXYs(regions, **kwargs)
        XYs1d(data=xys, axes=regions.axes)
        return XYs1d(data=xys, axes=regions.axes)

    def toXML_strList(self, indent="", **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        xml = ['%s<%s>' % (indent, self.moniker)]
        for child in (self.resolvedRegion, self.unresolvedRegion, self.fastRegion):
            if child is not None:
                xml += child.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        resolved = unresolved = fast = None
        for child in element:
            if child.tag == ResolvedRegion.moniker:
                resolved = ResolvedRegion.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == UnresolvedRegion.moniker:
                unresolved = UnresolvedRegion.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == FastRegion.moniker:
                fast = FastRegion.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise ValueError("Encountered unexpected element '%s' in %s" % (child.tag, element.tag))

        background_ = cls(resolved, unresolved, fast)

        xPath.pop()

        return background_


class ResonancesWithBackground(BaseCrossSectionForm):
    """
    This class stores cross-sections that include resonances along with a background contribution.
    Contains a link to the resonances, and the 'background' which consists of up to three terms:
    resolved, unresolved and fast regions.
    The full XYs1d cross-section can be obtained by first reconstructing resonances
    and then adding the background contribution (users should use the reactionSuite.reconstructResonances method).
    """

    moniker = 'resonancesWithBackground'
    ancestryMembers = ('resonances', 'background', 'uncertainty')

    def __init__(self, label, resonances, background, uncertainty=None):

        BaseCrossSectionForm.__init__(self)

        if not isinstance(label, str): raise TypeError('label must be a string')
        self.__label = label

        self.resonances = resonances
        self.background = background
        self.uncertainty = uncertainty

    @property
    def resonances(self): return self.__resonances

    @resonances.setter
    def resonances(self, value):
        if not isinstance(value, ResonanceLink):
            raise TypeError("Expected ResonanceLink instance, got %s instead" % type(value))
        self.__resonances = value
        self.__resonances.setAncestor(self)

    @property
    def background(self): return self.__background

    @background.setter
    def background(self, value):
        if not isinstance(value, Background):
            raise TypeError("Expected background instance, got %s instead" % type(value))
        self.__background = value
        self.__background.setAncestor(self)

    @property
    def uncertainty(self):
        return self.__uncertainty

    @uncertainty.setter
    def uncertainty(self, data):
        if isinstance(data, uncertaintiesModule.Uncertainty):
            data.setAncestor(self)
        elif data is not None:
            raise TypeError("Expected uncertainty instance, got %s insteand" % type(data))
        self.__uncertainty = data

    @property
    def label(self): return self.__label

    @property
    def domainMin(self):

        return self.background.domainMin

    @property
    def domainMax(self):

        return self.background.domainMax

    @property
    def domainUnit(self):

        return self.background.domainUnit

    def backgroundAsRegions1d(self):
        '''This methods returns the backgroup cross sections as a Regions1d instance.'''

        return self.background.toRegions()
        
    def convertUnits(self, unitMap):
        """See documentation for reactionSuite.convertUnits."""

        self.background.convertUnits(unitMap)
        # FIXME skipping uncertainties for now

# BRB FIXME This is a kludge until we fix production/crossSection/reference
# For example, see n + U236 in ENDF/B-7.1.
# <productions>
#    <production label="0" outputChannel="U235" ENDF_MT="16">
#      <crossSection>
#        <reference label="eval" ...
    def processMultiGroup(self, style, tempInfo, indent):

        return self.ancestor.processMultiGroup(style, tempInfo, indent)

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        Calls the **fixDomains** for the *background* member.
        """

        return self.background.fixDomains(energyMin, energyMax)

    def toPointwise_withLinearXYs(self, **kwargs):

        _component = self.findClassInAncestry(Component)
        reconstructed = _component.getStyleOfClass(stylesModule.CrossSectionReconstructed)
        if reconstructed is self:
            reconstructed = None  # Happens when there is no reconstructed data, probably because reconstruction failed.
        if reconstructed is not None:
            return reconstructed.toPointwise_withLinearXYs(**kwargs)     # Return a copy.
        raise Exception('resonancesWithBackground cross section has not been reconstructed via reactionSuite.reconstructResonances: %s' % self.toXLink())

    def  asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):

        return self.toPointwise_withLinearXYs(accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps, biSectionMax=biSectionMax)

    def toXML_strList(self, indent="", **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlList = ['%s<%s label="%s">' % (indent, self.moniker, self.label)]
        xmlList.append(self.resonances.toXML(indent2, **kwargs))
        xmlList += self.background.toXML_strList(indent2, **kwargs)
        if self.uncertainty is not None:
            xmlList += self.uncertainty.toXML_strList(indent2, **kwargs)
        xmlList[-1] += '</%s>' % self.moniker
        return xmlList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        link = ResonanceLink.parseNodeUsingClass(element.find(ResonanceLink.moniker), xPath, linkData, **kwargs)
        background = Background.parseNodeUsingClass(element.find(Background.moniker), xPath, linkData, **kwargs)
        uncertainty = None
        uncertaintyNode = element.find(uncertaintiesModule.Uncertainty.moniker)
        if uncertaintyNode is not None:
            uncertainty = uncertaintiesModule.Uncertainty()
            uncertainty.parseNode(uncertaintyNode, xPath, linkData, **kwargs)
        resWithBack = cls(element.get('label'), link, background, uncertainty)

        xPath.pop()

        return resWithBack


class Reference(linkModule.Link, BaseCrossSectionForm):
    """This cross section form consists of a reference to another cross section."""

    moniker = 'reference'

    def __init__(self, link=None, root=None, path=None, label=None, relative=False):

        linkModule.Link.__init__(self, link=link, root=root, path=path, label=label, relative=relative)
        BaseCrossSectionForm.__init__(self)

    @property
    def crossSection(self):

        if self.link is None: raise Exception("Unresolved link!")
        return self.link

    @property
    def domainMin(self):

        return self.crossSection.domainMin

    @property
    def domainMax(self):

        return self.crossSection.domainMax

    @property
    def domainUnit(self):

        return self.crossSection.domainUnit

    def convertUnits(self, unitMap):
        """See documentation for reactionSuite.convertUnits."""

        pass

    def fixDomains(self, labels, energyMin, energyMax):
        """This method does nothing."""

        return 0

    def getReference(self):

        return self.crossSection

    def processMultiGroup(self, style, tempInfo, indent):

        addToComponent = tempInfo.get('addToComponent', None)
        tempInfo['addToComponent'] = False
        multiGroup = self.crossSection.processMultiGroup(style, tempInfo, indent)
        tempInfo['addToComponent'] = addToComponent
        if addToComponent is None: del tempInfo['addToComponent']
        return multiGroup

    def toPointwise_withLinearXYs(self, **kwargs):

        return self.crossSection.toPointwise_withLinearXYs(**kwargs)

    def check(self):

        from fudge import warning
        warnings = []
        if self.rootAncestor != self.getReference().rootAncestor:
            warnings.append(warning.BadCrossSectionReference())
        return warnings


class CoulombPlusNuclearElastic(Reference):
    """Special type of link to doubleDifferentialCrossSection/chargedParticleElastic."""

    moniker = "CoulombPlusNuclearElastic"


class ThermalNeutronScatteringLaw1d(Reference):
    """Special type of link to doubleDifferentialCrossSection/thermalNeutronScatteringLaw."""

    moniker = 'thermalNeutronScatteringLaw1d'


class URR_probabilityTables1d(BaseCrossSectionForm, xDataBaseModule.XDataFunctional):
    """
    FIXME: name is misleading, these are actually cross section pdfs/cdfs rather than probability tables.
    Stores P(sigma | E) and CDF(sigma | E) in an XYs2d container, where incident energy E spans the unresolved region.
    """

    moniker = 'URR_probabilityTables1d'
    dimension = 1

    def __init__(self, label, URR_probabilityTables):

        BaseCrossSectionForm.__init__(self)
        # FIXME Doesn't really make sense to have this inherit XDataFunctional...
        xDataBaseModule.XDataFunctional.__init__(self, axes=None, label=label)

        if not isinstance(URR_probabilityTables, URR_probabilityTablesModule.XYs2d):
            raise TypeError('Invalid URR probability tables instance.')
        self.__data = URR_probabilityTables

    @property
    def data(self):

        return self.__data

    def convertUnits(self, unitMap):
        """
        Converts each axis units.
        unitMap is a dictionary of mapping old units to new units (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        self.__data.convertUnits(unitMap)

    def toXML_strList(self, indent="", **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlList = ['%s<%s label="%s">' % (indent, self.moniker, self.label)]
        xmlList += self.__data.toXML_strList(indent=indent2, **kwargs)
        xmlList[-1] += '</%s>' % self.moniker

        return xmlList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        items = {}
        for child in element:
            for cls2 in [URR_probabilityTablesModule.XYs2d]:
                if child.tag == cls2.moniker:
                    items[child.tag] = cls2.parseNodeUsingClass(child, xPath, linkData, **kwargs)
                    break

        function2d = items[list(items.keys())[0]]
        URR_probability_tables1d = cls(element.get('label'), function2d)

        xPath.pop()

        return URR_probability_tables1d


def chargeParticle_changeInterpolationSubFunction(threshold, x, x1, y1, x2, y2):

    B = math.log(x2 * y2 / (x1 * y1)) / (1. / math.sqrt(x1 - threshold) - 1. / math.sqrt(x2 - threshold))
    A = x1 * y1 * math.exp(B / math.sqrt(x1 - threshold))
    y = A * math.exp(-B / math.sqrt(x - threshold)) / x
    return y


#
# crossSection component
#
class Component(abstractClassesModule.Component):

    moniker = 'crossSection'

    def __init__(self):

        abstractClassesModule.Component.__init__(
                self, allowedClasses=(Ys1d, XYs1d, Regions1d, Gridded1d, ResonancesWithBackground,
                                      Reference, URR_probabilityTables1d, CoulombPlusNuclearElastic,
                                      ThermalNeutronScatteringLaw1d))

    def domainUnitConversionFactor(self, unitTo):

        if unitTo is None: return 1.
        return PQUModule.PQU('1 ' + self.domainUnit).getValueAs(unitTo)

    @property
    def domainMin(self):

        return self.evaluated.domainMin

    @property
    def domainMax(self):

        return self.evaluated.domainMax

    @property
    def domainUnit(self):

        return self.evaluated.domainUnit

    def diff(self, other, diffResults, **kwargs):

        threshold_epsilon = kwargs.get('threshold_epsilon', 1e-6)
        lowerEps = kwargs.get('lowerEps', 1e-8)
        upperEps = kwargs.get('upperEps', 1e-8)

        crossSection1 = self.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)
        crossSection2 = other.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)
        if (abs(crossSection2.domainMin - crossSection1.domainMin) >
                threshold_epsilon * (crossSection1.domainMin + crossSection1.domainMin)):
            relativeError = ((crossSection2.domainMin - crossSection1.domainMin) /
                             (crossSection1.domainMin + crossSection1.domainMin))
            diffResults.append('Cross section thresholds differ', 'relative error = %.6e: %.12e vs %.12e %s' %
                               (relativeError, crossSection1.domainMin, crossSection2.domainMin,
                                crossSection1.domainUnit), self.toXLink(), other.toXLink())

    def effectiveThreshold(self):
        """
        Some reactions have a cross section with multiple zero values at their beginning.
        In such a case, the effective threshold is the point before the first non-zero cross-section value.
        If an effective threshold cannot be determined, then domainMin is returned.
        """

        evaluated = self.evaluated
        if hasattr(evaluated, 'effectiveThreshold'):
            try:
                return evaluated.effectiveThreshold()
            except:
                pass

        return self.domainMin

    def hasLinearForm(self):

        for form in self:
            if isinstance(form, XYs1d):
                if form.nf_pointwiseXY.isInterpolationLinear(): return form
        return None

    def heat(self, style, EMin, lowerlimit=None, upperlimit=None, interpolationAccuracy=0.001,
             heatAllPoints=False, doNotThin=True, heatBelowThreshold=True, heatAllEDomain=True,
             setThresholdToZero=False, addToSuite=False):
        """
        Returns the result of self.toPointwise_withLinearXYs( ).heat( ... ).
        See method crossSection.XYs1d.heat for more information.
        If setThresholdToZero is True and self's cross section at the first point is 0., the heated cross-section's
        first value will also be 0.
        """

        styles = self.findAttributeInAncestry('styles')
        currentstyle = styles[style.derivedFrom]
        if not isinstance(currentstyle, (stylesModule.Evaluated, stylesModule.Heated)):
            TypeError('Form to heat is not heatable form: its style moniker is "%s"' % style.label)
        currentTemperature = PQUModule.PQU(currentstyle.temperature.getValueAs('K'), 'K') * PQUModule.PQU('1 k')
        currentTemperature = currentTemperature.getValueAs(self.domainUnit)

        newTemperature = PQUModule.PQU(style.temperature.getValueAs('K'), 'K') * PQUModule.PQU('1 k')
        newTemperature = newTemperature.getValueAs(self.domainUnit)

        reactionSuite = self.rootAncestor
        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        projectileMass = projectile.getMass('amu')

        targetID = reactionSuite.target
        if targetID in reactionSuite.PoPs.aliases: targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]

        if projectileMass == 0.0:
            massRatio = 0.0
        elif isinstance(target, chemicalElementPoPsModule.ChemicalElement):
            massRatio = 0.0
        else:
            massRatio = target.getMass('amu') / projectileMass

        linear = self.toLinear(label=style.derivedFrom, accuracy=1e-5, upperEps=1e-8)
        heated = linear.heat(currentTemperature, newTemperature, massRatio, EMin, lowerlimit, upperlimit,
                             interpolationAccuracy, heatAllPoints, doNotThin, heatBelowThreshold, heatAllEDomain,
                             setThresholdToZero=setThresholdToZero)
        heated.label = style.label
        if addToSuite: self.add(heated)
        return heated

    def check(self, info):
        """
        Check cross section data for correct threshold, negative cross-sections, etc.
        Returns a list of any warnings encountered during checking.

        :param dict info:  A Python dictionary containing the parameters that control the cross-section checking.
        :keyword boolean CoulombReaction: True if target and projectile are both charged particles,
                or if two or more products are charged particles.
        :keyword float Q: a parameter of `info`: if Q is positive (and CoulombReaction=False), cross-section
                must start at crossSectionEnergyMin, otherwise cross-section threshold must agree with Q-value
        :keyword float kinematicFactor: a parameter of `info`: equal to (targetMass + projectileMass) / targetMass
        :keyword string dThreshold: a parameter of `info`:
                allowable threshold energy mismatch as a string suitable for the PQU class
        :keyword float crossSectionEnergyMin: a parameter of `info`:
                non-threshold cross-section must start at this limit (usually 1e-5 eV)
        :keyword float crossSectionEnergyMax: a parameter of `info`:
                the cross-section must extend up to limit (usually 20 MeV)
        """

        from fudge import warning

        warnings = []

        # FIXME: domain is giving incorrect answers for threshold reactions with resonance contributions!
        lower = PQUModule.PQU(self.domainMin, self.domainUnit)
        upper = PQUModule.PQU(self.domainMax, self.domainUnit)
# BRB6 hardwired
        thresh = PQUModule.PQU(0, 'eV')
        # FIXME skipping Q-value vs. threshold check for continuum channels (e.g. MT=91) since GNDS only stores final Q
        if 'Q' in info and not info['ContinuumOutputChannel']:
            thresh = -info['Q'] * info['kinematicFactor']
            if not isinstance(info['Q'], PQUModule.PQU): thresh = PQUModule.PQU(thresh, 'eV')

            if not info['CoulombChannel']:
                # if Q is positive, cross-section must start at info['crossSectionEnergyMin'] (usually 1e-5 eV)
                # otherwise, cross-section threshold must agree with Q-value:
                if thresh.value > 0:
                    if abs(thresh-lower) > PQUModule.PQU(info['dThreshold']):
                        warnings.append(warning.Threshold_mismatch(lower, thresh, self))
                elif lower != PQUModule.PQU(info['crossSectionEnergyMin']):
                    reaction = self.ancestor
                    hasFissionGenre = False
                    if hasattr(reaction, 'fissionGenre'):       # Ignore 2nd, 3rd and 4th-chance fission (they have Q>0 but are still threshold reactions).
                        hasFissionGenre = reaction.fissionGenre in (enumsModule.FissionGenre.none, enumsModule.FissionGenre.total, 
                                enumsModule.FissionGenre.firstChance)
                    if not hasattr(reaction, 'outputChannel') or hasFissionGenre:
                        warnings.append(warning.Threshold_mismatch(lower, PQUModule.PQU(info['crossSectionEnergyMin']), self))
            else:
                # charged-particle reaction generally doesn't 'turn on' right at threshold due to Coulomb barrier.
                # In this case, just ensure the cross-section stays zero until at or above threshold:
                if lower < thresh:
                    warnings.append(warning.Coulomb_threshold_mismatch(lower, thresh, self))

        # cross-section must extend up to limit (usually 20 MeV):
        if upper < PQUModule.PQU(info['crossSectionEnergyMax']):
            warnings.append(warning.GapInCrossSection(upper, PQUModule.PQU(info['crossSectionEnergyMax']), self))

        evaluatedCrossSection = self.evaluated
        if isinstance(evaluatedCrossSection, ResonancesWithBackground):
            # ensure no gaps between resonance and background:
            resonances = info['reactionSuite'].resonances
            if resonances.resolved is not None:
                resDomainMax = resonances.resolved.domainMax
            if resonances.unresolved is not None and not resonances.unresolved.evaluated.useForSelfShieldingOnly:
                resDomainMax = resonances.unresolved.domainMax
            bckDomainMin = evaluatedCrossSection.background.domainMin
            if bckDomainMin > resDomainMax:
                warnings.append(warning.GapInCrossSection(resDomainMax, bckDomainMin, self))
            linearCrossSection = self[info['reconstructedStyleName']]

        elif (isinstance(evaluatedCrossSection, Reference) and
              isinstance(evaluatedCrossSection.link, ResonancesWithBackground)):
            linearCrossSection = evaluatedCrossSection.link.ancestor[info['reconstructedStyleName']]

        elif isinstance(evaluatedCrossSection, CoulombPlusNuclearElastic):
            return warnings     # already checked as part of doubleDifferentialCrossSection

        else:
            # can it be converted to XYs1d linear?
            try:
                linearCrossSection = evaluatedCrossSection.toPointwise_withLinearXYs(accuracy=1e-5, upperEps=1e-8)
            except Exception as e:
                warnings.append( warning.ExceptionRaised(e, self))
                if info['failOnException']: raise
                return warnings

        # test for negative values, and for non-zero cross-section at threshold
        if linearCrossSection.rangeMin < 0:
            for i, (en, xsc) in enumerate(linearCrossSection):
                if xsc < 0:
                    warnings.append(
                        warning.NegativeCrossSection(PQUModule.PQU(en, linearCrossSection.axes[-1].unit), i, self))

        if thresh.value > 0:
            if linearCrossSection[0][1] != 0:
                warnings.append(warning.NonZero_crossSection_at_threshold(linearCrossSection[0][1], self))

        interpolations = []
        if isinstance(evaluatedCrossSection, XYs1d):
            interpolations = [evaluatedCrossSection.interpolation]
        elif isinstance(evaluatedCrossSection, Regions1d):
            interpolations = [region.interpolation for region in evaluatedCrossSection]
        if any([interp is xDataEnumsModule.Interpolation.flat for interp in interpolations]):
            warnings.append(warning.CrossSectionFlatInterpolation())

        return warnings

    def findEntity(self, entityName, attribute=None, value=None):

        for form in self:
            if form.moniker == entityName: return form
        return ancestryModule.Ancestry.findEntity(self, entityName, attribute, value)

    def evaluate(self, energyIn):

        if (energyIn < self.domainMin) or (energyIn > self.domainMax): return 0.0

        pointwise = self.hasLinearForm()
        if pointwise is None: pointwise = self.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)
        return pointwise.evaluate(energyIn)

    def evaluateWithUncertainty(self, E, useCovariance=True, covariance=None, covarianceSuite=None):
        """

        :param E:
        :param useCovariance: use this to override covarance usage
        :type useCovariance: bool
        :param covariance:
        :param covarianceSuite:
        :return:
        """

        if not isinstance(E, PQUModule.PQU): raise TypeError( "E must be an PQUModule.PQU instance")

        pointwise = self.hasLinearForm()
        if pointwise is None: pointwise = self.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)

        EinDomainUnit = E.getValueAs(pointwise.domainUnit)
        if (EinDomainUnit < pointwise.domainMin) or (EinDomainUnit > pointwise.domainMax):
            meanValue = 0.0
        else:
            meanValue = pointwise.evaluate(EinDomainUnit)

        # We might be done
        if not useCovariance: return PQUModule.PQU(meanValue, unit=pointwise.rangeUnit)

        # Get that the covariance goes with the data.
        covariance = self.getMatchingCovariance(covariance, covarianceSuite)

        if covariance is None:
            return PQUModule.PQU(meanValue, unit=pointwise.rangeUnit)
        else:
            try:
                uncertainty = covariance.getUncertaintyVector(self.evaluated, relative=False).evaluate(EinDomainUnit)
                return PQUModule.PQU(meanValue, unit=pointwise.rangeUnit, uncertainty=uncertainty)
            except IndexError as err:  # FIXME: a kludge until cov routines working again, try it on 27Al(n,tot)
                print("WARNING: Could not extract uncertainty, got error %s" % str(err))
                return PQUModule.PQU(meanValue, unit=pointwise.rangeUnit)

#        if covariance is None:
#            if hasattr(self.evaluated, 'uncertainty') and self.evaluated.uncertainty:
#                if hasattr(self.evaluated.uncertainty.data,'link'):
#                    covariance = self.evaluated.uncertainty.data.link['eval']
#                elif covarianceSuite is not None:
#                    covariance = self.evaluated.uncertainty.data.follow(startNode=covarianceSuite)['eval']
#                else: return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit() )
#            else: return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit() )
#        else:
#            if covariance is not( isinstance( covariance, covModule.CovarianceMatrix ) ):
#                raise TypeError( 'covariance must be of type CovarianceMatrix, got %s'%str(type(covariance)))

#        # Compute uncertainty on convolution given mean value and covariance.
#        try:
#            uncertainty = covariance.getUncertaintyVector( self.evaluated, relative = False ).evaluate( EinDomainUnit )
#        except:
#            print("WARNING: could not get uncertainty in evaluateWithUncertainty for %s" % str(covariance.__class__))
#            return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit() )
#        if uncertainty is None:

#        # We haven't changed units on uncertainty or the meanValue, so we can reconstruct the physical quantity easily.
#        return( PQUModule.PQU( meanValue, unit = ptwise.rangeUnit(), uncertainty = uncertainty.evaluate(EinDomainUnit) ) )

    def getMatchingCovariance(self, covariance = None, covarianceSuite = None):
        """
        Retrieve the uncertainty that goes with this cross-section, if we can find it
        :return:
        """
        import fudge.covariances.covarianceMatrix as covModule

        # Get the covariance that goes with the data.
        if covariance is None:
            # FIXME ugly logic here since cross section 'forms' handle uncertainty differently. Need to make them consistent!
            if (hasattr(self.evaluated, 'uncertainty') and self.evaluated.uncertainty is not None and
                    self.evaluated.uncertainty.data is not None):
                if (isinstance(self.evaluated.uncertainty.data, linkModule.Link) and
                        self.evaluated.uncertainty.data.linkWithoutUpdating is not None):
                    covariance = self.evaluated.uncertainty.data.link
                elif covarianceSuite is not None:
                    if self.evaluated.uncertainty.data is not None:
                        covariance = self.evaluated.uncertainty.data.follow(startNode=covarianceSuite)
        else:
            if not isinstance(covariance, covModule.CovarianceMatrix):
                raise TypeError('covariance must be of type CovarianceMatrix, got %s' % type(covariance))

        return covariance

    def integrateTwoFunctionsWithUncertainty(self, f2, domainMin=None, domainMax=None, useCovariance=True,
                                             covariance=None, covarianceSuite=None, normalize=False):
        r"""
        Computes the spectrum (i.e. weighted) integral of self with the spectrum (i.e. weight)
        specified by ``spectrum``  optionally from Emin to Emax.  If the covariance is
        provided, the uncertainty on the spectrum integral will be computed.

        :param f2: spectrum over which to average
        :type f2: XYs1d instance
        :param domainMin: Lower integration limit.
            If None (the default), then the lower limit of self's domain is used
        :type domainMin: PQU or None
        :param domainMax: Upper integration limit.
            If None (the default), then the upper limit of self's domain is used
        :type domainMax: PQU or None
        :param useCovariance: use this to override covarance usage
        :type useCovariance: bool
        :param covariance: covariance to use when computing uncertainty on the spectral average.
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :param covarianceSuite: covarianceSuite to use for looking up covariance if necessary
        :type covarianceSuite: CovarianceSuite
        :param normalize: if set to True, normalize integral of the spectrum over the interval so we are doing a spectrum average
        :type normalize: bool
        :rtype: PQU


        :How does it work?:
            Given a  weighting spectrum :math:`\phi(E)`, we intend to average the cross section in self as follows:

            .. math::
                < \sigma > = \int_{E_{min}}^{E_{max}} dE \phi(E) \sigma(E)

            To compute the uncertainty, we resort to the basis function expansion of the covariance
            as follows:

            .. math::
                \Delta^2\sigma(E,E') = \sum_{ij} \Delta^2\sigma_{ij} B_i(E) B_j(E')

            In the ENDF format, all covariances are assumed to be grouped in energy so the
            basis functions are simple window functions.   With this,

            .. math::
                \delta <\sigma> = \sqrt{ \int dE dE' \phi(E)\phi(E') \Delta^2\sigma(E,E') }
                                = \sqrt{ \sum_{ij} \Delta^2\sigma_{ij} <B_i> <B_j> }

        """
        # Check that the inputs are of the correct type
        if not isinstance(f2, XYs1dModule.XYs1d): raise TypeError("spectrum must be an XYs1d instance")

        # Convert the cross-section toXYs1d
        ptwise = self.hasLinearForm()
        if ptwise is None: ptwise = self.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)

        # Check domains
        if domainMin is None:
            domainMin = PQUModule.PQU(max(ptwise.domainMin, f2.domainMin), ptwise.domainUnit)
        elif not isinstance(domainMin, PQUModule.PQU):
            raise TypeError("domainMin must be an PQUModule.PQU instance")

        if domainMax is None:
            domainMax = PQUModule.PQU(min(ptwise.domainMax, f2.domainMax), ptwise.domainUnit)
        elif not isinstance(domainMax, PQUModule.PQU):
            raise TypeError( "domainMax must be an PQUModule.PQU instance")

        # Convolve spectrum and self to get mean value
        meanValue = ptwise.integrateTwoFunctions(f2,domainMin,domainMax)
        meanValue = PQUModule.PQU(meanValue, ptwise.axes[0].unit)

        # Normalize the mean if we're averaging
        if normalize and meanValue.value != 0.0:
            norm = f2.integrate(domainMin,domainMax)
            if norm == 0.0:
                raise ValueError('zero norm (%s) while integrating function with %s ' %
                                 (str(norm), str(self.ancestor.ancestor)))
            meanValue = meanValue/norm

        # We might be done
        if not useCovariance: return meanValue # already a PQU

        # Get that the covariance goes with the data.
        covariance = self.getMatchingCovariance(covariance, covarianceSuite)
        if covariance is None: return meanValue

        try:
            theCovariance = covariance.toCovarianceMatrix() # may be redundant, but at least we'll get a usable data type
        except Exception as err:
            print(
                "WARNING: could not get covariance in integrateTwoFunctionsWithUncertainty for form %s, "
                "got error message '%s'" % (str(covariance.__class__), err))
            return meanValue

        # Compute weighting vector from spectrum and possibly the cross-section (if covariance is relative)
        grid = theCovariance.matrix.axes[-1].values.values
        gridUnit = theCovariance.matrix.axes[-1].unit
        covGroupBdries = list(grid)
        if theCovariance.type == covarianceEnumsModule.Type.absolute:
            phi = f2.group(covGroupBdries, norm=None)
        elif theCovariance.type == covarianceEnumsModule.Type.relative:
            phi = f2.group(covGroupBdries, ptwise, norm=None)
        else:
            raise ValueError("Unknown covariance type: %s" % str(type(covariance)))

        # Compute the uncertainty itself
        import numpy
        phi = numpy.asmatrix( phi )
        theArray = theCovariance.matrix.array.constructArray()
        try:
            coco = (phi * theArray * phi.T)[0,0]
            if coco < 0.0: print ("WARNING: covariance of spectrum integral is %s < 0.0" % str(coco))
            uncertainty = PQUModule.PQU(math.sqrt(max(coco, 0.0)), unit=meanValue.unit)
            if normalize and meanValue.value != 0.0: uncertainty = uncertainty/norm

            # it worked, return result
            return PQUModule.PQU(meanValue.value, unit=meanValue.unit, uncertainty=uncertainty.getValueAs(meanValue.unit))
        except Exception as err:
            print("WARNING: could not compute uncertainty in integrateTwoFunctionsWithUncertainty for form %s, "
                  "got error message '%s'" % (str(covariance.__class__), err))
            return meanValue

    def processGriddedCrossSections(self, style, verbosity=0, indent='', incrementalIndent='  ', isPhotoAtomic=False):

        crossSection = style.findFormMatchingDerivedStyle(self)
        ys = crossSection.ysMappedToXs(Ys1d, style.grid, label=style.label, extendToEnd=True)
        if isPhotoAtomic:
            ys.interpolation = xDataEnumsModule.Interpolation.loglog
        if isPhotoAtomic and (ys.Ys.start > 0): ys.Ys.values[0] = 1e-10
        self.add(ys)
