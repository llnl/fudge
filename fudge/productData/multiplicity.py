# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

try :
    from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
    floatToShortestString = pointwiseXY_CModule.floatToShortestString
except :
    from pqu import PQU as PQUModule
    floatToShortestString = PQUModule.floatToShortestString

from PoPs import IDs as IDsPoPsModule

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import constant as constantModule
from xData import XYs1d as XYs1dModule
from xData import series1d as series1dModule
from xData import regions as regionsModule
from xData import gridded as griddedModule
from xData import link as linkModule

from fudge import abstractClasses as abstractClassesModule

from fudge.reactions import base as reactionBaseModule

def defaultAxes( energyUnit ) :

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( 'multiplicity', 0, '' )
    axes[1] = axesModule.Axis( 'energy_in', 1, energyUnit )
    return( axes )

class BaseMultiplicityForm( abstractClassesModule.Form ) :

    keyName = 'label'

    @property
    def product( self ) :

        return( self.ancestor.ancestor )

class Unspecified( BaseMultiplicityForm ) :

    moniker = 'unspecified'

    def __init__( self, label ) :

        BaseMultiplicityForm.__init__( self )

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def convertUnits( self, unitMap ) :

        pass

    def evaluate( self, E ) :

        return( 0. )

    @property
    def label( self ) :

        return( self.__label )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """This method does nothing."""

        return 0

    def toPointwise_withLinearXYs( self, **kwargs ) :

        raise Exception( 'Linear-linear XYs1d data for multiplicity form %s does not make sense' % self.moniker )

    def toXML_strList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s/>' % ( indent, self.moniker, attributeStr ) ] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        val = cls( element.get('label') )

        xPath.pop()

        return val

class Constant1d( BaseMultiplicityForm, constantModule.Constant1d ) :

    def __init__(self, multiplicity, domainMin, domainMax, axes=None, label=None, index=None, outerDomainValue=None):

        BaseMultiplicityForm.__init__( self )
        constantModule.Constant1d.__init__(self, multiplicity, domainMin, domainMax, axes=axes, label=label)

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ EMin, EMax ] )

    def isConstant( self ) :

        return( True )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """This method returns the multiplicity as linear-linear XYs1d data which spans self's domain."""

        return( XYs1d( data = [ [ self.domainMin, self.value ], [ self.domainMax, self.value ] ], axes = self.axes ) )

class XYs1d( BaseMultiplicityForm, XYs1dModule.XYs1d ) :

    def __init__( self, **kwargs ) :

        BaseMultiplicityForm.__init__( self )
        XYs1dModule.XYs1d.__init__( self, **kwargs )

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.domainMin, self.domainMax )

    def evaluate( self, E ) :

        if( E < self[0][0] ) : E = self[0][0]
        if( E > self[-1][0] ) : E = self[-1][0]
        return( XYs1dModule.XYs1d.evaluate( self, E ) )

    def isConstant( self ) :

        return self.rangeMin == self.rangeMax

    def processMultiGroup(self, style, tempInfo, indent):

        from fudge.processing import miscellaneous as miscellaneousModule

        if tempInfo['verbosity'] > 2: print('%sMulti-grouping XYs1d multiplicity' % indent)

        return miscellaneousModule.groupFunctionCrossSectionAndFlux(Gridded1d, style, tempInfo, self)

class Regions1d( BaseMultiplicityForm, regionsModule.Regions1d ) :

    def __init__( self, **kwargs ) :

        BaseMultiplicityForm.__init__( self )
        regionsModule.Regions1d.__init__( self, **kwargs )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ self[0][0][0], self[-1][-1][0] ] )

    def evaluate( self, E ) :

        if( E < self[0][0][0] ) : E = self[0][0][0]
        for region in self :
            if( E < region[-1][0] ) : break
        return( region.evaluate( E ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps =  1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """See regionsXYs.toPointwise_withLinearXYs on the use of lowerEps, upperEps."""

        kwargs['removeOverAdjustedPoints'] = True
        xys = regionsModule.Regions1d.toPointwise_withLinearXYs( self, **kwargs )
        return( XYs1d( data = xys, axes = xys.axes ) )

class Reference( linkModule.Link, BaseMultiplicityForm ) :

    moniker = 'reference'

    def __init__( self, link = None, root = None, path = None, label = None, relative = False ) :

        linkModule.Link.__init__( self, link = link, root = root, path = path, label = label, relative = relative )
        BaseMultiplicityForm.__init__( self )

    @property
    def reference( self ):
        if self.link is None:
            raise Exception("Unresolved link!")
        return self.link

    @property
    def domainMin( self ) :

        return( self.reference.domainMin )

    @property
    def domainMax( self ) :

        return( self.reference.domainMax )

    @property
    def domainUnit( self ) :

        return( self.reference.domainUnit )

    def convertUnits( self, unitMap ) :

        pass

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """This method does nothing."""

        return 0

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.reference.getEnergyLimits( EMin, EMax ) )

    def evaluate( self, E ) :

        return( self.reference.evaluate( E ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.reference.toPointwise_withLinearXYs( **kwargs ) )

class Polynomial1d( BaseMultiplicityForm, series1dModule.Polynomial1d ) :

    def __init__( self, coefficients, domainMin, domainMax, **kwargs ) :

        BaseMultiplicityForm.__init__( self )
        series1dModule.Polynomial1d.__init__( self, coefficients, domainMin, domainMax, **kwargs )

    @property
    def domainUnit( self ) :

        return( self.findClassInAncestry( reactionBaseModule.Base_reaction ).domainUnit )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ EMin, EMax ] )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        xys = series1dModule.Polynomial1d.toPointwise_withLinearXYs( self, accuracy = 1e-6 )
        return( XYs1d( data = xys, axes = self.axes ) )

class PartialProduction( BaseMultiplicityForm ) :

    moniker = 'partialProduction'

    def __init__( self, label, data = 1 ) :

        BaseMultiplicityForm.__init__( self )
        self.data = data

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def evaluate( self, E ) :

        return( self.data )

    @property
    def label( self ) :

        return( self.__label )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        print(type(self.data))
        raise NotImplementedError("partialProduction.toPointwise_withLinearXYs")

    def toXML_strList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s value="%s"/>' % ( indent, self.moniker, attributeStr, self.data ) ] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        label = element.get( 'label' )
        mult = float( element.get( 'value' ) )
        if( mult.is_integer( ) ) : mult = int( mult )
        multiplicity = cls( label, mult )

        xPath.pop( )

        return( multiplicity )

class Gridded1d(BaseMultiplicityForm, griddedModule.Gridded1d):

    def __init__( self, axes, array, **kwargs ) :

        BaseMultiplicityForm.__init__(self)
        griddedModule.Gridded1d.__init__(self, axes, array, **kwargs)

class Branching1d( BaseMultiplicityForm ) :

    moniker = 'branching1d'

    def __init__( self, label ) :

        BaseMultiplicityForm.__init__( self )

        self.__label = label

    @property
    def label( self ) :

        return( self.__label )

    def convertUnits( self, unitMap ) :

        pass

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """The method does nothing."""

        return 0

    def processMultiGroup( self, style, tempInfo, indent ) :

        def countGammas( product, probability ) :

            count = 0.0
            for decayMode in product.decayData.decayModes :
                branchingRatio = decayMode.probability[0].value * probability
                count += branchingRatio

                decayPath = decayMode.decayPath[0]
                residual = decayPath.products[0].pid
                if( residual == IDsPoPsModule.photon ) : residual = decayPath.products[1].pid
                count += countGammas( PoPs[residual], branchingRatio )

            return( count )

        verbosity = tempInfo['verbosity']
        crossSection = tempInfo['crossSection']
        domainMin = crossSection.domainMin
        domainMax = crossSection.domainMax
        
        PoPs = tempInfo['reactionSuite'].PoPs

        nuclide = PoPs[self.product.parentProduct.pid]

        multipliticy = countGammas( nuclide, 1.0 )
        multipliticy = XYs1d( data = [ [ domainMin, multipliticy ], [ domainMax, multipliticy ] ], axes = defaultAxes( tempInfo['incidentEnergyUnit'] ) )

        return( multipliticy.processMultiGroup( style, tempInfo, indent ) )

    def toXML_strList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label

        xmlString = [ '%s<%s%s/>' % ( indent, self.moniker, attributeStr ) ]

        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        multiplicity = cls(element.get('label'))

        xPath.pop( )
        return( multiplicity )

def parseNodeUsingClass(cls, multElement, xPath, linkData, **kwargs):
    """Translate an xml <multiplicity> element into fudge."""

    xPath.append( multElement.tag )

    multiplicityComponent = cls()

    formClasses = {}
    for formClass in [ Unspecified, Constant1d, XYs1d, Regions1d, Reference, Polynomial1d, PartialProduction, Gridded1d, Branching1d ] :
        formClasses[formClass.moniker] = formClasses
    for form in multElement :
        formClass = formClasses.get( form.tag )
        if( formClass is None ) : raise Exception( "encountered unknown multiplicity form: %s" % form.tag )
        newForm = formClass.parseNodeUsingClass(form, xPath, linkData, **kwargs)
        multiplicityComponent.add( newForm )

    xPath.pop( )

    return( multiplicityComponent )

#
# multiplicity component
#
class Component(abstractClassesModule.Component):

    moniker = 'multiplicity'

    def __init__(self):

        abstractClassesModule.Component.__init__(
            self, (Unspecified, Constant1d, XYs1d, Regions1d, Reference, Polynomial1d,
                   PartialProduction, Gridded1d, Branching1d))

    def isConstant(self):

        return isinstance(self.evaluated, Constant1d)

    def getConstant(self):

        if self.isConstant(): return self.evaluated.evaluate(0)
        raise Exception('multiplicity type = "%s" is not single valued' % self.evaluated.moniker)

    @property
    def domainMin(self):

        return self.evaluated.domainMin

    @property
    def domainMax(self):

        return self.evaluated.domainMax

    @property
    def domainUnit(self):

        return self.evaluated.domainUnit

    def getEnergyLimits(self, EMin, EMax):

        return self.evaluated.getEnergyLimits(EMin, EMax)
        
    def evaluate(self, E):
        """
        Return the evaluated multiplicity at incident energy E.
        """

        return self.evaluated.evaluate(E)

    def isSpecified(self):
        """
        Returns **True** if *self* has data and if the evaluated form is not **Unspecified**.
        """

        if len(self) == 0:
            return False
        if isinstance(self.evaluated, Unspecified):
            return False

        return True

    def check(self, info):

        from fudge import warning
        warnings = []

        if self.isConstant():
            if self.getConstant() < 1:
                warnings.append(warning.NegativeMultiplicity(self.getConstant(), self))
        else:   # energy-dependent mult.
            for form in self:
                domain = None
                if hasattr(form, 'domain'):
                    domain = form.domainMin, form.domainMax
                if domain is not None and domain != info['crossSectionDomain']:
                    # photon multiplicities can start above cross section domainMin, but domainMax should agree:
                    if self.ancestor.pid == IDsPoPsModule.photon:
                        startRatio = domain[0] / info['crossSectionDomain'][0]
                        endRatio = domain[1] / info['crossSectionDomain'][1]
                        if (startRatio < 1 - standardsModule.Floats.epsilon or endRatio < 1 - standardsModule.Floats.epsilon
                            or endRatio > 1 + standardsModule.Floats.epsilon):
                            warnings.append(warning.Domain_mismatch(
                                *(domain + info['crossSectionDomain']), obj=self))
                    # For all other products both domainMin and domainMax should agree within epsilon
                    else:
                        for idx, (e1, e2) in enumerate(zip(domain, info['crossSectionDomain'])):
                            ratio = e1 / e2
                            if (ratio < 1 - standardsModule.Floats.epsilon or ratio > 1 + standardsModule.Floats.epsilon):
                                if ( idx==0 and domain[0] > info['crossSectionDomain'][0] and form.evaluate( domain[0] ) == 0 ):
                                    continue    # multiplicity starting after xsc is ok if first multiplicity point == 0
                                warnings.append(warning.Domain_mismatch( *(domain + info['crossSectionDomain']), obj=self))
                                break

                if hasattr(form, 'rangeMin') and form.rangeMin < 0:
                    warnings.append(warning.NegativeMultiplicity(form.rangeMin, obj=form))

        return warnings
