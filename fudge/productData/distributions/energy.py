# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
This module containes the energuy distribution classes for storing :math:`P(E')` and :math:`P(E'|E)` distributions
for uncorrelated distributions, as well as the :math:`P(E'|E)` when the
correlated distribution is represented as :math:`P(E'|E)` * :math:`P(\mu|E,E')`.
    
This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | XYs1d                             | Class to Store :math:`P(E')` as a list of :math:`Ei'_i`, :math:`P_i`. |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Regions1d                         | Stores :math:`P(E')` as a list of other 1d energy distribuitions      |
    |                                   | that abut in :math:`E'` and are contiguous.                           |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Xs_pdf_cdf1d                      | Stores :math:`P(E')` (i.e., :math:`pdf(E')`) and :math:`cdf(E')`      |
    |                                   | as a list of :math:`E'_i`, :math:`pdf_i` and :math:`cdf_i`.           |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Subform                           | Base class for all 2d energuy distributions.                          |
    +-----------------------------------+-----------------------------------------------------------------------+
    | DiscretePrimaryGamma              | Base class for :py:class:`DiscreteGamma` and :py:class:`PrimaryGamma`.|
    +-----------------------------------+-----------------------------------------------------------------------+
    | DiscreteGamma                     | Specifies that the product is a discrete (i.e., secondary) photon.    |
    +-----------------------------------+-----------------------------------------------------------------------+
    | PrimaryGamma                      | Specifies that the product is a capture primary photon.               |
    +-----------------------------------+-----------------------------------------------------------------------+
    | XYs2d                             | Stores :math:`P(E'|E)` as a list of :math:`E`, :math:`P(E')`.         |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Regions2d                         | Stores :math:`P(E'|E)` as a list of other 2d energy distributions,    |
    |                                   | that abut to form a contiguous distribution in E.                     |
    +-----------------------------------+-----------------------------------------------------------------------+
    | EnergyFunctionalData              | Abstract class for some parameters in many of the 2d energy           |
    |                                   | functionals.                                                          |
    +-----------------------------------+-----------------------------------------------------------------------+
    | A                                 | Class for storing the A function for the :py:class:`Watt` energy      |
    |                                   | functional.                                                           |
    +-----------------------------------+-----------------------------------------------------------------------+
    | B                                 | Class for storing the B function for the :py:class:`Watt` energy      |
    |                                   | functional.                                                           |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Theta                             | Class for storing the :math:`\theta` function used in several         |
    |                                   | 2d energy functionals.                                                |
    +-----------------------------------+-----------------------------------------------------------------------+
    | G                                 | Class for storing the G function for the                              |
    |                                   | :py:class:`GeneralEvaporation` functional.                            |
    +-----------------------------------+-----------------------------------------------------------------------+
    | T_M                               | Class for storing the :math:`T_M` function for the                    |
    |                                   | :py:class:`MadlandNix` functional.                                    |
    +-----------------------------------+-----------------------------------------------------------------------+
    | FunctionalBase                    | Abstract base class for the 2d energy functional distributions.       |
    +-----------------------------------+-----------------------------------------------------------------------+
    | GeneralEvaporation                | Class for storing the 2d general evaporation energy distribution.     |
    +-----------------------------------+-----------------------------------------------------------------------+
    | SimpleMaxwellianFission1d         | Class for storing the 1d simple Maxwellian fission energy             |
    |                                   | distribution.                                                         |
    +-----------------------------------+-----------------------------------------------------------------------+
    | SimpleMaxwellianFission           | Class for storing the 2d simple Maxwellian fission distribution.      |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Evaporation1d                     | Class for storing the 1d evaporation energy distribution.             |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Evaporation                       | Class for storing the 2d evaporation energy distribution.             |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Watt                              | Class for storing the 2d Watt energy distribution.                    |
    +-----------------------------------+-----------------------------------------------------------------------+
    | MadlandNix                        | Class for storing the Madland/Nix energy distribution.                |
    +-----------------------------------+-----------------------------------------------------------------------+
    | NBodyPhaseSpace                   | Class for storing the 2d N-body phase space energy distribution.      |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Weighted1d                        | Class for storing an instance of a :py:class:`Weighted1d`             |
    |                                   | at a specific projectile energy.                                      |
    +-----------------------------------+-----------------------------------------------------------------------+
    | WeightedFunctionals1d             | Class for storing an instance of a :py:class:`WeightedFunctionals`    |
    |                                   | at a specific projectile energy (this is a 1d energy functional).     |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Weighted                          | Class for storing one 2d energy distribution in the                   |
    |                                   | :py:class: WeightedFunctionals instance.                              |
    +-----------------------------------+-----------------------------------------------------------------------+
    | WeightedFunctionals               | Class for storing a list of 2d energy functional distributions as     |
    |                                   | :py:class:`Weighted` instances that are weighted to form a            |
    |                                   | complete 2d energy distribution.                                      |
    +-----------------------------------+-----------------------------------------------------------------------+
    | EnergyFunctionalDataToPointwise   | For internal use. This class is used by the methods                   |
    |                                   | toPointwise_withLinearXYs of the 2d functional energy distribution    |
    |                                   | to convert the functional into an :py:class:`XYs2d` instance.         |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Form                              | This class stores the enegy distribution part (i.e., :math:`P(E'|E)`) |
    |                                   | for an uncorrelated distribution.                                     |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

import math

from LUPY import ancestry as ancestryModule

from pqu import PQU as PQUModule
from fudge.core.utilities import brb

from xData import enums as xDataEnumsModule
from xData import base as xDataBaseModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule

from fudge import physicalQuantity as physicalQuantityModule

from . import base as baseModule
from . import probabilities as probabilitiesModule
from . import miscellaneous as miscellaneousModule

def defaultAxes( energyUnit ) :
    """
    Returns an :py:class:`Axes` instance for energy data (i.e., :math:`P(E'|E)`.

    :param energyUnit:  Unit for the energy data.

    :return:            An :py:class:`Axes` instance.
    """

    axes = axesModule.Axes(3)
    axes[2] = axesModule.Axis( 'energy_in',  2, energyUnit )
    axes[1] = axesModule.Axis( 'energy_out', 1, energyUnit )
    axes[0] = axesModule.Axis( 'P(energy_out|energy_in)', 0, '1/' + energyUnit )
    return( axes )

class XYs1d( XYs1dModule.XYs1d ) :
    """
    Class to Store :math:`P(E')` as a list of pairs of :math:`Ei'_i`, :math:`P_i`.
    """

    def averageEnergy( self ) :
        """
        The method calculates the average outgoing energy of *self*.
        """

        allowedInterpolations = [xDataEnumsModule.Interpolation.linlin, xDataEnumsModule.Interpolation.flat]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations, XYs1dModule.defaultAccuracy )
        return( xys.integrateWithWeight_x( ) )

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )

class Regions1d( regionsModule.Regions1d ) :
    """
    Stores :math:`P(E')` as a list of other 1d energy distribuitions that abut in :math:`E'` and are contiguous. 
    """

    def averageEnergy( self ) :
        """
        The method calculates the average outgoing energy of *self*.
        """

        averageEnergy = 0
        for region in self : averageEnergy += region.averageEnergy( )
        return( averageEnergy )

    def integrateWithWeight_x( self ) :
        r"""
        Returns the value of the integral :math:`\int dE' E' P(E').

        :return:        A float.
        """

        sum = 0
        for region in self : sum += region.integrateWithWeight_x( )
        return( sum )

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )

    @staticmethod
    def allowedSubElements():
        """
        This method returns the list of classes that can be sub-nodes (i.e., 2-d function) of an :py:class:`Regions1d` instance.
        """

        return( XYs1d, )

class Xs_pdf_cdf1d( xs_pdf_cdfModule.Xs_pdf_cdf1d ) :
    """
    This class stores both the :math:`P(E'|E)` (i.e., pdf) and its cumlative distribution function cdf as three lists: one list
    for the :math:`E'` values, and one list for the associated pdf value and one list for the associated cdf value for each :math:`E'` value,
    """

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )

class Subform( baseModule.Subform ) :
    """Abstract base class for 2d energy subforms."""

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns None.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                None.
        """


        return( None )

class DiscretePrimaryGamma( Subform ) :
    """
    Abstract base class for classes :py:class:`DiscreteGamma` and :py:class:`PrimaryGamma`.

    :param value:       Discrete photon energy or nuclear binding energy needed to calcualte the ergy for a primary photon.
    :param domainMin:   Minimum projectile energy for which the photon is emitted.
    :param domainMax:   Maximum projectile energy for which the photon is emitted.
    :param axes:        Axes instance for the data.
    """

    dimension = 2

    def __init__( self, value, domainMin, domainMax, axes = None ) :

        Subform.__init__( self )

        if( isinstance( value, int ) ) : value = float( value )
        if( not( isinstance( value, float ) ) ) : raise TypeError( 'value must be a float.' )
        self.value = value

        if( isinstance( domainMin, int ) ) : domainMin = float( domainMin )
        if( not( isinstance( domainMin, float ) ) ) : raise TypeError( 'domainMin must be a float.' )
        self.__domainMin = domainMin

        if( isinstance( domainMax, int ) ) : domainMax = float( domainMax )
        if( not( isinstance( domainMax, float ) ) ) : raise TypeError( 'domainMax must be a float.' )
        self.__domainMax = domainMax

        if( axes is None ) :
            self.__axes = None
        else :
            if( not( isinstance( axes, axesModule.Axes ) ) ) : raise TypeError( 'axes is not an axes instance' )
            if( len( axes ) <= self.dimension ) : raise Exception( 'len( axes ) = %d != ( self.dimension + 1 ) = %d' % ( len( axes ), ( self.dimension + 1 ) ) )
            self.__axes = axes.copy( )
            self.__axes.setAncestor( self )

    @property
    def axes( self ) :
        """Returns a return to the axes of *self*."""

        return( self.__axes )

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy for *self*.
        """

        return( self.__domainMin )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy for *self*.
        """

        return( self.__domainMax )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile.
        """

        return( self.__axes[-1].unit )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        factors = self.axes.convertUnits( unitMap )
        self.value *= factors[1]
        self.__domainMin *= factors[2]
        self.__domainMax *= factors[2]

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        return self.__class__( self.value, self.__domainMin, self.__domainMax, self.axes )

    __copy__ = copy

    def energySpectrumAtEnergy( self, energy, discreteGammaResolution = 1e-2 ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energy*,

        :param energyIn:                    Energy of the projectile.
        :param discreteGammaResolution:     Relative width of triangle presenting the energy spectrum.

        :return:                    XYs1d instance for the energy spectrum.
        """

        if( ( self.__domainMin > energy ) or ( self.__domainMax < energy ) ) : return( XYs1d( axes = defaultAxes( self.domainUnit ) ) )

        photonEnergy = self.energyAtEnergy( energy )
        energy1 = photonEnergy * ( 1.0 - discreteGammaResolution )
        energy2 = photonEnergy * ( 1.0 + discreteGammaResolution )
        height = 2.0 / ( energy2 - energy1 )
        return( XYs1d( data = [ [ energy1, 0.0 ], [ photonEnergy, height ], [ energy2, 0.0 ] ], axes = defaultAxes( self.domainUnit ) ) )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Sets *domainMin* and *domainMax* per the arguments.
        """

        OldDomainMin = self.domainMin
        OldDomainMax = self.domainMax

        domainMin = max(domainMin, self.domainMin)
        domainMax = min(domainMax, self.domainMax)
        if fixToDomain == xDataEnumsModule.FixDomain.lower:
            if domainMin > self.__domainMin: self.__domainMin = domainMin
        elif fixToDomain == xDataEnumsModule.FixDomain.upper:
            self.__domainMax = domainMax
        else:
            self.__domainMin = domainMin
            self.__domainMax = domainMax

        if OldDomainMin == self.domainMin and OldDomainMax == self.domainMax: return 0
        return 1

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        extraAttribute = ''
        if isinstance(self, PrimaryGamma):
            if self.finalState is not None:
                extraAttribute = ' finalState="%s"' % self.finalState

        XMLStringList = ['%s<%s value="%s" domainMin="%s" domainMax="%s"%s' % (indent, self.moniker, PQUModule.floatToShortestString(self.value, 12), 
                PQUModule.floatToShortestString(self.__domainMin, 12), PQUModule.floatToShortestString(self.__domainMax, 12), extraAttribute)]

        if( self.axes is None ) : 
            XMLStringList[-1] += '/>'
        else :
            XMLStringList[-1] += '>'
            XMLStringList += self.axes.toXML_strList( indent=indent2, **kwargs )
            XMLStringList[-1] += '</%s>' % self.moniker

        return XMLStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        value = float(node.get('value'))
        domainMin = float(node.get('domainMin'))
        domainMax = float(node.get('domainMax'))
        finalState = node.get('finalState', None)

        axes = None
        for child in node:
            if( child.tag == axesModule.Axes.moniker ) :
                axes = axesModule.Axes.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise Exception('Invalid child node with tag = "%s"' % child.tag)

        if issubclass(cls, PrimaryGamma):
            return cls(value, domainMin, domainMax, axes=axes, finalState=finalState)
        else:
            return cls(value, domainMin, domainMax, axes=axes)

    def getEnergyArray(self, EMin=None, EMax=None):
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        return [self.domainMin, self.domainMax]

class DiscreteGamma(DiscretePrimaryGamma):
    """
    This class stores the data for a discrete (i.e., secondary) gamma that results from the decay of a nucleus
    from one nuclear energy level to a lower nuclear energy level,
    """

    moniker = 'discreteGamma'

    def check(self, info):
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        if self.value <= 0: warnings.append(warning.NegativeDiscreteGammaEnergy())
        return warnings

    def averageEp(self, E):
        """
        Calculated the average energy of the outgoing photon which is just its energy.

        :param E:       The energy of the projectile energy. This is not needed for a discrete photon.

        :return:                Python float.
        """

        return self.energyAtEnergy(E)

    def energyAtEnergy(self, energyIn):
        """
        Returns the energy of the outgoing photon which is just its energy.

        :param energyIn:        The energy of the projectile energy. This is not needed for a discrete photon.

        :return:                Python float.
        """

        if energyIn < self.domainMin:
            return 0
        return self.value

    def integrate(self, energyIn, energyOut):
        """
        This meethod integrates the energu distribution at projectile energy over the specified outgoing energy range.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*.
        
        :param energyIn:            The energy of the projectile.
        :param energyOut:           The outgoing energy range to integrate over.

        :return:                    A float representing the value of the integration.
        """

        if self.domainMin <= energyIn <= self.domainMax:
            domainMin, domainMax = miscellaneousModule.domainLimits(energyOut, self.value, self.value)
            if domainMin <= self.value <= domainMax: return 1.0

        return 0.0


class PrimaryGamma(DiscretePrimaryGamma):

    moniker = 'primaryGamma'

    def __init__(self, value, domainMin, domainMax, axes=None, finalState=None):
        """Constructor.

        :param value:           The neutron binding energy.
        :param domainMin:       The minimum projectile energy that the primary photon is emitted.
        :param domainMax:       The maximum projectile energy that the primary photon is emitted.
        :param axes:            The axes for *self*.
        :param finalState:      The nuclear state the compound is in after the primary photon (gamma) is emitted.
        """

        DiscretePrimaryGamma.__init__(self, value, domainMin, domainMax, axes=axes)
        self.__massRatio = None     # In ENDF lingo this is AWR / ( AWR + 1 ).

        self.finalState = finalState

    @property
    def finalState(self):
        """Returns the finalState."""

        return self.__finalState

    @finalState.setter
    def finalState(self, value):
        """
        Set the final state member to *value*.

        :param value:   The value to set finalState to.
        """

        if not isinstance(value, (str, type(None))):
            raise TypeError('Invalid finalState type: got "%s".' % type(value))

        self.__finalState = value

    @property
    def massRatio(self):
        """
        Returns the ratio of the target mass to the sum of the projectile and target masses.

        :return:        Python float.
        """

        if self.__massRatio is None:
            self.__massRatio = self.findAttributeInAncestry("getMassRatio")()
        return self.__massRatio

    def check(self, info):
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        Qvalue = self.findAttributeInAncestry('getQ')('eV')
        if isinstance(self.value, PQUModule.PQU):
            testValue = self.value.getValueAs('eV')
        else:
            testValue = self.value
        if testValue > Qvalue:
            warnings.append(warning.PrimaryGammaEnergyTooLarge(self.value, 100 * testValue / Qvalue))
        return warnings

    def averageEp(self, energyIn):
        """
        Calculated the average energy of the outgoing photon at projectile energy *energyIn*.

        :param energyIn:       The energy of the projectile energy.

        :return:                Python float.
        """

        return self.energyAtEnergy(energyIn)

    def energyAtEnergy(self, energyIn):
        """
        Returns the energy of the outgoing photon for projectile energy *energyIn*.

        :param energyIn:        The energy of the projectile energy.

        :return:                Python float.
        """

        if energyIn < self.domainMin:
            return 0
        return float(self.value) + self.massRatio * energyIn

    def integrate(self, energyIn, energyOut):
        """
        This meethod integrates the energu distribution at projectile energy over the specified outgoing energy range.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*.
        
        :param energyIn:            The energy of the projectile.
        :param energyOut:           The outgoing energy range to integrate over.

        :return:                    Python float.
        """

        gammaEnergy = float(self.value) + self.massRatio * energyIn
        if self.domainMin <= energyIn <= self.domainMax:
            domainMin, domainMax = miscellaneousModule.domainLimits(energyOut, gammaEnergy, gammaEnergy)
            if domainMin <= gammaEnergy <= domainMax: return 1.0

        return 0.0


class XYs2d( Subform, probabilitiesModule.PofX1GivenX2 ) :
    """
    Class to store :math:`P(E'|E)` as a list of pairs of :math:`E`, :math:`P(E')`. The :math:`P(E')` can be any of the
    1d energy distributions.

    :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
    """

    def __init__( self, **kwargs ):

        probabilitiesModule.PofX1GivenX2.__init__( self, **kwargs )
        Subform.__init__( self )

    def evaluate(self, domainValue, extrapolation=xDataEnumsModule.Extrapolation.none, epsilon = 0, **kwargs):
        """
        Returns a :math:`P(E')` at projectile energy *domainValue*.

        :param domainValue:     The energy of the projectile.
        :param extrapolation:   Specifies how to extrapolate *self* for *domainValue* outside the domain of the projectile's energy.

        :return:                A 1d energy distribution.
        """

        # FIXME silently converting InterpolationQualifier.none to unitbase.
        interpolationQualifier = self.interpolationQualifier
        if interpolationQualifier is xDataEnumsModule.InterpolationQualifier.none:
            interpolationQualifier = xDataEnumsModule.InterpolationQualifier.unitBase

        return probabilitiesModule.PofX1GivenX2.evaluate(self, domainValue, extrapolation = extrapolation, epsilon = epsilon,
                                                         interpolationQualifier=interpolationQualifier, **kwargs)


    def getAtEnergy( self, energy ) :
        """This method is deprecated, use getSpectrumAtEnergy."""

        return( self.getSpectrumAtEnergy( energy ) )

    def energySpectrumAtEnergy( self, energy ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energy*,

        :param energy:              Energy of the projectile.

        :return:                    XYs1d instance for the energy spectrum.
        """

        spectrum = self.evaluate(energy, extrapolation=xDataEnumsModule.Extrapolation.flat)
        if( isinstance( spectrum, Regions1d ) ) : spectrum = spectrum.toPointwise_withLinearXYs( lowerEps = 1e-6, upperEps = 1e-6 )
        return( spectrum )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the results of calling **energySpectrumAtEnergy**. This method is deprecated."""

        return( self.energySpectrumAtEnergy( energy ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        Es = [ data.outerDomainValue for data in self ]
        if( EMin is not None ) :
            if( EMin < ( 1.0 - 1e-15 ) * Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def averageEp( self, energy ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *energy*.

        :param energy:       The energy of the projectile energy.

        :return:                Python float.
        """

        code, lower, upper, frac, interpolation, interpolationQualifier = self.getBoundingSubFunctions( energy )
        if( code is None ) :
            raise Exception( 'No distribution' )
        elif( code == '' ) :
            EpLowerMin, EpLowerMax = lower.domainMin, lower.domainMax
            EpUpperMin, EpUpperMax = upper.domainMin, upper.domainMax
            EpMidMin = ( 1 - frac ) * EpLowerMin + frac * EpUpperMin
            EpMidMax = ( 1 - frac ) * EpLowerMax + frac * EpUpperMax
            EpLower = ( lower.integrateWithWeight_x( ) - EpLowerMin ) / ( EpLowerMax - EpLowerMin ) * \
                    ( EpMidMax - EpMidMin ) + EpMidMin
            EpUpper = ( upper.integrateWithWeight_x( ) - EpUpperMin ) / ( EpUpperMax - EpUpperMin ) * \
                    ( EpMidMax - EpMidMin ) + EpMidMin
            return( ( 1 - frac ) * EpLower + frac * EpUpper )
        else :
            return( lower.integrateWithWeight_x( ) )

    def check(self, info):
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []

        if self.interpolation is xDataEnumsModule.Interpolation.flat:
            warnings.append(warning.FlatIncidentEnergyInterpolation())

        if self.interpolationQualifier is xDataEnumsModule.InterpolationQualifier.none:
            domains = set([(xys1d.domainMin, xys1d.domainMax) for xys1d in self])
            if len(domains) != 1:
                transportable = self.findAttributeInAncestry('pid') in info['transportables']
                warnings.append(warning.MissingInterpolationQualifier(transportable))

        for idx in range(len(self)):
            integral = self[idx].integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                incidentEnergy = PQUModule.PQU(self[idx].outerDomainValue, self.axes[-1].unit)
                warnings.append(warning.UnnormalizedDistribution(incidentEnergy, idx, integral, self[idx]))

            if self[idx].rangeMin < 0.0:
                incidentEnergy = PQUModule.PQU(self[idx].outerDomainValue, self.axes[-1].unit)
                warnings.append(warning.NegativeProbability(self[idx].rangeMin, incidentEnergy, obj=self[idx]))

        return warnings

    def sqrtEp_AverageAtE( self, E ) :
        r"""
        This method returns the value of the integral :math:`\int dE' \sqrt{E'} P(E') at projectile energy *E*.

        :param E:   Energy of the projectile.

        :return:    Python float.
        """

        return( self.energySpectrumAtEnergy( E ).integrateWithWeight_sqrt_x( ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        return( multiD_XYsModule.XYs2d.toPointwise_withLinearXYs( self, cls = XYs2d, **kwargs ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* as a :py:class:`XYs2d` representation with the
        P(E') data as :py:class:`Xs_pdf_cdf1d` instances.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
                    
        :return:                :py:class:`XYs2d`
        """     

        linear = self
        for xys in self :
            if( isinstance( xys, XYs1d ) ) :
                if xys.interpolation not in [xDataEnumsModule.Interpolation.linlin, xDataEnumsModule.Interpolation.flat]:
                    linear = self.toPointwise_withLinearXYs( accuracy = XYs1dModule.defaultAccuracy, upperEps = 1e-8 )
                    break
            else :
                linear = self.toPointwise_withLinearXYs( accuracy = XYs1dModule.defaultAccuracy, upperEps = 1e-8 )
                break
        subform = XYs2d( axes = self.axes, interpolation = self.interpolation, 
                interpolationQualifier = self.interpolationQualifier )
        for xys in linear:
            subform.append(Xs_pdf_cdf1d.fromXYs(xys, xys.outerDomainValue, thinEpsilon=1e-14))

        return subform

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 2-d function) of an :py:class:`XYs2d` instance.
        """

        return( ( XYs1d, Regions1d, Xs_pdf_cdf1d ) )

class Regions2d( Subform, regionsModule.Regions2d ) :
    """
    Stores :math:`P(E'|E)` as a list of other 2d energy distributions, that abut to form a contiguous distribution in E. 

    :param kwargs:              A dictionary that contains data to control the way this method acts.
    """

    def __init__( self, **kwargs ):

        regionsModule.Regions2d.__init__( self, **kwargs )
        Subform.__init__( self )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        for idx, region in enumerate( self ):
            regionWarnings = region.check( info )
            if regionWarnings:
                warnings.append( warning.Context("Region %d:" % idx, regionWarnings) )
        return warnings

    def energySpectrumAtEnergy(self, energyIn):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn*,

        :param energyIn:                    Energy of the projectile.
        :param discreteGammaResolution:     Relative width of triangle presenting the energy spectrum.

        :return:                    XYs1d instance for the energy spectrum.
        """

        for region in self:
            if energyIn < region.domainMax:
                return region.energySpectrumAtEnergy(energyIn)

        return XYs1d([], axes=defaultAxes(energyUnit=self.domainUnit))

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        return( regionsModule.Regions2d.toPointwise_withLinearXYs( self, cls = XYs2d, **kwargs ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* as a :py:class:`Regions2d` representation with the
        P(E') data as :py:class:`Xs_pdf_cdf1d` instances.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
                    
        :return:                :py:class:`Regions2d`
        """     

        _regions2d = Regions2d( axes = self.axes )
        for region in self : _regions2d.append( region.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( _regions2d )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 2-d function) of an :py:class:`Regions2d` instance.
        """

        return( ( XYs2d, ) )

class EnergyFunctionalData( ancestryModule.AncestryIO ) :
    """
    Abstract class for some parameters in many of the 2d energy functionals.

    The following table list the primary members of this class:
    
    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | data      | This member stores the actual data.                       |           
    +-----------+-----------------------------------------------------------+

    :param data:    The data for the function.
    """

    ancestryMembers = ( 'data', )

    def __init__( self, data ) :

        ancestryModule.AncestryIO.__init__( self )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.data.convertUnits( unitMap )

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        return self.__class__( self.data.copy( ) )

    __copy__ = copy

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXML_strList( indent + '  ', **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )
        subClass = {
            XYs1d.moniker           : XYs1dModule.XYs1d,
            Regions1d.moniker       : regionsModule.Regions1d,
            Xs_pdf_cdf1d.moniker    : xs_pdf_cdfModule.Xs_pdf_cdf1d
        }.get( element[0].tag )
        if( subClass is None ) : raise Exception( "encountered unknown energy functional subform: %s" % element[0].tag )
        EFD = cls(subClass.parseNodeUsingClass(element[0], xPath, linkData, **kwargs))
        xPath.pop()
        return EFD

class A( EnergyFunctionalData ) :
    """
    Class for storing the A function for the :py:class:`Watt` energy functional. 
    """

    moniker = 'a'

class B( EnergyFunctionalData ) :
    """
    Class for storing the B function for the :py:class:`Watt` energy functional. 
    """

    moniker = 'b'

class Theta( EnergyFunctionalData ) :
    """
    Class for storing the :math:`\theta` function for several of the energy functionals.
    """

    moniker = 'theta'

class G( EnergyFunctionalData ) :
    """
    Class for storing the G function for the :py:class:`GeneralEvaporation` functional.  
    """

    moniker = 'g'

class T_M( EnergyFunctionalData ) :
    """
    Class for storing the :math:`T_M` function for the :py:class:`MadlandNix` functional.
    """

    moniker = 'T_M'

class FunctionalBase( Subform ) :
    """
    Abstract base class for the 2d energy functional distributions.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | LF            | This member is the ENDF-6 LF flag.                                |
    +---------------+-------------------------------------------------------------------+
    | U             | Constant that defines the upper energy limit for the              |
    |               | secondary particle so that 0 ≤ E' ≤ E - U.                        |
    +---------------+-------------------------------------------------------------------+
    | parameter1    | 1d function needed by most functionals.                           |
    +---------------+-------------------------------------------------------------------+
    | parameter2    | Addition 1d function needed by :py:class:`GeneralEvaporation`     |
    |               | and :py:class:`MadlandNix` instances.                             |
    +---------------+-------------------------------------------------------------------+

    :param LF:              The data for the function.
    :param U:               Constant that defines the upper energy limit for the secondary particle so that 0 ≤ E' ≤ E - U.
    :param parameter1:      For a :py:class:`GeneralEvaporation`, :py:class:`SimpleMaxwellianFission` or :py:class:`Evaporation` 
                            instance this is the :math:`\theta(E)` function. For a :py:class:`Watt` instance this is the :py:class:`A`
                            function and for a :py:class:`MadlandNix` instance this is the :py:class:`T_M` function.
    :param parameter2:      For a :py:class:`GeneralEvaporation` instance this is the :py:class:`B` function and for a :py:class:`Watt` 
                            instance this is the :py:class:`B` function.
    """

    ancestryMembers = ( 'parameter1', 'parameter2' )

    def __init__( self, LF, U, parameter1, parameter2 = None ) :

        Subform.__init__( self )

        if( U is not None ) :
            if( not( isinstance( U, physicalQuantityModule.U ) ) ) : raise TypeError( 'Invalid U type' )
        self.U = U
        self.LF = LF

        self.parameter1 = parameter1
        self.parameter1.setAncestor( self )

        self.parameter2 = parameter2
        if( parameter2 is not None ) : self.parameter2.setAncestor( self )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if( self.U is not None ) : self.U.convertUnits( unitMap )
        self.parameter1.convertUnits( unitMap )
        if( self.parameter2 is not None ) : self.parameter2.convertUnits( unitMap )

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        U = self.U
        if U is not None: U = self.U.copy()
        if self.parameter2 is None:
            return self.__class__(U, self.parameter1.copy())
        else :
            return self.__class__(U, self.parameter1.copy(), self.parameter2.copy())

    __copy__ = copy

    def check( self, info ):
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        if( ( self.domainMin - self.U.value ) < 0 ) :
            warnings.append( warning.EnergyDistributionBadU( self ) )
        return( warnings )

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy for *self*.
        """

        return( self.parameter1.data.domainMin )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy for *self*.
        """

        return( self.parameter1.data.domainMax )

    @property
    def domainUnit(self):
        """
        Returns the energy unit of the projectile.
        """

        return self.parameter1.data.axes[0].unit

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        if( isinstance( self.parameter1.data, regionsModule.Regions1d ) ) :
            Es = []
            for region in self.parameter1.data :
                Es = Es[:-1] + [ E for E, p in region ]
        else :
            Es = [ E for E, p in self.parameter1.data ]
        if( EMin is not None ) :
            if( EMin < ( 1.0 - 1e-15 ) * Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        Calls the *fixDomains* for the **parameter1** member and if *self* is a **Watt** spectrum the **parameter2** member.
        """

        numberOfFixes = self.parameter1.data.fixDomains(energyMin, energyMax, domainToFix)
        if isinstance(self, Watt): numberOfFixes += self.parameter2.data.fixDomains(energyMin, energyMax, domainToFix)

        return numberOfFixes

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        if( self.LF == 12 ) : 
            xmlString += self.EFL.toXML_strList( indent2, **kwargs )
            xmlString += self.EFH.toXML_strList( indent2, **kwargs )
        else :
            xmlString += self.U.toXML_strList( indent2, **kwargs )
        xmlString += self.parameter1.toXML_strList( indent2, **kwargs )
        if self.parameter2 is not None: xmlString += self.parameter2.toXML_strList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class GeneralEvaporation( FunctionalBase ) :
    r"""
    Class for storing the 2d general evaporation energy distribution. The distribution is of the form
    :math:`g(E' / theta(E))`.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | U             | Constant that defines the upper energy limit for the              |
    |               | secondary particle so that 0 ≤ E' ≤ E - U.                        |
    +---------------+-------------------------------------------------------------------+
    | thetas        | The :math:`theta(E)` function.                                    |
    +---------------+-------------------------------------------------------------------+
    | gs            | The :math:`g(x)` function.                                        |
    +---------------+-------------------------------------------------------------------+

    :param U:               Constant that defines the upper energy limit for the secondary particle so that 0 ≤ E' ≤ E - U.
                            This is not used/needed so should be removed.
    :param thetas:          The :math:`theta(E)` function.
    :param gs:              The :math:`g(x)` function.
    """

    moniker = 'generalEvaporation'

    def __init__( self, U, thetas, gs ) :

        FunctionalBase.__init__( self, 5, U, thetas, gs )

    def convertUnits( self, unitMap ) :
        """
        Overriding method from base class due to a processing error
        (axes for 'g' parameter not copied during processing).
        """
        #FIXME: This method should be removed once processed libraries are replaced.

        if( self.U is not None ) : self.U.convertUnits( unitMap )
        self.parameter1.convertUnits( unitMap )
        # skip 'g': it's unitless and may not have axes defined
        #if( self.parameter2 is not None ) : self.parameter2.convertUnits( unitMap )

    def averageEp( self, E ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *E*.

        :param E:       The energy of the projectile energy.

        :return:                Python float.
        """

        return( self.parameter1.data.evaluate( E ) * self.parameter2.data.integrateWithWeight_x( ) )

    def sqrtEp_AverageAtE( self, E ) :
        r"""
        This method returns the value of the integral :math:`\int dE' \sqrt{E'} P(E') at projectile energy *E*.

        :param E:   Energy of the projectile.

        :return:    Python float.
        """

        return( math.sqrt( self.parameter1.data.evaluate( E ) ) * self.parameter2.data ).integrateWithWeight_sqrt_x( )

    def energySpectrumAtEnergy( self, energyIn ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn*,

        :param energyIn:            Energy of the projectile.

        :return:                    XYs1d instance for the energy spectrum.
        """

        theta = self.parameter1.data.evaluate(energyIn)
        spectrum = self.parameter2.data.toPointwise_withLinearXYs( accuracy = 1e-5, lowerEps = 1e-6, upperEps = 1e-6 )
        spectrum = spectrum.nf_pointwiseXY.scaleOffsetXAndY(xScale = theta, yScale = 1.0 / theta )
        return( XYs1d( spectrum, axes = defaultAxes( energyUnit = self.domainUnit ) ) )

    def isLinear( self, qualifierOk = False, flatIsOk = False ) :
        r"""
        Returns the results of isLinear called on the axes of :math:`g(E'|E)`.
        """

        return( self.parameter2.axes.isLinear( qualifierOk = qualifierOk, flatIsOk = flatIsOk ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* with its g data as a :py:class:`Xs_pdf_cdf1d` instance.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                :py:class:`GeneralEvaporation`
        """

        _gs = G(Xs_pdf_cdf1d.fromXYs(self.parameter2.data, thinEpsilon=1e-14))
        _gs.data.axes = self.parameter2.data.axes.copy()
        _form = GeneralEvaporation( self.U, thetas = self.parameter1.copy( ), gs = _gs )
        return( _form )

    def toPointwise_withLinearXYs( self, **kwargs  ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        pwl = XYs2d( axes = defaultAxes( self.parameter1.data.domainUnit ) )
        thetas = self.parameter1.data.toPointwise_withLinearXYs( **kwargs )
        gs = self.parameter2.data.toPointwise_withLinearXYs( **kwargs )
        for E_in, theta in thetas :
            data = [ [ theta * x, y / theta ] for x, y in gs ]
            data = XYs1d( data, outerDomainValue = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )
        theta_ = Theta.parseNodeUsingClass(element.find(Theta.moniker), xPath, linkData, **kwargs)
        g_ = G.parseNodeUsingClass(element.find(G.moniker), xPath, linkData, **kwargs)
        U = physicalQuantityModule.U.parseNodeUsingClass(element.find( 'U' ), xPath, linkData, **kwargs)
        GES = cls( U, theta_, g_ )
        xPath.pop()
        return GES

class SimpleMaxwellianFission1d:       # FIXME, needs units
    r"""
    This class represents a :py:class:`SimpleMaxwellianFission` instance evaluated at the specified projectile energy.
    The 1d function is :math:`P(E') = \sqrt{E'} \, \exp(-E' /  \theta) / I` with :math:`\theta = \theta(E)` where E is the projectile energy
    for this function.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | energy        | The projectile energy this function is evaluated at.              |
    +---------------+-------------------------------------------------------------------+
    | U             | Constant that defines the upper energy limit for the              |
    |               | secondary particle so that 0 ≤ E' ≤ E - U.                        |
    +---------------+-------------------------------------------------------------------+
    | theta         | The value of :math:`theta(E)` at E = *energy*.                    |
    +---------------+-------------------------------------------------------------------+

    :param energy:          The projectile energy which :math:`theta = theta(E)` is evaluated at.
    :param U:               Constant that defines the upper energy limit for the secondary particle so that 0 ≤ E' ≤ E - U.
                            This is not used/needed so should be removed.
    :param theta:           The value of :math:`theta = theta(E)` evaluated for E = energy.
    """

    def __init__( self, energy, theta, U ) :

        self.energy = energy
        self.theta = theta
        self.U = U

        self.norm = self.evaluateIndefiniteIntegral( energy - U )

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy for *self*.
        """

        return( 0.0 )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy for *self*.
        """

        return( self.U )

    def evaluate( self, energy ) :
        """
        Returns the value of *self* evaluated at outgoing particle energy *energy*.

        :param energy:      The energy of the outgoing particle.

        :return:            Python float.
        """

        return( math.sqrt( energy ) * math.exp( -energy / self.theta ) / self.norm )

    def evaluateIndefiniteIntegral( self, energy ) :
        """
        Returns the value of the integral of *self* from E' = 0 to E' = *energy*.

        :param energy:      The energy of the outgoing particle.

        :return:            Python float.
        """

        X = energy / self.theta
        sqrtX = math.sqrt( X )
        return( self.theta**1.5 * ( 0.5 * math.sqrt( math.pi ) * math.erf( sqrtX ) - sqrtX * math.exp( -X ) ) )

    def integrate( self, energyMin, energyMax ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param energyMin:           The minimum value for the outgoing energy range.
        :param energyMax:           The maximum value for the outgoing energy range.

        :return:                    A float representing the value of the integration.
        """

        energyMin = max( 0.0, min( energyMin, self.U ) ) / self.theta
        energyMax = max( 0.0, min( energyMax, self.U ) ) / self.theta

        return( ( self.evaluateIndefiniteIntegral( energyMax ) - self.evaluateIndefiniteIntegral( energyMin ) ) / self.norm )

class SimpleMaxwellianFission( FunctionalBase ) :
    r"""
    Class for storing the 2d simple Maxwellian fission distribution. The distribution is defined as
    :math:`P(E') = \sqrt{E'} \, \exp(-E' /  \theta(E)) / I`.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | U             | Constant that defines the upper energy limit for the              |
    |               | secondary particle so that 0 ≤ E' ≤ E - U.                        |
    +---------------+-------------------------------------------------------------------+
    | thetas        | The :math:`theta(E)` function.                                    |
    +---------------+-------------------------------------------------------------------+

    :param U:               Constant that defines the upper energy limit for the secondary particle so that 0 ≤ E' ≤ E - U.
                            This is not used/needed so should be removed.
    :param thetas:          The :math:`theta(E)` function.
    """

    moniker = 'simpleMaxwellianFission'

    def __init__( self, U, thetas ) :

        FunctionalBase.__init__( self, 7, U, thetas )

    @property
    def theta( self ) :
        """Returns a reference to the 1d :math:`\theta(E)` function."""

        return( self.parameter1 )

    def averageEp( self, E ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *E*.

        :param E:       The energy of the projectile energy.

        :return:                Python float.
        """

        theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 1575. - a * ( 180. + 8 * a ) ) / 2625. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        erf_sqrt_a = math.sqrt( math.pi ) * math.erf( sqrt_a )
        return( theta * ( 0.75 * erf_sqrt_a - sqrt_a * ( 1.5 + a ) * exp_a ) / ( 0.5 * erf_sqrt_a - sqrt_a * exp_a ) )

    def sqrtEp_AverageAtE( self, E ) :
        r"""
        This method returns the value of the integral :math:`\int dE' \sqrt{E'} P(E') at projectile energy *E*.

        :param E:   Energy of the projectile.

        :return:    Python float.
        """

        theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( math.sqrt( theta * a ) * ( 2100. - a * ( 140. + 9. * a ) ) / 2800. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        erf_sqrt_a = math.sqrt( math.pi ) * math.erf( sqrt_a )
        return( math.sqrt( theta ) * ( 1 - ( 1. + a ) * exp_a ) / ( 0.5 * erf_sqrt_a - sqrt_a * exp_a ) )

    def energySpectrumAtEnergy( self, energyIn ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn*,

        :param energyIn:            Energy of the projectile.

        :return:                    XYs1d instance for the energy spectrum.
        """

        def A( energyOut, self ) :

            return( self.evaluate( energyOut ) )

        spectrum1 = self.evaluate( energyIn )
        spectrum2 = XYs1dModule.pointwiseXY_C.createFromFunction( [ 0.0, energyIn - self.U.value ], A, spectrum1, 1e-3, 12 )
        return( XYs1d( spectrum2, axes = defaultAxes( energyUnit = self.domainUnit ) ) )

    def evaluate( self, energy ) :
        """
        Returns a 1d function of *self* evaluated at projectile energy *energy*.

        :param energy:      The energy of the outgoing particle.

        :return:            :py:class:`SimpleMaxwellianFission1d` instance.
        """

        return( SimpleMaxwellianFission1d( energy, self.theta.data.evaluate( energy ), self.U.value ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        def evaluateAtX( self, x ) :
            """
            A function needed by :py:class:`EnergyFunctionalDataToPointwise` to evaluate a simple Maxwellian fission
            distribution at the outgoing particle energy *x* for :math:`\theta(E)` = self.p1.

            :param x:       The energy of the outgoing particle.

            :return:        Python float.
            """

            return( math.sqrt( x ) * math.exp( -x / self.p1 ) )

        ef = EnergyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( **kwargs ) )

    @classmethod
    def parseNodeUsingClass(cls, MFelement, xPath, linkData, **kwargs):
        """
        Parse *MFelement* into an instance of *cls*.

        :param cls:         Form class to return.
        :param MFelement:   Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *MFelement*.
        """

        xPath.append( MFelement.tag )
        theta_ = Theta.parseNodeUsingClass(MFelement.find(Theta.moniker), xPath, linkData, **kwargs)
        U = physicalQuantityModule.U.parseNodeUsingClass(MFelement.find( 'U' ), xPath, linkData, **kwargs)
        SMF = cls( U, theta_ )
        xPath.pop()
        return SMF

class Evaporation1d:       # FIXME, needs units
    r"""
    This class represents an :py:class:`Evaporation` instance evaluated at the specified projectile energy.
    The 1d function is :math:`P(E') = E' \, \exp(-E' /  \theta) / I` with :math:`\theta = \theta(E)` where E is the projectile energy
    for this function.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | energy        | The projectile energy this function is evaluated at.              |
    +---------------+-------------------------------------------------------------------+
    | theta         | The value of :math:`theta(E)` at E = *energy*.                    |
    +---------------+-------------------------------------------------------------------+
    | U             | Constant that defines the upper energy limit for the              |
    |               | secondary particle so that 0 ≤ E' ≤ E - U.                        |
    +---------------+-------------------------------------------------------------------+

    :param energy:          The projectile energy which :math:`theta = theta(E)` is evaluated at.
    :param theta:           The value of :math:`theta = theta(E)` evaluated for E = energy.
    :param U:               Constant that defines the upper energy limit for the secondary particle so that 0 ≤ E' ≤ E - U.
                            This is not used/needed so should be removed.
    """

    def __init__( self, energy, theta, U ) :

        self.energy = energy
        self.theta = theta
        self.U = U

        x = ( energy - U ) / theta
        self.norm = theta**2 * ( 1.0 - math.exp( -x ) * ( 1 + x ) )

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy for *self*.
        """

        return( 0.0 )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy for *self*.
        """

        return( self.U )

    def evaluate( self, energy ) :
        """
        Returns the value of *self* evaluated at outgoing particle energy *energy*.

        :param energy:      The energy of the outgoing particle.

        :return:            Python float.
        """

        return( energy * math.exp( -energy / self.theta ) / self.norm )

    def integrate( self, energyMin, energyMax ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param energyMin:           The minimum value for the outgoing energy range.
        :param energyMax:           The maximum value for the outgoing energy range.

        :return:                    A float representing the value of the integration.
        """

        energyMin = max( 0.0, min( energyMin, self.U ) ) / self.theta
        energyMax = max( 0.0, min( energyMax, self.U ) ) / self.theta

        value1 = ( 1 + energyMin ) * math.exp( -energyMin )
        value2 = ( 1 + energyMax ) * math.exp( -energyMax )
        integral = self.theta**2 * ( value1 - value2 ) / self.norm

        return( integral )

class Evaporation( FunctionalBase ) :
    r"""
    Class for storing the 2d evaporation energy distribution. The distribution is of the form
    :math:`E' \, \exp(-E'/theta(E)) / I`.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | U             | Constant that defines the upper energy limit for the              |
    |               | secondary particle so that 0 ≤ E' ≤ E - U.                        |
    +---------------+-------------------------------------------------------------------+
    | thetas        | The :math:`theta(E)` function.                                    |
    +---------------+-------------------------------------------------------------------+

    :param U:               Constant that defines the upper energy limit for the secondary particle so that 0 ≤ E' ≤ E - U.
                            This is not used/needed so should be removed.
    :param thetas:          The :math:`theta(E)` function.
    """

    moniker = 'evaporation'

    def __init__( self, U, thetas ) :

        FunctionalBase.__init__( self, 9, U, thetas )

    @property
    def theta( self ) :
        """Returns a reference to the 1d :math:`\theta(E)` function."""

        return( self.parameter1 )

    def averageEp( self, E ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *E*.

        :param E:       The energy of the projectile energy.

        :return:                Python float.
        """

        if( isinstance( self.parameter1.data, regionsModule.Regions1d ) ) :
            for region in self.parameter1.data :
                if( E <= region[-1][0] ) : break
            theta = region.evaluate( E )
        else :
            theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 180. - a * ( 15. + a ) ) / 270. )
        exp_a = math.exp( -a )
        return( theta * ( 2. - a**2 * exp_a / ( 1. - ( 1. + a ) * exp_a ) ) )

    def evaluate( self, energy ) :
        """
        Returns a 1d function of *self* evaluated at projectile energy *energy*.

        :param energy:      The energy of the outgoing particle.

        :return:            :py:class:`Evaporation1d` instance.
        """

        return( Evaporation1d( energy, self.theta.data.evaluate( energy ), self.U.value ) )

    def energySpectrumAtEnergy( self, energyIn ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn*,

        :param energyIn:            Energy of the projectile.

        :return:                    XYs1d instance for the energy spectrum.
        """

        def A( energyOut, self ) :

            return( self.evaluate( energyOut ) )

        spectrum1 = self.evaluate( energyIn )
        spectrum2 = XYs1dModule.pointwiseXY_C.createFromFunction( [ 0.0, energyIn - self.U.value ], A, spectrum1, 1e-3, 12 )
        return( XYs1d( spectrum2, axes = defaultAxes( energyUnit = self.domainUnit ) ) )

    def sqrtEp_AverageAtE( self, E ) :
        r"""
        This method returns the value of the integral :math:`\int dE' \sqrt{E'} P(E') at projectile energy *E*.

        :param E:   Energy of the projectile.

        :return:    Python float.
        """

        theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( math.sqrt( theta * a ) * ( 252. - a * ( 12. + a ) ) / 315. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        return( math.sqrt( theta ) * ( 1.32934038817913702 * math.erf( sqrt_a ) - sqrt_a * ( 1.5 + a ) * exp_a ) / ( 1. - ( 1. + a ) * exp_a ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        def evaluateAtX( self, x ) :
            """
            A function needed by :py:func:`EnergyFunctionalDataToPointwise` to evaluate an evaporation 
            distribution at the outgoing particle energy *x* for :math:`\theta(x)` = self.p1.

            :param x:       The energy of the outgoing particle.

            :return:        Python float.
            """


            return( x * math.exp( -x / self.p1 ) )

        ef = EnergyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( **kwargs ) )

    @classmethod
    def parseNodeUsingClass(cls, evapElement, xPath, linkData, **kwargs):
        """
        Parse *evapElement* into an instance of *cls*.

        :param cls:         Form class to return.
        :param evapElement: Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *evapElement*.
        """

        xPath.append( evapElement.tag )
        theta_ = Theta.parseNodeUsingClass(evapElement.find(Theta.moniker), xPath, linkData, **kwargs)
        U = physicalQuantityModule.U.parseNodeUsingClass(evapElement.find( 'U' ), xPath, linkData, **kwargs)
        ES = Evaporation( U, theta_ )
        xPath.pop()
        return ES

class Watt( FunctionalBase ) :
    r"""
    Class for storing the 2d Watt energy distribution. The distribution is of the form
    :math:`E' \, \exp(-E'/a(E)) \, \sinh(\sqrt{b(E) \, E'})/ I` where :math:`a(E)` and :math:`b(E)`
    are 1d function of the projectile energy :math:`E`.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | U             | Constant that defines the upper energy limit for the              |
    |               | secondary particle so that 0 ≤ E' ≤ E - U.                        |
    +---------------+-------------------------------------------------------------------+
    | a             | The :math:`a(E)` function.                                        |
    +---------------+-------------------------------------------------------------------+
    | b             | The :math:`b(E)` function.                                        |
    +---------------+-------------------------------------------------------------------+

    :param U:               Constant that defines the upper energy limit for the secondary particle so that 0 ≤ E' ≤ E - U.
                            This is not used/needed so should be removed.
    :param a:               The :math:`a(E)` function.
    :param b:               The :math:`b(E)` function.
    """

    moniker = 'evaporation'

    moniker = 'Watt'

    def __init__( self, U, a, b ) :

        FunctionalBase.__init__( self, 11, U, a, b )

    def averageEp( self, E ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *E*.

        :param E:       The energy of the projectile energy.

        :return:                Python float.
        """

        a, b = self.parameter1.data.evaluate( E ), self.parameter2.data.evaluate( E )
        domainMax_a  = ( E - self.U.value ) / a
        domainMax_b  = math.sqrt( b * ( E - self.U.value ) )
        ab = a * b
        sqrt_ab = math.sqrt( ab )
        I = 0.25 * math.sqrt( math.pi * ab ) * math.exp( 0.25 * ab ) * \
            ( math.erf( math.sqrt( domainMax_a ) - 0.5 * sqrt_ab ) + math.erf( math.sqrt( domainMax_a ) + 0.5 * sqrt_ab ) ) \
            - math.exp( -domainMax_a ) * math.sinh( domainMax_b )
        EI = a * math.sqrt( math.pi * ab ) * ( ab + 6 ) * math.exp( 0.25 * ab ) * \
            ( math.erf( math.sqrt( domainMax_a ) - 0.5 * sqrt_ab ) + math.erf( math.sqrt( domainMax_a ) + 0.5 * sqrt_ab ) ) \
            - 0.25 * a * math.exp( -domainMax_a ) * math.sinh( domainMax_b ) * ( 2. * domainMax_b + ab + 4. + 4. * domainMax_a )
        return( EI / ( 16 * I ) )

    def energySpectrumAtEnergy( self, energyIn ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn*,

        :param energyIn:            Energy of the projectile.

        :return:                    XYs1d instance for the energy spectrum.
        """

        return( self.evaluate( energyIn ) )

    def evaluate(self, energyIn, extrapolation=xDataEnumsModule.Extrapolation.none):
        """
        This method returns a 1d function of *self* evaluated at projectile energy *energyIn*.

        :param energyIn:        The projectile energy to evaluate *self* at.
        :param extrapolation:   This argument is not used.

        :return:                Instance of :py:class:`XYs1d`.
        """

        def A( energyOut, parameters ) :

            a_parameter, b_parameter = parameters
            return( math.exp( -energyOut / a_parameter ) * math.sinh( math.sqrt( b_parameter * energyOut ) ) )

        energyUnit = self.parameter1.data.domainUnit
        energyOut = PQUModule.PQU( 1e-8, 'MeV' ).getValueAs( energyUnit )
        energyOutMax = energyIn - self.U.value
        energyOuts = [ 0.0 ]
        while( 2 * energyOut < energyOutMax ) :
            energyOuts.append( energyOut )
            energyOut *= 10.0
        energyOuts.append( energyOutMax )

        a_parameter = self.parameter1.data.evaluate( energyIn )
        b_parameter = self.parameter2.data.evaluate( energyIn )
        xys1d = XYs1d.createFromFunction( defaultAxes( energyUnit = energyUnit ), energyOuts, A, [ a_parameter, b_parameter ], 1e-3, 10 )

        xys1d.outerDomainValue = energyIn
        xys1d.normalize( insitu = True )

        return( xys1d )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        aMin, aMax = self.parameter1.data.domainMin, self.parameter1.data.domainMax
        if( EMin is None ) : EMin = aMin
        if( EMax is None ) : EMax = aMax
        if( EMin < aMin ) : EMin = aMin
        if( EMax < aMax ) : EMax = aMax
        Es = [ EMin, EMax ]
        for E, a in self.parameter1.data :
            if( E <= EMin ) : continue
            if( E >= EMax ) : continue
            if( E not in Es ) : Es.append( E )
        for E, b in self.parameter2.data :
            if( E <= EMin ) : continue
            if( E >= EMax ) : continue
            if( E not in Es ) : Es.append( E )
        Es.sort( )
        return( Es )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        def evaluateAtX( self, x ) :
            r"""
            A function needed by :py:func:`EnergyFunctionalDataToPointwise` to evaluate a Watt distribution
            distribution at the outgoing particle energy *x* for :math:`\A(E)` = self.p1 and :math:`\B(E)` = self.p2.

            :param x:       The energy of the outgoing particle.

            :return:        Python float.
            """


            return( math.exp( -x / self.p1 ) * math.sinh( math.sqrt( self.p2 * x ) ) )

        ef = EnergyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( **kwargs ) )

    @classmethod
    def parseNodeUsingClass(cls, WattElement, xPath, linkData, **kwargs):
        """
        Parse *WattElement* into an instance of *cls*.

        :param cls:         Form class to return.
        :param WattElement: Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *WattElement*.
        """

        xPath.append( WattElement.tag )
        _a = A.parseNodeUsingClass(WattElement.find(A.moniker), xPath, linkData, **kwargs)
        _b = B.parseNodeUsingClass(WattElement.find(B.moniker), xPath, linkData, **kwargs)
        U = physicalQuantityModule.U.parseNodeUsingClass(WattElement.find( 'U' ), xPath, linkData, **kwargs)
        WS = cls( U, _a, _b )
        xPath.pop()
        return WS

class MadlandNix( FunctionalBase ) :
    """
    Class for storing the 2d Madland/Nix energy distribution.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | EFL           | Average kinetic energy per nucleon of the light fission fragment. |
    +---------------+-------------------------------------------------------------------+
    | EFH           | Average kinetic energy per nucleon of the heavy fission fragment. |
    +---------------+-------------------------------------------------------------------+
    | Ts            | The :math:`T_M(E)` function.                                      |
    +---------------+-------------------------------------------------------------------+

    :param EFL:             Average kinetic energy per nucleon of the light fission fragment.
    :param EFH:             Average kinetic energy per nucleon of the heavy fission fragment.
    :param Ts:              The :math:`T_M(E)` function.
    """

    moniker = 'MadlandNix'

    def __init__( self, EFL, EFH, Ts ) :

        FunctionalBase.__init__( self, 12, None, Ts )
        self.EFL = EFL
        self.EFH = EFH

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        return MadlandNix( self.EFL.copy( ), self.EFH.copy( ), self.parameter1.copy( ) )

    __copy__ = copy

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        if self.EFL.value <= 0 or self.EFH.value <= 0 or self.parameter1.data.rangeMin <= 0:
            warnings.append( warning.MadlandNixBadParameters( self.EFL, self.EFH, self.parameter1.data.rangeMin, self ) )

        return warnings

    def averageEp( self, E ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *E*.

        :param E:       The energy of the projectile energy.

        :return:                Python float.
        """

        unit = self.parameter1.data.axes[-1].unit
        return( 0.5 * ( self.EFL.getValueAs( unit ) + self.EFH.getValueAs( unit ) ) + 4. * self.parameter1.data.evaluate( E ) / 3. )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        FunctionalBase.convertUnits( self, unitMap )
        self.EFL.convertUnits( unitMap )
        self.EFH.convertUnits( unitMap )

    def energySpectrumAtEnergy( self, energyIn ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn*,

        :param energyIn:            Energy of the projectile.

        :return:                    XYs1d instance for the energy spectrum.
        """

        return( self.evaluate( energyIn ) )

    def evaluate(self, energyIn, extrapolation=xDataEnumsModule.Extrapolation.none):
        """
        This method returns a 1d function of *self* evaluated at projectile energy *energyIn*.

        :param energyIn:        The projectile energy to evaluate *self* at.
        :param extrapolation:   This argument is not used.

        :return:                Instance of :py:class:`XYs1d`.
        """

        from numericalFunctions import specialFunctions

        def MadlandNixFunc( Ep, parameters ) :

            def g( Ep, E_F, T_M ) :

                u1 = ( math.sqrt( Ep ) - math.sqrt( E_F ) )
                u1 *= u1 / T_M
                u2 = ( math.sqrt( Ep ) + math.sqrt( E_F ) )
                u2 *= u2 / T_M
                E1 = 0                      # u1^3/2 * E1 is zero for u1 = 0. but E1 is infinity, whence, the next test.
                if( u1 != 0 ) : E1 = specialFunctions.exponentialIntegral( 1, u1 )
                E2 = specialFunctions.exponentialIntegral( 1, u2 )
                complementary = ( u1 > 2. )
                gamma1 = specialFunctions.incompleteGamma( 1.5, u1, complementary )
                gamma2 = specialFunctions.incompleteGamma( 1.5, u2, complementary )
                signG = 1
                if( complementary ) : signG = -1
                return( ( u2 * math.sqrt( u2 ) * E2 - u1 * math.sqrt( u1 ) * E1 + signG * ( gamma2 - gamma1 ) ) / ( 3 * math.sqrt( E_F * T_M ) ) )

            EFL, EFH, T_M = parameters
            return( 0.5 * ( g( Ep, EFL, T_M ) + g( Ep, EFH, T_M ) ) )

        energyUnit = self.parameter1.data.domainUnit

        energyOut = PQUModule.PQU( 1e-8, 'MeV' ).getValueAs( energyUnit )
        energyOutMax = PQUModule.PQU( 40.0, 'MeV' ).getValueAs( energyUnit )
        energyOuts = [ 0.0 ]
        while( 2 * energyOut < energyOutMax ) :
            energyOuts.append( energyOut )
            energyOut *= 10.0
        energyOuts.append( energyOutMax )

        EFL = self.EFL.getValueAs( energyUnit )
        EFH = self.EFH.getValueAs( energyUnit )
        T_M = self.parameter1.data.evaluate( energyIn )
        parameters = [ EFL, EFH, T_M ]
        xys1d = XYs1d.createFromFunction( defaultAxes( energyUnit = energyUnit ), energyOuts, MadlandNixFunc, parameters, 1e-3, 10 )

        xys1d.outerDomainValue = energyIn
        xys1d.normalize( insitu = True )

        return( xys1d )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        return( [ x for x, y in self.parameter1.data ] )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        pwl = XYs2d( axes = defaultAxes( energyUnit = self.parameter1.data.axes[0].unit ) )
        for E, T_M in self.parameter1.data :        # This logic ignores the interpolation of parameter1 as the only two subforms in ENDF/B-VII shows 
            pwl.append( self.evaluate( E ) )        # that linear-linear is better than the 'log-log' given in the ENDF/B-VII/data.

        return( pwl )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* as a :py:class:`XYs2d` representation with the
        P(E') data as :py:class:`Xs_pdf_cdf1d` instances.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
                    
        :return:                :py:class:`XYs2d`
        """

        linear = self.toPointwise_withLinearXYs( )
        return( linear.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

    @classmethod
    def parseNodeUsingClass(cls, MNelement, xPath, linkData, **kwargs):
        """
        Parse *MNelement* into an instance of *cls*.

        :param cls:         Form class to return.
        :param MNelement:   Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *MNelement*.
        """

        xPath.append( MNelement.tag )
        tm = T_M.parseNodeUsingClass(MNelement.find(T_M.moniker), xPath, linkData, **kwargs)
        EFL = physicalQuantityModule.EFL.parseNodeUsingClass(MNelement.find( "EFL" ), xPath, linkData, **kwargs)
        EFH = physicalQuantityModule.EFH.parseNodeUsingClass(MNelement.find( "EFH" ), xPath, linkData, **kwargs)
        MN = cls( EFL, EFH, tm )
        xPath.pop()
        return MN

class NBodyPhaseSpace( Subform ) :
    """
    Class for storing the 2d N-body phase space energy distribution in the center-of-mass frame.

    The following table list the primary members of this class:
    
    +-------------------+-------------------------------------------------------------------+
    | Member            | Description                                                       |
    +===================+===================================================================+
    | numberOfProducts  | Number of outgoing particles.                                     |
    +-------------------+-------------------------------------------------------------------+
    | mass              | Mass of the *numberOfProducts* outgoing particle.                 |
    +-------------------+-------------------------------------------------------------------+

    :param numberOfProducts:    Number of outgoing particles.
    :param mass:                Mass of the *numberOfProducts* outgoing particle.
    """

    moniker = 'NBodyPhaseSpace'

    def __init__( self, numberOfProducts, mass ) :

        Subform.__init__( self )
        self.numberOfProducts = numberOfProducts
        if( not( isinstance( mass, physicalQuantityModule.Mass ) ) ) :
            TypeError( 'mass instance must be a physicalQuantityModule.Mass' )
        self.mass = mass 

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """
        pass

    def copy( self ) :
        """
        Returns a copy of *self*.
        """

        return NBodyPhaseSpace( self.numberOfProducts, self.mass.copy( ) )

    __copy__ = copy

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        #if ( abs( self.mass - info['reactionSuite'].products['n'].getMass('amu') * self.numberOfProducts ) >
        #        self.mass * 0.1 ) :    # return warning?

        return []

    def averageEp( self, E, massUnit, projectileMass, targetMass, productMass, Q ) :
        """
        Calculate the average energy of the product in the center-of-mass frame for projectile energy E. This method only works for
        a one-step reaction.

        :param E:                   The energy of the projectile.
        :param massUnit:            The unit for *projectileMass*, *targetMass* and *productMass*.
        :param projectileMass:      The mass of the projectile in unit of *massUnit*.
        :param targetMass:          The mass of the target in unit of *massUnit*.
        :param productMass:         The mass of the product in unit of *massUnit*.
        :param Q:                   The Q-value for the reaction.

        :return:                Python float.
        """

        M = self.mass.getValueAs( massUnit )
        Ea = targetMass / ( targetMass + projectileMass ) * E + Q
        return( Ea * ( M - productMass ) / ( M * ( self.numberOfProducts - 1 ) ) )

    def evaluate(self, energy_in, **kwargs):
        """
        This method returns the P(E) at projectile energy *energy_in*.
        The accuracy of the returned P(E) can be set by the value of the 'accuracy' key in *kwargs*.

        :param energy_in:       The energy of the projectile.
        :param kwargs:          A dictionary that contains data to control the way this method acts.
        """

        from fudge import product as productModule
        from fudge import reactionSuite as reactionSuiteModule
        from fudge import outputChannel as outputChannelModule
        from fudge.reactions import reaction as reactionModule
        from fudge.core.math import fudgemath

        class Tester:
            """
            This class is used by :py:func:`fudgemath.thickenXYList` to add points until the desired accuracy is reached.
            """

            def __init__(self, relativeTolerance, absoluteTolerance, n, EMax_i):
                """
                :param relativeTolerance:   The relative desired accuracy.
                :param absoluteTolerance:   The absolute desired accuracy.
                :param n:
                :param EMax_i:              The maximum energy of the outgoing product.
                """

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance
                self.n = n
                self.EMax_i = EMax_i

            def evaluateAtX(self, energy_out):

                return math.sqrt(energy_out) * math.pow((self.EMax_i - energy_out), 1.5 * self.n - 4.0)

        accuracy = kwargs.get('accuracy', 1e-3)

        mass, massUnit = self.mass.value, self.mass.unit
        product = self.findClassInAncestry(productModule.Product)
        productMass = product.getMass(massUnit)

        reactionSuite = self.findClassInAncestry(reactionSuiteModule.ReactionSuite)
        projectile = reactionSuite.PoPs.final(reactionSuite.projectile)
        projectileMass = projectile.getMass(massUnit)

        targetID = reactionSuite.target
        target = reactionSuite.PoPs.final(reactionSuite.target)
        targetMass = target.getMass(massUnit)

        outputChannel = self.findClassInAncestry(outputChannelModule.OutputChannel)
        Q = outputChannel.Q.getConstant()

        Ea = targetMass / (targetMass + projectileMass) * energy_in + Q
        EMax_i = Ea * (mass - productMass) / mass
        if EMax_i < 0:
            EMax_i = 1e-5                                    # This is a kludge

        tester = Tester(accuracy, 1e-10, self.numberOfProducts, EMax_i)
        tester.absoluteTolerance = 1e-10 * tester.evaluateAtX(0.5 * EMax_i)
        data = fudgemath.thickenXYList([[0., 0.], [EMax_i, 0.]], tester, biSectionMax=10)
        data = XYs1d(data, outerDomainValue=energy_in)
        data.normalize(insitu=True)

        return data

    def fixDomains(self, labels, energyMin, energyMax):
        """This method does nothing."""

        return 0

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        from fudge import product as productModule
        from fudge import reactionSuite as reactionSuiteModule
        from ...reactions import reaction as reactionModule
        from ... import outputChannel as outputChannelModule
        from fudge.core.math import fudgemath

        class Tester :

            def __init__( self, relativeTolerance, absoluteTolerance, n ) :

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance
                self.n = n
                self.setEMax_i( 1 )

            def evaluateAtX(self, x):

                return math.sqrt(x) * math.pow((self.EMax_i - x), 1.5 * self.n - 4.0)

            def setEMax_i( self, EMax_i ) :

                self.EMax_i = EMax_i

        p = self.findClassInAncestry(productModule.Product)
        mass, massUnit = self.mass.value, self.mass.unit
        productMass = p.getMass(massUnit)

        r = self.findClassInAncestry(reactionModule.Reaction)
        EMin = r.domainMin
        EMax = r.domainMax
        energyUnit = r.domainUnit

        reactionSuite = self.findClassInAncestry(reactionSuiteModule.ReactionSuite)
        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        projectileMass = projectile.getMass(massUnit)
        targetID = reactionSuite.target
        if targetID in reactionSuite.PoPs.aliases: targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]
        targetMass = target.getMass(massUnit)

        c = self.findClassInAncestry(outputChannelModule.OutputChannel)
        Q = c.Q.getConstant( )

        axes = defaultAxes( energyUnit )
        pwl = XYs2d( axes=axes )

        accuracy = kwargs.get( 'accuracy', XYs1dModule.defaultAccuracy )

        t = Tester( accuracy, 1e-10, self.numberOfProducts )
        n = 21
        f = math.pow( EMax / EMin, 1. / n )
        E_ins = [ EMin * f**idx for idx in range( n ) ]
        E_ins[-1] = EMax        # Fix possible round off issue.
        for E_in in E_ins:
            Ea = targetMass / ( targetMass + projectileMass ) * E_in + Q
            EMax_i = Ea * ( mass - productMass ) / mass
            if( EMax_i < 0 ) : EMax_i = 1e-5                # This is a kludge
            t.setEMax_i( EMax_i )
            t.absoluteTolerance = 1e-10 * t.evaluateAtX( 0.5 * EMax_i )
            data = fudgemath.thickenXYList( [ [ 0., 0. ], [ EMax_i, 0. ] ], t, biSectionMax = 10 )
            data = XYs1d( data, outerDomainValue = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return pwl

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent, 
                extraAttributesAsStrings = ' numberOfProducts="%s"' % ( self.numberOfProducts ), emptyTag = False ) ]
        xmlString += self.mass.toXML_strList( indent2, **kwargs )
        xmlString[-1] += "</%s>" % self.moniker
        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        mass = physicalQuantityModule.Mass.parseNodeUsingClass(element.find( 'mass' ), xPath, linkData, **kwargs)
        return cls( int( element.get( "numberOfProducts" ) ), mass )

class Weighted1d :
    """
    This class represents an :py:class:`Weighted` instance evaluated at the specified projectile energy.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | weight        | A Python float the multiplies the functional.                     |
    +---------------+-------------------------------------------------------------------+
    | functional    | A 1d energy distribution of :math:`P(E')`.                        |
    +---------------+-------------------------------------------------------------------+

    :param weight:              The weight to apply to the functional to get the final :math:`f(E')` which my not be nornalized,
    :param functional:          A 1d energy distribution of :math:`P(E')`.
    """

    def __init__( self, weight, functional ) :

        self.weight = weight
        self.functional = functional

    @property
    def domainMin( self ) :
        """
        Returns the minimum output particle energy for *self*.
        """

        return( self.functional.domainMin )

    @property
    def domainMax( self ) :
        """
        Returns the maximum output particle energy for *self*.
        """

        return( self.functional.domainMax )

    def evaluate( self, energy ) :
        """
        Returns the value of *self* evaluated at outgoing particle energy *energy*.

        :param energy:      The energy of the outgoing particle.

        :return:            Python float.
        """

        return( self.weight * self.functional.evaluate( energy ) )

    def integrate( self, energyMin, energyMax ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param energyMin:           The minimum value for the outgoing energy range.
        :param energyMax:           The maximum value for the outgoing energy range.

        :return:                    A float representing the value of the integration.
        """

        return( self.weight * self.functional.integrate( energyMin, energyMax ) )

class WeightedFunctionals1d :
    """
    This class represents an :py:class:`WeightedFunctionals` instance evaluated at the specified projectile energy.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | weighted1d    | The list of :py:class:`Weighted1d` instances that make the        |
    |               | complete, normalized :math:`P(E')`.                               |
    +---------------+-------------------------------------------------------------------+

    :param weighted1d:      The list of :py:class:`Weighted1d` instances that make the complete, normalized :math:`P(E')`.
    """

    def __init__( self ) :

        self.weighted1d = []

    @property
    def domainMin( self ) :
        """
        Returns the minimum outgoing particle energy for *self*.
        """

        return( min( weighted1d1.domainMin for weighted1d1 in self.weighted1d ) )

    @property
    def domainMax( self ) :
        """
        Returns the maximum outgoing particle energy for *self*.
        """

        return( max( weighted1d1.domainMax for weighted1d1 in self.weighted1d ) )

    def append( self, _weighted1d ) :
        """
        Appends a :py:class:`Weighted1d` instance to *self*.

        :param _weighted1d: A :py:class:`Weighted1d` instance.
        """

        self.weighted1d.append( _weighted1d )

    def evaluate( self, energy ) :
        """
        Returns the value of *self* evaluated at outgoing particle energy *energy*.

        :param energy:      The energy of the outgoing particle.

        :return:            Python float.
        """

        value = 0.0
        for weighted1d in self.weighted1d : value += weighted1d.evaluate( energy )
        return( value )

    def integrate( self, energyMin, energyMax ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param energyMin:           The minimum value for the outgoing energy range.
        :param energyMax:           The maximum value for the outgoing energy range.

        :return:                    A float representing the value of the integration.
        """

        value = 0.0
        for weighted1d in self.weighted1d : value += weighted1d.integrate( energyMin, energyMax )
        return( value )

class Weighted(ancestryModule.AncestryIO):
    r"""
    This class represents 2d energy distribution math:`P(E'|E)` with a projectile energy dependent
    weight function :math:`w(E)`. That is, the 2d energy distribution is :math:`w(E) \, P(E'|E)`.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | weight        | The projectile energy dependent function :math:`w(E)`.            |
    +---------------+-------------------------------------------------------------------+
    | functional    | The 2d energy distribution :math:`P(E'|E)`.                       |
    +---------------+-------------------------------------------------------------------+

    :param weight:              The weight function.
    :param functional:          The 2d energy distribution :math:`P(E')`.
    """

    moniker = 'weighted'
    ancestryMembers = ( 'weight', 'functional' )

    def __init__( self, weight, functional ) :

        ancestryModule.AncestryIO.__init__( self )

        self.weight = weight
        self.weight.setAncestor( self )
# BRB6 weight should be its own class so that label need not be set here.
        self.weight.label = 'weight'

        self.functional = functional
        self.functional.setAncestor( self )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.weight.convertUnits( unitMap )
        self.functional.convertUnits( unitMap )

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        return Weighted( self.weight.copy( ), self.functional.copy( ) )

    __copy__ = copy

    def checkProductFrame( self ) :
        "Calls checkProductFrame on self's functional."

        self.functional.checkProductFrame( )

    def averageEp( self, E ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *E*.

        :param E:           The energy of the projectile energy.

        :return:            Python float.
        """

        return( self.weight.evaluate( E ) * self.functional.averageEp( E ) )

    def evaluate( self, energy ) :
        """
        Returns the value of *self* evaluated at projectile energy *energy*.

        :param energy:      The energy of the projectile.

        :return:            A :py:class:`Weighted1d` instance..
        """

        return( Weighted1d( self.weight.evaluate( energy ), self.functional.evaluate( energy ) ) )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        Calls the **fixDomains** for the **weight** and **functional** members.
        """

        numberOfFixes  = self.weight.fixDomains(energyMin, energyMax, domainToFix)
        numberOfFixes += self.functional.fixDomains(energyMin, energyMax, domainToFix)

        return  numberOfFixes

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        energyArray = self.functional.getEnergyArray( EMin, EMax )
        for x, y in self.weight :
            if( x not in energyArray ) : energyArray.append( x )
        energyArray.sort( )
        return( energyArray )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s>' %(  indent, self.moniker ) ]
        xmlString += self.weight.toXML_strList( indent2, **kwargs )
        xmlString += self.functional.toXML_strList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        weights, functional = node
        _weight = XYs1dModule.XYs1d.parseNodeUsingClass(weights, xPath, linkData, **kwargs)
        subformClass = {
                XYs2d.moniker                   : XYs2d,
                GeneralEvaporation.moniker      : GeneralEvaporation,
                Watt.moniker                    : Watt,
                MadlandNix.moniker              : MadlandNix,
                SimpleMaxwellianFission.moniker : SimpleMaxwellianFission,
                Evaporation.moniker             : Evaporation,
            }.get( functional.tag )
        functional = subformClass.parseNodeUsingClass(functional, xPath, linkData, **kwargs)
        instance = Weighted(_weight, functional)

        xPath.pop( )

        return instance

class WeightedFunctionals( Subform ) :
    """
    This class represents the :math:`P(E'|E)` distribution as a list of :py:class:`Weighted` instances.

    The following table list the primary members of this class:
    
    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | weights       | The list of :py:class:`Weighted` instances that make the          |
    |               | complete, normalized distribution :math:`P(E'|E)`.                |
    +---------------+-------------------------------------------------------------------+
    """

    moniker = 'weightedFunctionals'
    ancestryMembers = ( 'weights', )

    def __init__( self ) :

        Subform.__init__( self )
        self.weights = []

    def __len__( self ) :
        """Returns the number of :py:class:`Weighted` instances in *self*."""

        return( len( self.weights ) )

    def __getitem__( self, i ) :
        """
        Returns the (i-1)^th :py:class:`Weighted` instance in *self*.

        :param i:       The index of the :py:class:`Weighted` instance to return.

        :return:        A :py:class:`Weighted` instance.
        """

        return( self.weights[i] )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        for weight in self.weights : weight.convertUnits( unitMap )

    def copy( self ) :
        """
        Returns a copy of *self*.
        """

        newWF = WeightedFunctionals( )
        for weight in self : newWF.addWeight( weight.copy( ) )
        return( newWF )

    __copy__ = copy

    def addWeight( self, weight ) :
        """
        Appends a :py:class:`Weighted` instance to *self*.

        :param weight:  A :py:class:`Weighted` instance.
        """

# BRB. either we remove this class (WeightedFunctionals) or self.weights should be a suite like object and it should be
# the ancestor of weight (i.e., self.weights[-1].setAncestor( self.weights ) and not weight.setAncestor( self ).
        if( not( isinstance( weight, Weighted ) ) ) : raise Exception( 'Invalid weight of type %s' % brb.getType( weight ) )
        self.weights.append( weight )
        weight.setAncestor( self )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        totalWeight = sum( [w.weight for w in self.weights] )
        if (totalWeight.rangeMin != 1.0) and (totalWeight.rangeMax != 1.0):
            warnings.append( warning.WeightsDontSumToOne( obj=self ) )

        for weight in self:
            warnings += weight.functional.check( info )

        return warnings

    def averageEp( self, E ) :
        """
        Calculated the average energy of the outgoing photon at projectile energy *E*.

        :param E:           The energy of the projectile energy.

        :return:            Python float.
        """

        Ep = 0
        for weight in self : Ep += weight.averageEp( E )
        return( Ep )

    def evaluate( self, energy ) :
        """
        This method returns a 1d functional representation of *self* evaluated at projectile energy *energy*.

        :param energy:      The projectile energy to evaluate *self* at.

        :return:            A :py:class:`WeightedFunctionals1d` instance.
        """

        weightedFunctionals1d1 = WeightedFunctionals1d( )
        for weight in self : weightedFunctionals1d1.append( weight.evaluate( energy ) )
        return( weightedFunctionals1d1 )

    def energySpectrumAtEnergy( self, energyIn ) :
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn*,

        :param energyIn:            Energy of the projectile.

        :return:                    XYs1d instance for the energy spectrum.
        """

        xys1d = XYs1d( axes = defaultAxes( energyUnit = self.weights[0].functional.domainUnit ) )
        for weight in self :
            spectrum1 = weight.weight.evaluate( energyIn ) * weight.functional.energySpectrumAtEnergy( energyIn )
            xys1d += spectrum1
        return( xys1d )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        Calls the **fixDomains** for each **weight** of *self*.
        """

        numberOfFixes = 0
        for weight in self: numberOfFixes += weight.fixDomains(energyMin, energyMax, domainToFix)

        return  numberOfFixes

    def sqrtEp_AverageAtE( self, E ) :
        """
        This method has not been implemented. It returns None so the method uncorrelated.calculateAverageProductData will still work
        when calculating the momentum deposition for isotropic scattering in the lab frame, but will execute a raise otherwise.

        :param E:       The energy of the projectile.

        :return:        Python None
        """

        return( None )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        energyArray = []
        for weight in self :
            for energy in weight.getEnergyArray( EMin, EMax ) :
                if( energy not in energyArray ) : energyArray.append( energy )
        energyArray.sort( )
        return( energyArray )

    def to_xs_pdf_cdf1d(self, style, tempInfo, indent):

        count = sum([1 if isinstance(weighted.functional, XYs2d) else 0 for weighted in self.weights])
        if count == 0:
            return None

        if count != len(self.weights):
            raise Exception('This needs work!.')

        weightedFunctionals = WeightedFunctionals()
        for weighted in self.weights:
            functional = weighted.functional.to_xs_pdf_cdf1d(style, tempInfo, indent)
            weightedFunctionals.addWeight(Weighted(weighted.weight.copy( ), functional))

        return weightedFunctionals

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        if( len( self ) > 2 ) : raise Exception( 'more than two weights currently not supported' )
        E_ins, data = [], []
        for weighted_ in self :
            w = weighted_.weight.toPointwise_withLinearXYs( **kwargs )
            e = weighted_.functional.toPointwise_withLinearXYs( **kwargs )
            data.append( [ w, e ] )
            for x, y in w :
                if( x not in E_ins ) : E_ins.append( x )
            for x in e :
                if( x.outerDomainValue not in E_ins ) : E_ins.append( x.outerDomainValue )
        E_ins.sort( )
        pwl = XYs2d( axes = defaultAxes( self[0].weight.domainUnit ) )
        for E_in in E_ins :
            wv1, ev1 = data[0]
            wv2, ev2 = data[1]
            w1, w2 = wv1.evaluate( E_in ), wv2.evaluate( E_in )
            if( w1 == 1 ) :
                e = ev1.getSpectrumAtEnergy( E_in )
            elif( w2 == 1 ) :
                e = ev2.getSpectrumAtEnergy( E_in )
            else :
                raise Exception( 'One weight must be zero: %s, %s' % ( w1, w2 ) )
            e.outerDomainValue = E_in
            e.normalize( insitu = True )
            pwl.append( e )
        return( pwl )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        for _subform in self.weights : xmlString += _subform.toXML_strList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        WF = cls()

        for child in node: WF.addWeight(Weighted.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop( )

        return WF

class EnergyFunctionalDataToPointwise :
    """
    For internal use. This class is used by the methods toPointwise_withLinearXYs of the 2d functional energy distribution
    to convert the functional into an :py:class:`XYs2d` instance.

    :param data:            An instance of the 2d functional energy distribution.
    :param evaluateAtX:     A function that evaluates the 2d functional energy distribution at a specific projectile energy.
    """

    def __init__( self, data, evaluateAtX ) :

        self.data = data
        self.evaluateAtX = evaluateAtX

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        from fudge.core.math import fudgemath

        class Tester :

            def __init__( self, relativeTolerance, absoluteTolerance, evaluateAtX ) :

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance
                self.evaluateAtX_ = evaluateAtX
                self.setParameter1( 1. )
                self.setParameter2( 1. )

            def evaluateAtX( self, x ) :

                return( self.evaluateAtX_( self, x ) )

            def setParameter1( self, p1 ) :

                self.p1 = p1

            def setParameter2( self, p2 ) :

                self.p2 = p2

        accuracy = xDataBaseModule.getArguments( kwargs, { 'accuracy' : XYs1dModule.defaultAccuracy } )['accuracy']

        parameter1 = self.data.parameter1.data.toPointwise_withLinearXYs( **kwargs )
        axes = defaultAxes( self.data.domainUnit )
        pwl = XYs2d( axes = axes )
# BRB6 hardwired
        one_eV = PQUModule.PQU( '1 eV' ).getValueAs( axes[-1].unit )

        t = Tester( accuracy, 1e-10, self.evaluateAtX )
        parameter2 = self.data.parameter2
        if( parameter2 is not None ) : parameter2 = parameter2.data.toPointwise_withLinearXYs( **kwargs )
        for index, E_in_p1 in enumerate( parameter1 ) :
            E_in, p1 = E_in_p1
            EpMax = E_in - self.data.U.value
            EpMin = 0.              # Only used with debugging.
            if( EpMax == 0. ) :
                if( index != 0 ) : raise Exception( "i = %d, E_in = %s, U = %s" % ( index, E_in, self.data.U.value ) )
                EpMax = E_in * 1e-6
                if( EpMax == 0. ) :
                    EpMax = parameter1[1][1] * 1e-6         # This and the next line are arbitary.
                if( EpMax > 1e-3 ) : EpMax = 2e-5
                data = [ [ EpMin, 1. ], [ EpMax, 1e-10 ] ]
            else :
                t.setParameter1( p1 )
                if( parameter2 is not None ) : t.setParameter2( parameter2.evaluate( E_in ) )
                level, data = 0, [ [ EpMin, t.evaluateAtX( EpMin ) ], [ EpMax, t.evaluateAtX( EpMax ) ] ]
                while( data[1][1] < 1e-10 ) :               # Fill in some tail points while they are < 1e-10
                    level, EpMax = level + 1, 0.5 * EpMax
                    data.insert( 1, [ EpMax, t.evaluateAtX( EpMax ) ] )
                    if( level > 10 ) : break
                if( data[0][0] < one_eV ) :
                    if( data[1][0] > 1e6 * one_eV ) : data.insert( 1,  [ 1e6 * one_eV, t.evaluateAtX( 1e6 * one_eV ) ] )
                    if( data[1][0] > 1e3 * one_eV ) : data.insert( 1,  [ 1e3 * one_eV, t.evaluateAtX( 1e3 * one_eV ) ] )
                    if( data[1][0] > one_eV ) : data.insert( 1,  [ one_eV, t.evaluateAtX( one_eV ) ] )
                data = fudgemath.thickenXYList( data, t, biSectionMax = 10 )
            data = XYs1d( data = data, outerDomainValue = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

class Form( baseModule.Form ) :
    """
    This class stores the enegy distribution part (i.e., :math:`P(E'|E)`) for an uncorrelated distribution.

    The following table list the primary members of this class:
    
    +-------------------+-----------------------------------------------------------+
    | Member            | Description                                               |
    +===================+===========================================================+
    | energySubform     | 2d energy distribution.                                   |
    +-------------------+-----------------------------------------------------------+
    
    :param label:               The label for this form.
    :param productFrame:        The frame the product data are specified in.
    :param energySubform:       An :py:class:`Subform` representing the energy data.
    """

    moniker = 'energy'
    subformAttributes = ( 'energySubform', )

    def __init__( self, label, productFrame, energySubform ) :

        if( not( isinstance( energySubform, Subform ) ) ) : raise TypeError( 'instance is not an energy subform' )
        baseModule.Form.__init__( self, label, productFrame, ( energySubform, ) )

    @classmethod
    def parseNodeUsingClass(cls, energyElement, xPath, linkData, **kwargs):
        """
        Parse *energyElement* into an instance of *cls*.

        :param cls:             Form class to return.
        :param energyElement:   Node to parse.
        :param xPath:           List containing xPath to current node, useful mostly for debugging.
        :param linkData:        dict that collects unresolved links.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *energyElement*.
        """

        subformClass = {
                DiscreteGamma.moniker :             DiscreteGamma,
                PrimaryGamma.moniker :              PrimaryGamma,
                XYs2d.moniker :                     XYs2d,
                Regions2d.moniker :                 Regions2d,
                GeneralEvaporation.moniker :        GeneralEvaporation,
                Watt.moniker :                      Watt,
                MadlandNix.moniker :                MadlandNix,
                SimpleMaxwellianFission.moniker :   SimpleMaxwellianFission,
                Evaporation.moniker :               Evaporation,
                WeightedFunctionals.moniker :       WeightedFunctionals,
                NBodyPhaseSpace.moniker :           NBodyPhaseSpace,
            }.get( energyElement.tag )
        if( subformClass is None ) : raise Exception( "encountered unknown energy subform: %s" % energyElement.tag )

        energyForm = subformClass.parseNodeUsingClass(energyElement, xPath, linkData, **kwargs)

        return( energyForm )
