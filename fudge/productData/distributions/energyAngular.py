# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Energy/angular double differential distribution classes.

This module contains the following classes:
         
    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+
    | XYs1d         | Class used to store f(mu) at a specified E'.                          |
    +---------------+-----------------------------------------------------------------------+
    | Legendre      | Class used to store f(mu) at a specified E' as a Legendre series.     |
    +---------------+-----------------------------------------------------------------------+
    | XYs2d         | Class used to store P(mu,E') at a specified E.                        |
    +---------------+-----------------------------------------------------------------------+
    | Subform       | Base class for :py:class:`XYs3d` and :py:class:`Regions3d`.           |
    +---------------+-----------------------------------------------------------------------+
    | XYs3d         | Class used to store P(E',mu|E).                                       |
    +---------------+-----------------------------------------------------------------------+
    | Regions3d     | Class used to store P(E',mu|E) with multiple regions.                 |
    +---------------+-----------------------------------------------------------------------+
    | Form          | Class representing GNDS eneary/angular data (i.e., P(E',mu|E)).       |
    +---------------+-----------------------------------------------------------------------+
"""

import math

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import series1d as series1dModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule

from PoPs import IDs as PoPsIDsModule

from . import base as baseModule
from . import energy as energyModule
from . import miscellaneous as miscellaneousModule


def defaultAxes( energyUnit ) :
    """
    Returns an :py:class:`axesModule.Axes` instance for energy, angular data (i.e., P(E',mu|E)).

    :param energyUnit:  Unit for the energy data.

    :return:            An :py:class:`Axes` instance.
    """

    axes = axesModule.Axes(4)
    axes[0] = axesModule.Axis( 'P(energy_out,mu|energy_in)', 0, '1/' + energyUnit )
    axes[1] = axesModule.Axis( 'mu', 1, '' )
    axes[2] = axesModule.Axis( 'energy_out', 2, energyUnit )
    axes[3] = axesModule.Axis( 'energy_in', 3, energyUnit )
    return( axes )

class XYs1d( XYs1dModule.XYs1d ) :    # FIXME: BRB, Should this class inherit from angular.XYs1d?
    """
    Class used to store a tabulated f(mu) at a specified E'.
    """

    def averageMu( self ) :
        """
        This method returns the average value of mu for *self*.
        """

        allowedInterpolations = [xDataEnumsModule.Interpolation.linlin, xDataEnumsModule.Interpolation.flat]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations, XYs1dModule.defaultAccuracy )
        return( xys.integrateWithWeight_x( ) )

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.

        :return:            Reference to the :py:class`XYs1d`.
        """

        return( XYs1d )

class Legendre( series1dModule.LegendreSeries ) :
    """
    Class used to store a f(mu) at a specified E' as a Legendre series.
    """

    def averageMu( self ) :
        """
        This function returns the average value of mu for *self*.
        """

        return( self.getCoefficientSafely( 1 ) )

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.

        :return:            Reference to the :py:class`XYs1d`.
        """

        return( XYs1d )

class XYs2d( multiD_XYsModule.XYs2d ) :
    """
    Class used to store a P(E',mu) at a specified E.
    """

    def averageEnergy( self ) :
        """
        The method calculates the average outgoing energy of *self*.
        """

        integralOfMu = [ [ pdfOfMuAtEp.outerDomainValue, pdfOfMuAtEp.integrate() ] for pdfOfMuAtEp in self ]

        return XYs1dModule.XYs1d(integralOfMu, interpolation = self.interpolation, axes=self.axes).integrateWithWeight_x()

    def averageForwardMomentum( self, sqrt_2_ProductMass ) :
        """
        The method calculates the average outgoing momentem of *self* in the direction of the projectile.

        :param sqrt_2_ProductMass:      Twice the mass of the produce in units of energy per c**2.

        :return:                        A float.
        """

        averageMu = [ [ pdfOfMuAtEp.outerDomainValue, pdfOfMuAtEp.averageMu( ) ] for pdfOfMuAtEp in self ]

        return sqrt_2_ProductMass * XYs1dModule.XYs1d(averageMu, interpolation = self.interpolation, axes=self.axes).integrateWithWeight_sqrt_x()

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 2-d data of *self*.

        :return:            Reference to the :py:class`XYs2d`.
        """

        return( XYs2d )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs2d` instance.
        """

        return( ( XYs1d, Legendre ) )

class Subform( baseModule.Subform ) :
    """Abstract base class for :py:class`Form` subforms."""

    pass

class XYs3d( Subform, multiD_XYsModule.XYs3d ) :
    """
    Class used to store a P(E',mu|E).
    """

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )
        Subform.__init__( self )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        verbosity = kwargs.get( 'verbosity', 0 )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        multiplicity = kwargs['multiplicity']
        projectileMass = kwargs['projectileMass']
        targetMass = kwargs['targetMass']
        productMass = kwargs['productMass']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        productFrame = kwargs['productFrame']

        sqrt_2_ProductMass = math.sqrt( 2 * productMass )
#
# FIXME, still need to fill in between energy points.
#
        aveEnergy = []
        aveMomentum = []
        if productFrame == xDataEnumsModule.Frame.lab:
            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * pdfOfEpMuAtE.averageEnergy( ) ] )

            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                momemtum = multiplicity.evaluate( energy ) * pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
                aveMomentum.append( [ energy, momemtum ] )

        else :                              # Center-of-mass.
            const1 = math.sqrt( 2 * projectileMass ) / ( projectileMass + targetMass )
            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                vCOM = const1 * math.sqrt( energy )
                EpCOM = pdfOfEpMuAtE.averageEnergy( )
                ppCOM = pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                productLabEnergy = 0.5 * productMass * vCOM * vCOM + EpCOM + vCOM * ppCOM
                aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * productLabEnergy ] )

            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                vCOM = const1 * math.sqrt( energy )
                productLabForwardMomentum = productMass * vCOM + pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                aveMomentum.append( [ energy, multiplicity.evaluate( energy ) * productLabForwardMomentum ] )

        return( [ aveEnergy ], [ aveMomentum ] )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning
        from pqu import PQU
        warnings = []
        axes = axesModule.Axes(2)
        axes[0] = self.axes[-2]
        axes[1] = self.axes[-1]

        for index, energy_in in enumerate(self):
            integral = energy_in.integrate()
            if abs(integral-1.0) > info['normTolerance']:
                warnings.append(warning.UnnormalizedDistribution(PQU.PQU(energy_in.outerDomainValue, axes[0].unit),
                                                                 index, integral, self.toXLink()))
            energy_outs = {}
            for xy in energy_in:
                if isinstance(xy, Legendre):
                    if xy[0] == 0 and any([coefficient != 0 for coefficient in xy[1:]]):
                        incident = PQU.PQU(energy_in.outerDomainValue, axes[0].unit)
                        outgoing = PQU.PQU(xy.outerDomainValue, axes[1].unit)

                        warnings.append(warning.InvalidLegendreSeries(incident, outgoing, xy))

                energy_outs[xy.toPointwise_withLinearXYs(accuracy=XYs1dModule.defaultAccuracy, upperEps=1e-8).rangeMin] = xy.outerDomainValue
            minProb = min(energy_outs.keys())
            if minProb < 0:
                warnings.append(
                    warning.NegativeProbability(minProb, PQU.PQU(energy_in.outerDomainValue, axes[0].unit),
                                                energy_out=PQU.PQU(energy_outs[minProb], axes[1].unit), obj=energy_in))
            if self.interpolationQualifier is xDataEnumsModule.InterpolationQualifier.unitBase:
                # check for more than one outgoing distribution integrating to 0 at high end of incident energies
                integrals = [eout.integrate() for eout in energy_in]
                for idx, integral in enumerate(integrals[::-1]):
                    if integral != 0: break
                if idx > 1:
                    warnings.append(warning.ExtraOutgoingEnergy(PQU.PQU(energy_in.outerDomainValue, axes[-1].unit),
                                                                obj=energy_in))
        return warnings

    def energySpectrumAtEnergy(self, energyIn, **kwargs):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn* for frame *frame*,

        :param energyIn:            Energy of the projectile.
        :param kwargs:              A dictionary that contains data to control the way this method acts.
        
        :return:                    XYs1d instance for the energy spectrum.
        """

        muMin = kwargs.get('muMin', -1.0)
        muMax = kwargs.get('muMax',  1.0)

        energyUnit = self.domainUnit
        flag, functional2d1, functional2d2, frac, interpolation, interpolationQualifier = self.getBoundingSubFunctions(energyIn)
        if flag in ['<', '>', None]:
            return energyModule.XYs1d(axes=energyModule.defaultAxes(energyUnit))
        data = [[functional1d.outerDomainValue, functional1d.integrate(domainMin=muMin, domainMax=muMax)] for functional1d in functional2d1]
        functional1d = energyModule.XYs1d(data, interpolation = functional2d1.interpolation, axes = energyModule.defaultAxes(energyUnit))
        if flag != '=':
            data = [[functional1d.outerDomainValue, functional1d.integrate(domainMin=muMin, domainMax=muMax)] for functional1d in functional2d2]
            functional1d2 = energyModule.XYs1d(data, interpolation = functional2d2.interpolation, axes = energyModule.defaultAxes(energyUnit))
            frac = (energyIn - functional2d1.outerDomainValue) / (functional2d2.outerDomainValue - functional2d1.outerDomainValue)
            if False:
                functional1d = (1.0 - frac) * functional1d + frac * functional1d2
            else:
                functional1d = XYs1dModule.pointwiseXY_C.unitbaseInterpolate(energyIn, functional2d1.outerDomainValue, functional1d.nf_pointwiseXY,
                        functional2d2.outerDomainValue, functional1d2.nf_pointwiseXY, 1)
                functional1d = energyModule.XYs1d(functional1d, interpolation = functional1d.getInterpolation(), axes = energyModule.defaultAxes(energyUnit))

        return functional1d

    def normalize( self, insitu = True ) :
        """
        This methods returns a renormalized copy (or *self) of *self*.

        :param insitu:      If True, normalizes *self*, else copies *self* and return the normalized copy.

        :return:            Returns the normalzied function which can be *self* if *insitu* is True.
        """

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for E_MuEpPs in n :
            sum = E_MuEpPs.integrate()
            for muEpPs in E_MuEpPs : muEpPs.setData( muEpPs / sum )
        return( n )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        productFrame = tempInfo['productFrame']

        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))
        TM_1, TM_E = transferMatricesModule.EEpMuP_TransferMatrix( style, tempInfo, productFrame, tempInfo['crossSection'], self,
            tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs3d` instance.
        """

        return( multiD_XYsModule.XYs3d.toPointwise_withLinearXYs( self, cls = XYs3d, **kwargs ) )

    def to_xs_pdf_cdf1d(self, style, tempInfo, indent):
        """
        This method returns a copy of *self* as an :py:class:`energyAngularMCModule.Energy` representation
        and an :py:class:`energyAngularMCModule.EnergyAngular` representation. This is, P(E',mu|E) becomes
        P(E'|E) * P(mu|E,E') with the E' data in P(E'|E) and the mu data in P(mu|E,E') being represented with
        :py:class:`energyAngularMCModule.Energy` and :py:class:`energyAngularMCModule.EnergyAngular` instances, respectively.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                :py:class:`energyAngularMCModule.Energy` and :py:class:`energyAngularMCModule.EnergyAngular` instances.
        """

        from . import energyAngularMC as energyAngularMCModule

        energy = energyModule.XYs2d(axes=self.axes, interpolation=self.interpolation, interpolationQualifier=self.interpolationQualifier)
        for PofMuGivenEp in self:
            data = energyModule.XYs1d([[PofMu.outerDomainValue, PofMu.integrate()] for PofMu in PofMuGivenEp], interpolation=PofMuGivenEp.interpolation)
            energy.append(energyModule.Xs_pdf_cdf1d.fromXYs(data, outerDomainValue=PofMuGivenEp.outerDomainValue, thinEpsilon=1e-14))

        xys3d = energyAngularMCModule.XYs3d(axes=self.axes)
        for index, PofMuGivenEp in enumerate(self):
            xys2d = energyAngularMCModule.XYs2d(outerDomainValue=PofMuGivenEp.outerDomainValue)
            EPrimes = energy[index].xs.values.values        # The outgoing energies have been thinned to thinEpsilon in energyModule.Xs_pdf_cdf1d.fromXYs above.
            for PofMu in PofMuGivenEp:
                if PofMu.outerDomainValue not in EPrimes:
                    continue
                _PofMu = PofMu
                if isinstance(PofMu, Legendre):
                    if PofMu[0] == 0:
                        _PofMu = XYs1d([[-1, 0.5], [1, 0.5]])
                    _PofMu = _PofMu.toPointwise_withLinearXYs(accuracy=XYs1dModule.defaultAccuracy, upperEps=1e-8)
                xys2d.append(energyAngularMCModule.Xs_pdf_cdf1d.fromXYs(_PofMu, PofMu.outerDomainValue, thinEpsilon=1e-14))
            xys3d.append(xys2d)

        return energyAngularMCModule.Energy(energy), energyAngularMCModule.EnergyAngular(xys3d)

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs3d` instance.
        """

        return( ( XYs2d, ) )

class Regions3d( Subform, regionsModule.Regions3d ) :
    """
    Class used to store a P(E',mu|E).
    """

    def __init__( self, **kwargs ) :

        regionsModule.Regions3d.__init__( self, **kwargs )
        Subform.__init__( self )

    def calculateAverageProductData(self, style, indent='', **kwargs):
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        aveEnergies = []
        aveMomenta = []
        for region in self:
            aveEnergy, aveMomentum = region.calculateAverageProductData(style, indent=indent, **kwargs)
            aveEnergies.append(aveEnergy)
            aveMomenta.append(aveMomentum)

        print(aveEnergies)
        return aveEnergies, aveMomenta

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        raise Exception( 'Not implemented' )

    def normalize( self, insitu = True ) :
        """
        This methods returns a renormalized copy (or *self) of *self*.

        :param insitu:      If True, normalizes *self*, else copies *self* and return the normalized copy.

        :return:            Returns the normalzied function which can be *self* if *insitu* is True.
        """

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for region in n : region.normalize( insitu = True )
        return( n )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Currently not implemented so executes a raise.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`energyAngularModule.XYs3d` instance.
        """

        raise Exception( 'Not implemented' )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`Regions3d` instance.
        """

        return( ( XYs3d, ) )

class Form( baseModule.Form ) :
    """
    This class stores a P(E',mu|E) distribution.
    The following table list the primary members of this class:

    +-----------------------+-----------------------------------------------------------+
    | Member                | Description                                               |
    +=======================+===========================================================+
    | energyAngularSubform  | This member stores the actual data as a 3d function.      |
    +-----------------------+-----------------------------------------------------------+
        
    :param label:                   The label for this form.
    :param productFrame:            The frame the product data are specified in.
    :param energyAngularSubform:    A :py:class:`Subform` instance.
    """

    moniker = 'energyAngular'
    subformAttributes = ( 'energyAngularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, energyAngularSubform ) :

        if( not( isinstance( energyAngularSubform, Subform ) ) ) : raise TypeError( 'instance is not an energyAngular subform' )
        baseModule.Form.__init__( self, label, productFrame, ( energyAngularSubform, ) )

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy given in the P(E',mu|E) data.
        """

        return( self.energyAngularSubform.domainMin )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy given in the P(E',mu|E) data.
        """

        return( self.energyAngularSubform.domainMax )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile for the data.
        """

        return( self.energyAngularSubform.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        kwargs['productFrame'] = self.productFrame
        return( self.energyAngularSubform.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn* for frame *frame*,

        :param energy_in:           Energy of the projectile.
        :param frame:               The frame to calculate the energy spectrum in. 
        :param kwargs:              A dictionary that contains data to control the way this method acts.
        
        :return:                    XYs1d instance for the energy spectrum.
        """

        if self.product.pid == PoPsIDsModule.photon:
            frame = self.productFrame       # Kludge for now.
        if self.productFrame == frame:
            return self.energyAngularSubform.energySpectrumAtEnergy(energyIn, **kwargs)
        else :
            if self.productFrame == xDataEnumsModule.Frame.lab:
                raise TypeError( 'Lab to center-of-mass translation not supported.' )

            muMin = kwargs.get('muMin', -1.0)
            muMax = kwargs.get('muMax',  1.0)
            xys2d = self.spectrumAtEnergy(energyIn, xDataEnumsModule.Frame.lab)
            data = [[xys1d.outerDomainValue, xys1d.integrate(domainMin=muMin, domainMax=muMax)] for xys1d in xys2d]
            return energyModule.XYs1d(data, axes=energyModule.defaultAxes(self.domainUnit), interpolation=xys2d.interpolation)

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        This method call *fixDomains* on the *energyAngularSubform* member.
        """

        return self.energyAngularSubform.fixDomains(energyMin, energyMax, domainToFix, tweakLower=True, epsilon=2e-2)
        
    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param reaction_suite:      The :py:class:`ReactionSuite` instance for this distribution.
        :param energyIn:            The energy of the projectile.
        :param energyOut:           The outgoing energy range to integrate over.
        :param muOut:               The outgoing mu range to integrate over.
        :param phiOut:              The outgoing phi range to integrate over.
        :param frame:               The frame the outgoing energy, mu and phi range specify.
        :param LegendreOrder:       The parameter is not used.

        :return:                    A float representing the value of the integration.
        """

# Need to check frame.

        data = self.energyAngularSubform
        if( energyIn < data.domainMin ) : return( 0.0 )
        if( energyIn > data.domainMax ) : energyIn = self.domainMax

        XYs2d = data.evaluate( energyIn )
        points = []
        for function1d in XYs2d :
            muOutMin, muOutMax = miscellaneousModule.domainLimits( muOut, function1d.domainMin, function1d.domainMax )
            if( muOutMax is None ) :
                muPartialIntegral = function1d.evaluate( muOutMin )
            else :
                muPartialIntegral = function1d.integrate(muOutMin, muOutMax)
            points.append( [ function1d.outerDomainValue, muPartialIntegral ] )
        xys = XYs1d( data = points )

        energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, xys.domainMin, xys.domainMax )
        if( energyOutMax is None ) :
            energyOutMuEvaluate = xys.evaluate( energyOutMin )
        else :
            energyOutMuEvaluate = xys.integrate( energyOutMin, energyOutMax )

        phiEvaluate = miscellaneousModule.muPhiEvaluate( None, phiOut )

        return( phiEvaluate * energyOutMuEvaluate )

    def processMC_cdf(self, style, tempInfo, indent):
        """
        This methods returns an :py:class:`energyAngularMCModule.Form` instance representing *self*.

    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                An instance of self.
        """     

        from . import energyAngularMC as energyAngularMCModule

        energy, energyAngular = self.energyAngularSubform.to_xs_pdf_cdf1d(style, tempInfo, indent)
        return( energyAngularMCModule.Form( style.label, self.productFrame, energy, energyAngular ) )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        tempInfo['productFrame'] = self.productFrame
        return( self.energyAngularSubform.processMultiGroup( style, tempInfo, indent ) )

    def spectrum( self, frame, **kwargs ) :
        """
        Returns an XYs3d instance representing *self*'s P(E',mu|E) in the requested frame.
        
        :param frame:               The frame the spectrum is returned in.
        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    A XYs3d instance.
        """

        if( not( isinstance( self.energyAngularSubform, XYs3d ) ) ) : raise TypeError( 'Form "%s" is not supported' % self.energyAngularSubform.moniker )

        if( frame == self.productFrame ) : return( self.energyAngularSubform.toPointwise_withLinearXYs( accuracy = 1e-3 ) )
        if frame == xDataEnumsModule.Frame.centerOfMass:
            raise ValueError( 'Conversion from frame "%s" to "%s" not supported' % ( self.productFrame, frame ) )

        xys3d = XYs3d( axes = defaultAxes( self.domainUnit ) )
        for energy in self.energyAngularSubform.domainGrid : xys3d.append( self.spectrumAtEnergy( energy, frame ) )

        return( xys3d )

    def spectrumAtEnergy( self, energyIn, frame ) :
        """
        Returns an XYs2d instance representing *self*'s P(E',mu|E=energyIn) at prjojectile energy *energyIn* and 
        in the requested frame.

        :param energyIn:            Energy of the projectile.
        :param frame:               The frame the spectrum is returned in.

        :return:                    A XYs2d instance.
        """

        class AngualarAtEnergyCOM :
            """
            The class needed by the function energyAngularSpectrumFromCOMSpectrumToLabAtEnergy to store data and evaluate the
            angular distribution in the center-of-mass.
            """

            def __init__( self, probability ) :

                self.probability = probability

            def probabilityCOM( self, energyPrimeCOM, muCOM ) :

                return( self.probability.evaluate( energyPrimeCOM ).evaluate( muCOM ) )

        _spectrumAtEnergy = self.energyAngularSubform.evaluate( energyIn ).toPointwiseLinear( lowerEps = 1e-7 ) # lowerEps should be user input and not hardwired.
        if( frame == self.productFrame ) : return( _spectrumAtEnergy )

        if frame == xDataEnumsModule.Frame.centerOfMass:
            raise ValueError( 'Conversion from frame "%s" to "%s" not supported' % ( self.productFrame, frame ) )

        energyProbability = [ [ probability.outerDomainValue, probability.integrate() ] for probability in _spectrumAtEnergy ]
        energyProbability = energyModule.XYs1d( energyProbability, axes = energyModule.defaultAxes( self.domainUnit ) )
        return( miscellaneousModule.energyAngularSpectrumFromCOMSpectrumToLabAtEnergy( self, energyIn, energyProbability, 
                AngualarAtEnergyCOM( _spectrumAtEnergy ), angularIsNormalized = False ) )

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
        subformElement = element[0]
        subformClass = {
                XYs3d.moniker: XYs3d,
                Regions3d.moniker: Regions3d,
                }.get( subformElement.tag )
        if subformClass is None: raise Exception( "encountered unknown energyAngular subform: %s" % subformElement.tag )
        subForm = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)
        energyAngular = cls( element.get( "label" ), element.get('productFrame'), subForm )
        xPath.pop()
        return energyAngular
