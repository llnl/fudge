# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This moudule contains the distribution suite class. This class stores the list for distribution forms for a product.

This module contains the following classes:
        
    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+ 
    | Component     | Class representing the distribution suite.                            |
    +---------------+-----------------------------------------------------------------------+
"""

from PoPs import IDs as IDsPoPsModule

from fudge import abstractClasses as abstractClassesModule
from fudge import styles as stylesModule
from xData import standards as standardsModule
from xData import enums as xDataEnumsModule

from . import angular as angularModule
from . import energy as energyModule
from . import energyAngular as energyAngularModule
from . import energyAngularMC as energyAngularMCModule
from . import angularEnergyMC as angularEnergyMCModule
from . import KalbachMann as KalbachMannModule
from . import angularEnergy as angularEnergyModule
from . import LLNL_angularEnergy as LLNL_angularEnergyModule
from . import uncorrelated as uncorrelatedModule
from . import Legendre as LegendreModule
from . import photonScattering as photonScatteringModule
from . import reference as referenceModule
from . import multiGroup as multiGroupModule
from . import unspecified as unspecifiedModule
from . import branching3d as branching3dModule

# probably missing stuff from photonScattering.py.

class Component( abstractClassesModule.Component ) :
    """
    Class representing the distribution suite. An instance of this class stores the list for distribution forms for a product.
    """

    moniker = 'distribution'

    def __init__( self ) :

        abstractClassesModule.Component.__init__( self, ( angularModule.Form, angularModule.TwoBody,
                KalbachMannModule.Form,
                energyAngularModule.Form, energyAngularMCModule.Form,
                angularEnergyModule.Form, angularEnergyMCModule.Form,
                LLNL_angularEnergyModule.LLNLAngularEnergyForm,
                uncorrelatedModule.Form, LegendreModule.Form, referenceModule.Form,
                referenceModule.CoulombPlusNuclearElastic, referenceModule.ThermalNeutronScatteringLaw,
                photonScatteringModule.CoherentPhotonScattering.Form, photonScatteringModule.IncoherentPhotonScattering.Form, photonScatteringModule.IncoherentBoundToFreePhotonScattering.Form,
                multiGroupModule.Form, unspecifiedModule.Form, branching3dModule.Form ) )

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn* for frame *frame*,

        :param energy_in:           Energy of the projectile.
        :param frame:               The frame to calculate the energy spectrum in.
        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    XYs1d instance for the energy spectrum.
        """

        styleLabel = kwargs.get('styleLabel', self.evaluated.label)
        form = self[styleLabel]
        if hasattr(form, 'energySpectrumAtEnergy'):
            if frame == xDataEnumsModule.Frame.centerOfMass:
                if form.productFrame == xDataEnumsModule.Frame.lab:
                    form = None
        else:
            form = None

        if form is not None:
            return form.energySpectrumAtEnergy(energyIn, frame, **kwargs)
        else:
            form = self[styleLabel]
            if hasattr(form, 'energySpectrumAtEnergy'):
                print('        WARNING: lab to center-of-mass translation not supported.')
            else:
                print('        WARNING: distribution "%s" does not have energySpectrumAtEnergy method.' % form.moniker)
            print('            %s' % self.toXLink())
            return energyModule.XYs1d(axes=energyModule.defaultAxes(form.domainUnit))

    def getSpectrumAtEnergy( self, energy ) :
        """This method is deprecated, use energySpectrumAtEnergy instead. Returns the energy spectrum for self at projectile energy."""

        return self.energySpectrumAtEnergy(energy, xDataEnumsModule.Frame.lab)

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.

        :return:        The product's average energy and momentum data.
        """

        form = style.findFormMatchingDerivedStyle( self )
        if( form is None ) : raise Exception( 'No matching style' )
        return( form.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def check( self, info ):
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning, enums
        warnings = []

        for form in self:
            if not isinstance(self.rootAncestor.styles[form.label], stylesModule.Evaluated):
                continue

            if info['interaction'] is enums.Interaction.atomic:
                if form.moniker not in (
                        photonScatteringModule.CoherentPhotonScattering.Form.moniker,
                        photonScatteringModule.IncoherentPhotonScattering.Form.moniker,
                        angularModule.TwoBody.moniker,
                        uncorrelatedModule.Form.moniker,
                        unspecifiedModule.Form.moniker):
                    warnings.append(warning.WrongDistributionComponent(form.moniker, '2-body'))

            elif info['isTwoBody']:
                if form.productFrame != xDataEnumsModule.Frame.centerOfMass:
                    warnings.append( warning.Wrong2BodyFrame( form ) )

                if form.moniker not in (angularModule.TwoBody.moniker,
                                        referenceModule.Form.moniker,
                                        referenceModule.CoulombPlusNuclearElastic.moniker,
                                        unspecifiedModule.Form.moniker):
                    warnings.append( warning.WrongDistributionComponent( form.moniker, '2-body' ) )
            else:
                if form.moniker in (angularModule.TwoBody.moniker,
                                    angularModule.Form.moniker,
                                    energyModule.Form.moniker):
                    warnings.append( warning.WrongDistributionComponent( form.moniker, 'N-body' ) )

            def checkSubform(subform, contextMessage):
                distributionErrors = []
                if hasattr(subform, 'domainMin') and (subform.domainMin, subform.domainMax) != info['crossSectionDomain']:
                    domain = (subform.domainMin, subform.domainMax)
                    # For gamma products, domainMin should be >= cross-section start, upper bounds should match.
                    if self.ancestor.pid == IDsPoPsModule.photon:
                        startRatio = subform.domainMin / info['crossSectionDomain'][0]
                        endRatio = subform.domainMax / info['crossSectionDomain'][1]
                        if (startRatio < 1-standardsModule.Floats.epsilon or endRatio < 1-standardsModule.Floats.epsilon
                                or endRatio > 1+standardsModule.Floats.epsilon):
                            distributionErrors.append( warning.Domain_mismatch(
                                    *(domain + info['crossSectionDomain']), obj=subform ) )
                    # For all other products, check lower and upper edges: only warn if they disagree by > eps
                    else:
                        for e1,e2 in zip(domain, info['crossSectionDomain']):
                            ratio = e1 / e2
                            if ratio < 1-standardsModule.Floats.epsilon or ratio > 1+standardsModule.Floats.epsilon:
                                distributionErrors.append(warning.Domain_mismatch(
                                        *(domain + info['crossSectionDomain']), obj=subform))
                                break
                if not hasattr(subform, 'check'):
                    distributionErrors.append(warning.NotImplemented(subform.moniker, subform))
                    if info['failOnException']:
                        raise NotImplementedError("Checking distribution form '%s'" % subform.moniker)
                else:
                    distributionErrors += subform.check( info )
                if distributionErrors:
                    warnings.append( warning.Context( contextMessage + " - %s:" % subform.moniker, distributionErrors) )

            if isinstance(form, uncorrelatedModule.Form):
                for subformName in ('angularSubform','energySubform'):
                    subform = getattr(form, subformName ).data
                    checkSubform( subform, 'uncorrelated - ' + subformName.replace('Subform','') )

            elif isinstance(form, KalbachMannModule.Form):
                checkSubform( form, form.moniker )

            elif isinstance(form, branching3dModule.Form):
                checkSubform(form, form.moniker)

            else:
                for subform in form.subforms:
                    checkSubform( subform, form.moniker )

        return warnings

    def diff( self, other, diffResults ) :
        """
        Adds information to *diffResults* if one of *self* or *other* has data and the other does not.

        :param diffResults:     List if diff results.
        """

        if( self.hasData( ) != other.hasData( ) ) :
            if( self.hasData( ) ) :
                diffResults.append( 'Distribution unspecified - 2', '', self.toXLink( ), other.toXLink( ) )
            else :
                diffResults.append( 'Distribution unspecified - 1', '', self.toXLink( ), other.toXLink( ) )

    def patch( self, other ) :

        pass

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides :py:meth:`ancestryModule.Ancestry.findEntity`.

        :param entityName:  Name is the desired entry in *self*.
        :param attribute:   Name of an attribute in *self*.
        :param value:       Value of the named attribute.

        :return:            The entry matching *entityName* or *attribute*/*value*.
        """

        if attribute is not None:
            for entity in self:
                if entity.moniker == entityName and getattr(entity,attribute) == value:
                    return entity
        else:
            for entity in self:
                if entity.moniker == entityName:
                    return entity
        return abstractClassesModule.Component.findEntity( self, entityName, attribute, value )

    def hasData( self ) :
        """
        Returns False if self's only has unspecified form; otherwise, returns True.
        """

        for form in self :
            if( not( isinstance( form, unspecifiedModule.Form ) ) ) : return( True )
        return( False )

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

        if( len( self ) > 0 ) :
            form = self[0]
#            if( form.productFrame == xDataEnumsModule.Frame.centerOfMass: return( 0.0 )
            if( hasattr( form, 'integrate' ) ) :
                return( form.integrate( reaction_suite, energyIn, energyOut = energyOut, muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder ) )
            else :
                print( 'missing integrate', type( form ) )

        return( 0.0 )

    def isSpecified(self):
        """
        Returns True if *self* has at least one distribution representation (i.e., form) and if the first representation 
        is not an :py:class:`unspecifiedModule.Form)` instance.

        :return:                Boolean.
        """

        if len(self) == 0:
            return False
        if isinstance(self[0], unspecifiedModule.Form):
            return False

        return True

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    ?
        """

        return( self.evaluated.toPointwise_withLinearXYs( **kwargs ) )
