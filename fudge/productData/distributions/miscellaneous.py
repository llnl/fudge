# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains various functions that do common calculations.

This module contains the following functions:
        
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | Function                                              | Description                                                           | 
    +=======================================================+=======================================================================+
    | energyAngularSpectrumFromCOMSpectrumToLabAtEnergy     | This function converts a center-of-mass (COM) P(E') and P(E',mu) into |
    |                                                       | a lab frame.                                                          |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | calculateDepositionEnergyFromEpP                      |                                                                       |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | calculateDepositionEnergyFromAngular_angularEnergy    |                                                                       |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | GaussQuadrature2                                      |                                                                       |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | GnG_adaptiveQuadrature                                |                                                                       |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | domainLimits                                          |                                                                       |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | integratePhi                                          |                                                                       |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
    | muPhiEvaluate                                         |                                                                       |
    +-------------------------------------------------------+-----------------------------------------------------------------------+
"""

import sys
import math

from xData import axes as axesModule
from xData import XYs1d as XYs1dModule

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from fudge.productData import averageProductEnergy as averageProductEnergyModule

def energyAngularSpectrumFromCOMSpectrumToLabAtEnergy( self, energyIn, energySpectrumAtEnergyCOM, energyAngualarAtEnergyCOM, angularIsNormalized = True ) :
    """
    This function converts a center-of-mass (COM) P(E') and P(E',mu) into a lab frame. 

    :param self:                            The distribution form from which energySpectrumAtEnergyCOM and energyAngualarAtEnergyCOM are derived.
    :param energySpectrumAtEnergyCOM:       P(E') in the center-of-mass frame.
    :param energyAngualarAtEnergyCOM:       P(E',mu) or P(mu|E') in the center-of-mass frame.
    :param angularIsNormalized:             If True, energyAngualarAtEnergyCOM is P(E',mu) otherwise it is P(mu|E').

    :return:                                An XYs2d instance representing P(E',mu) in the lab frame.
    """

    from . import energyAngular as energyAngularModule

    domainMin, domainMax = energySpectrumAtEnergyCOM.domainMin, energySpectrumAtEnergyCOM.domainMax

    def probabilityLab( energyPrimeLab, muLab, energyPrimeCOM, muCOM ) :
        """
        Calculates the probability in the lab frame for the specified inputs. For internal use only.

        :param energyPrimeLab:  The energy of the product in the lab frame.
        :param muLab:           The mu of the product in the lab frame.
        :energyPrimeCOM:        The energy of the product in the center-of-mass frame.
        :param muCOM:           The mu of the product in the center-of-mass frame.

        :return:                The list [muLab, P(muLab)].
        """

        if( energyPrimeCOM == 0.0 ) : energyPrimeCOM = domainMin + 1e-6 * ( domainMax - domainMin )
        PAtEnergyPrime = 1.0
        if( angularIsNormalized ) : PAtEnergyPrime = energySpectrumAtEnergyCOM.evaluate( energyPrimeCOM )
        return( [ muLab, math.sqrt( energyPrimeLab / energyPrimeCOM ) * PAtEnergyPrime * energyAngualarAtEnergyCOM.probabilityCOM( energyPrimeCOM, muCOM ) ] )

    numberOfEnergyPoints = 401
    numberOfAngularPoints = 201
    muMinMax = 0.999999

    energyUnit = self.domainUnit
    massUnit = energyUnit + '/c**2'

    reactionSuite = self.rootAncestor

    projectile = reactionSuite.PoPs[reactionSuite.projectile]
    projectileMass = projectile.getMass( massUnit )
    projectileZA = chemicalElementMiscPoPsModule.ZA( projectile )

    target = reactionSuite.PoPs[reactionSuite.target]
    targetMass = target.getMass( massUnit )
    targetZA = chemicalElementMiscPoPsModule.ZA( target )

    compoundMass = projectileMass + targetMass

    product = self.product
    productZA = chemicalElementMiscPoPsModule.ZA( product.particle )
    productMass = reactionSuite.PoPs[product.pid].getMass( massUnit )

    residualZA = projectileZA + targetZA - productZA
    residualMass = self.residualMass( reactionSuite.PoPs, residualZA, massUnit, compoundMass, product )

    if productMass == 0.0:
        raise Exception('Product "%s" is massless which is currently not supported.' % product.pid)

    if projectileMass == 0.0:
        COM_speed = energyIn / ( energyIn + targetMass )
    else:
        projectileSpeed = math.sqrt( 2.0 * energyIn / projectileMass )
        COM_speed = projectileMass / ( projectileMass + targetMass ) * projectileSpeed      # Speed of the center-of-mass.
    productEnergyWithCOM_speed = 0.5 * productMass * COM_speed * COM_speed

    spectrum = energyAngularModule.XYs2d( axes = energyAngularModule.defaultAxes( energyUnit ), outerDomainValue = energyIn )

    speedPrimeMax = math.sqrt( 2.0 * domainMax / productMass )
    if( domainMax >= productEnergyWithCOM_speed ) :
        energyPrimeLabMin = 0.0
        POfMu = [ [ -1.0, 0.0 ], [ 1.0, 0.0 ] ]
    else :
        if( domainMax <  ( 1.0 - muMinMax ) * productEnergyWithCOM_speed ) :        # Special treatment near threshold.
            numberOfEnergyPoints = 2
            speedPrimeMax = ( 1.0 - muMinMax ) * COM_speed
        energyPrimeLabMin = 0.5 * productMass * ( COM_speed - speedPrimeMax )**2
        POfMu = [ probabilityLab( energyPrimeLabMin, muMinMax, domainMax, -1.0 ), probabilityLab( energyPrimeLabMin, 1.0, domainMax, 1.0 ) ]

    spectrum.append( energyAngularModule.XYs1d( POfMu, axes = spectrum.axes, outerDomainValue = energyPrimeLabMin ) )

    if( numberOfEnergyPoints == 2 ) :
        POfMu = [ probabilityLab( productEnergyWithCOM_speed, muMinMax, domainMin, -1.0 ), probabilityLab( productEnergyWithCOM_speed, 1.0, domainMin, 1.0 ) ]
        spectrum.append( energyAngularModule.XYs1d( POfMu, axes = spectrum.axes, outerDomainValue = productEnergyWithCOM_speed ) )

    energyPrimeLabMax = 0.5 * productMass * ( COM_speed + speedPrimeMax )**2

    if( energyPrimeLabMin == 0.0 ) : energyPrimeLabMin = 1e-8 * energyPrimeLabMax
    fractionEnergyPrime = math.pow( energyPrimeLabMax / energyPrimeLabMin, 1.0 / ( numberOfEnergyPoints - 1 ) )
    energyPrimeLab = energyPrimeLabMin
    for energyPrimeLabIndex in range( 1, numberOfEnergyPoints - 1 ) :
        energyPrimeLab *= fractionEnergyPrime
        speedLab = math.sqrt( 2.0 * energyPrimeLab / productMass )

        muMin = ( productEnergyWithCOM_speed + energyPrimeLab - domainMax ) / ( 2.0 * math.sqrt( productEnergyWithCOM_speed * energyPrimeLab ) )
        if( muMin >= muMinMax ) : muMin = muMinMax
        muMin = max( -1.0, muMin )

        numberOfAngularPoints2 = numberOfAngularPoints
        if( muMin > 0.5 ) : numberOfAngularPoints2 //= 2
        if( muMin > 0.9 ) : numberOfAngularPoints2 //= 2

        sqrtEnergiesTimes2 = 2.0 * math.sqrt( energyPrimeLab * productEnergyWithCOM_speed )
        POfMu = []
        for muIndex in range( numberOfAngularPoints2 ) :
            fractionMu = float( muIndex ) / ( numberOfAngularPoints2 - 1 )
            muLab = ( 1.0 - fractionMu ) * muMin + fractionMu
            if( abs( muLab ) < 1e-8 ) : muLab = 0.0

            energyPrime = productEnergyWithCOM_speed + energyPrimeLab - muLab * sqrtEnergiesTimes2
            if( ( muIndex == numberOfAngularPoints2 - 1 ) ) :
                muLab = 1.0
                energyPrime = productEnergyWithCOM_speed + energyPrimeLab - ( muLab - 1e-6 ) * sqrtEnergiesTimes2

            muCOM = max( -1.0, min( 1.0, ( muLab * speedLab - COM_speed ) / math.sqrt( 2.0 * energyPrime / productMass ) ) )

            energyPrime = min( domainMax, max( energyPrime, domainMin ) )
            POfMu.append( probabilityLab( energyPrimeLab, muLab, energyPrime, muCOM ) )
        spectrum.append( energyAngularModule.XYs1d( POfMu, axes = spectrum.axes, outerDomainValue = energyPrimeLab ) )

    POfMu = [ probabilityLab( energyPrimeLabMax, muMinMax, domainMax, -1.0 ), probabilityLab( energyPrimeLabMax, 1.0, domainMax, 1.0 ) ]
    spectrum.append( energyAngularModule.XYs1d( POfMu, axes = spectrum.axes, outerDomainValue = energyPrimeLabMax ) )

    spectrum.normalize( insitu = True )
    return( spectrum )

def calculateDepositionEnergyFromEpP( E, EpP ) :
    """
`    Calculates the average product energy given P(E') in the lab frame.

    :param E:       Energy of the projectile.
    :param EpP:     List of [ E', P(E') ]
    """

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( 'a', 0, EpP.axes[0].unit )
    axes[1] = axesModule.Axis( 'b', 1, EpP.axes[1].unit )
    Ep = XYs1dModule.XYs1d( data = [ [ EpP[0][0], EpP[0][0] ], [ EpP[-1][0], EpP[-1][0] ] ], axes = axes )
    return EpP.integrateTwoFunctions(Ep)

def calculateDepositionEnergyFromAngular_angularEnergy( label, angular, energy, multiplicity, doingGammaMomentum = False, accuracy = 1e-6 ) :
    """
    This function calculates the average product energy (or average phton as product momentum if *doingGammaMomentum*
    is True) for angular distribution P(mu|E) and energy distribution P(E'|E,mu) in the lab frame. The distribution
    P(muy,E'|E) is the product P(mu|E) * P(E'|E,mu).

    :param label:               The label for the returned average product energy (momentum) form.
    :param angular:             The angular distribution P(mu|E).
    :param energy:              The energy distribution P(E'|E,mu).
    :param multiplicity:        The multiplicity of the product as a function of projectile energy.
    :param doingGammaMomentum:  If True, the product is a photon and its average momentum is to be calculated.
    :param accuracy:            The accuracy of the returned average energy (momentum) as a function of projectile energy.
    """

    energyUnit = energy.axes[0].unit
    energyPrimeUnit = energy.axes[2].unit
    momentumUnit = energyPrimeUnit + '/c'

    sqrtP6 = 0.77459666924148337704 / 2.      # sqrt( 0.6 ) / 2
    depEnergy = []
    for indexE, muEpPs in enumerate( energy ) :
        E = muEpPs.outerDomainValue
        I1MuP = angular[indexE]
        sum = 0.
        for indexMu, muEpP in enumerate( muEpPs ) :
            mu2, P2 = I1MuP[indexMu]
            mu2 = muEpP.outerDomainValue
            Ep2 = muEpP.integrateWithWeight_x( )
            if( indexMu != 0 ) :
                muMid = 0.5 * ( mu1 + mu2 )
                EpMid = muEpPs.interpolateAtValue( muMid, unitBase = True ).integrateWithWeight_x( )
                if( doingGammaMomentum ) :
                    dMu = sqrtP6 * ( mu2 - mu1 )
                    muG1 = muMid - dMu
                    muEpG1 = muG1 * muEpPs.interpolateAtValue( muG1, unitBase = True ).integrateWithWeight_x( )
                    muG2 = muMid + dMu
                    muEpG2 = muG2 * muEpPs.interpolateAtValue( muG2, unitBase = True ).integrateWithWeight_x( )
                    sum += ( mu2 - mu1 ) * ( 5. * ( muEpG1 + muEpG2 ) + 8. * muMid * EpMid )
                else :
                    sum += 3. * ( mu2 - mu1 ) * ( P1 * Ep1 + 2 * ( P1 + P2 ) * EpMid + P2 * Ep2 )   # 3 due to 18 instead of 6 below.
            P1 = P2
            mu1 = mu2
            Ep1 = Ep2
        depEnergy.append( [ E, multiplicity.evaluate( E ) * sum / 18. ] )

    if( doingGammaMomentum ) : return( depEnergy )

    axes = averageProductEnergyModule.defaultAxes( energyUnit = energyUnit )
    return( averageProductEnergyModule.XYs1d( data = depEnergy, axes = axes, label = label ) )

def GaussQuadrature2( function, parameters, a, b ) :
    """
    Returns the integral of *function* from *a* to *b*.

    :param function:    The function to integrate.
    :param parameters:  Paramters to pass to *function*.
    :param a:           Lower limit of the integral.
    :param b:           Upper limit of the integral.
    """

    if( a == b ) : return( 0. )
    xp, m, width = 0.57735026918962576451, 0.5 * ( a + b ), b - a                 # sqrt( 1. / 3. ), mid-point, width
    x1, x2 = m - 0.5 * width * xp, m + 0.5 * width * xp
    return( 0.5 * width * ( function( x1, parameters ) + function( x2, parameters ) ) )

def GnG_adaptiveQuadrature( function, a, b, parameters, quadrature, tolerance, maxEvaluations = 1000 ) :
    """
    This function uses adaptive quadrature to integrate *function* form *a* to *b*.

    :param function:            The function to integrate.
    :param a:                   The lower limit of the intrgral.
    :param b:                   The upper limit of the intrgral.
    :param parameters:          Parameters to pass to *function*.
    :param quadrature:          The integration quadrature to use.
    :param tolerance:           The requested accuracy of the integral.
    :param maxEvaluations:      The maximum number of evaluations to perform.
    """

    class QuadratureInfo :

        def __init__( self, function, parameters, quadrature, estimate, maxEvaluations ) :

            self.function = function
            self.parameters = parameters
            self.quadrature = quadrature
            self.estimate = estimate
            self.evaluations = 0
            self.totalEevaluations = 0
            self.maxEvaluations = maxEvaluations
            self.maxLevelReached = 0

        def GnG_adaptiveQuadrature_2( self, course, a, b, level ) :

            if( a == b ) : return( 0. )
            self.evaluations += 1
            if( self.evaluations > self.maxEvaluations ) : return( 0. )
            level += 1
            if( level > self.maxLevelReached ) : self.maxLevelReached = level
            m = 0.5 * ( a + b )
            l, r = self.quadrature( self.function, self.parameters, a, m ), self.quadrature( self.function, self.parameters, m, b )
            fine = l + r
            extrapolate = ( 16. * fine - course ) / 15.
            if( self.estimate + ( extrapolate - fine ) == self.estimate ) : return( fine )
            return( self.GnG_adaptiveQuadrature_2( l, a, m, level ) + self.GnG_adaptiveQuadrature_2( r, m, b, level ) )

        def __repr__( self ) :

            return( "evaluations = %d, totalEevaluations = %d, maxLevelReached = %d" % ( self.evaluations, self.totalEevaluations,
                self.maxLevelReached ) )

    if( a == b ) : return( 0., None )
    estimate = 0.
    for r in [ 0.9501, 0.2311, 0.6068, 0.4860, 0.8913 ] : estimate += function( a + ( b - a ) * r, parameters )
    estimate = 0.5 * ( estimate * 0.2 * ( b - a ) + quadrature( function, parameters, a, b ) )
    if( estimate == 0 ) : estimate = b - a
    if( tolerance < sys.float_info.epsilon ) : tolerance = sys.float_info.epsilon
    quadInfo = QuadratureInfo( function, parameters, quadrature, tolerance * estimate / sys.float_info.epsilon, maxEvaluations )
    course = quadrature( function, parameters, a, b )
    value = quadInfo.GnG_adaptiveQuadrature_2( course, a, b, 0 )
    r = value / estimate
    if( ( value != 0.0 ) and ( ( r < .1 ) or ( r > 10. ) ) ) :
        quadInfo.estimate = tolerance * value / sys.float_info.epsilon
        quadInfo.totalEevaluations = quadInfo.evaluations
        quadInfo.evaluations = 0
        quadInfo.maxLevelReached = 0
        value = quadInfo.GnG_adaptiveQuadrature_2( course, a, b, 0 )
    quadInfo.totalEevaluations += quadInfo.evaluations
    return( value, quadInfo )

def domainLimits( domain, domainMin, domainMax ) :
    """
    Returns the domain minimum and maximum given a domain.
    This is used to determine the limits of integration or, if the returned maximum is None, then the returned
    minimum is the point to evaluate a function at.
    Domain can be None, a number (int or float) or a list (or tuple) of length 1 or 2. All values of the 
    list must be a number.  If domain is None, domainMin and domainMax are returned. 
    If domain is a number or a list of length 1, the number and None are returned.
    Otherwise, domain is a list of length 2 with the first being the domain minimum and maximum.

    :param domain:      User specified domain limits.
    :param domainMin:   Minimum domain limit.
    :param domainMax:   Maximum domain limit.
    """

    if( domain is None ) : return( domainMin, domainMax )
    if( isinstance( domain, ( int, float ) ) ) : return( domain, None )        # Evaluate at domain.
    if( len( domain ) == 1 ) : return( domain[0], None )
    if( domain[1] == 0.0 ) : return( domain[0], None )

    return( domain[0], domain[1] )

def integratePhi( phiDomain ) :
    """
    Returns the integral of phi over the domain specified by *phiDomain*. The funtion :py:func:`domainLimits` is used to
    determine the phi domain from *phiDomain*. The full phi varies from 0 to 2 pi.

    :param phiDomain:   The phi domain.
    """

    twoPi = 2.0 * math.pi
    phiMin, phiMax = domainLimits( phiDomain, 0.0, twoPi )

    return( ( phiMax - phiMin ) / twoPi )

def muPhiEvaluate( muOut = None, phiOut = None ) :
    """
    Returns the integral over mu and phi over the limits specified by *muOut* and *phiOut*.

    :param muOut:   The mu domain.
    :param phiOut:  The phi domain.
    """

    muOutFactor = 1.0
    if( muOut is not None ) : muOutFactor = 0.5

    phiOutFactor = 1.0
    if( phiOut is not None ) : phiOutFactor = 1.0 / ( 2.0 * math.pi )

    return( muOutFactor * phiOutFactor )
