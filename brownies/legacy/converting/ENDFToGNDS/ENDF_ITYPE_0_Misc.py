# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Helper functions for reading ENDF ITYPE=0 data ('normal' reaction evaluations)
"""

epsilonExponent = -10
energyUnit = 'eV'

totalToken = 'total'
promptToken = 'prompt'
delayedToken = 'delayed'

import sys
import math

from pqu import PQU as PQUModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge.core.math import linearAlgebra as linearAlgebraModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import constant as constantModule
from xData import values as valuesModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule
from xData import link as linkModule
from xData import multiD_XYs as multiD_XYsModule
from xData import xDataArray as arrayModule
from xData import gridded as griddedModule
from xData import uncertainties as uncertaintiesModule

from PoPs import IDs as IDsPoPsModule
from PoPs import alias as PoPsAliasModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.fissionFragmentData import rate as rateModule
from PoPs.quantities import quantity as quantityModule

from fudge import enums as enumsModule
from fudge import physicalQuantity as physicalQuantityModule
from fudge import warning as warningModule
from fudge import externalFile as externalFileModule

import fudge.sums as sumsModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import \
    nuclearPlusInterference as nuclearPlusInterferenceModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import nuclearAmplitudeExpansion as nuclearAmplitudeExpansionModule

from fudge.resonances import resolved as resolvedModule
import fudge.resonances.resonances as resonancesModule
import fudge.resonances.scatteringRadius as scatteringRadiusModule
import fudge.resonances.resolved as resolvedResonanceModule
import fudge.resonances.unresolved as unresolvedResonanceModule
import fudge.resonances.common as commonResonanceModule

import fudge.covariances.covarianceMatrix as covarianceMatrixModule
import fudge.covariances.covarianceSuite as covarianceSuiteModule
import fudge.covariances.covarianceSection as covarianceSectionModule
import fudge.covariances.summed as covarianceSummedModule
import fudge.covariances.shortRangeSelfScalingVariance as shortRangeSelfScalingVarianceModule
import fudge.covariances.mixed as covarianceMixedModule
import fudge.covariances.modelParameters as covarianceModelParametersModule
import fudge.covariances.enums as covarianceEnumsModule

from fudge import outputChannel as outputChannelModule
from fudge.outputChannelData import Q as QModule
from fudge.outputChannelData.fissionFragmentData import fissionEnergyRelease as fissionEnergyReleaseModule
from fudge.outputChannelData.fissionFragmentData import delayedNeutron as delayedNeutronModule

import fudge.reactions.reaction as reactionModule
import fudge.reactions.production as productionModule
import fudge.reactionData.crossSection as crossSectionModule

import fudge.productData.multiplicity as multiplicityModule
from fudge.productData import averageProductEnergy as averageProductEnergyModule
import fudge.productData.distributions.unspecified as unspecifiedModule
import fudge.productData.distributions.angular as angularModule
import fudge.productData.distributions.energy as energyModule
import fudge.productData.distributions.uncorrelated as uncorrelatedModule
import fudge.productData.distributions.angularEnergy as angularEnergyModule
import fudge.productData.distributions.energyAngular as energyAngularModule
import fudge.productData.distributions.KalbachMann as KalbachMannModule
import fudge.productData.distributions.reference as referenceModule

from .. import endf_endl as endf_endlModule
from .. import toGNDSMisc as toGNDSMiscModule
from . import endfFileToGNDSMisc as endfFileToGNDSMiscModule

from . import metaStableData as metaStableDataModule

MTWithOnlyNeutronProducts = [2]
for MT in endf_endlModule.endfMTtoC_ProductLists:
    if (endf_endlModule.endfMTtoC_ProductLists[MT][IDsPoPsModule.neutron] > 0) or endf_endlModule.endfMTtoC_ProductLists[MT].isFission:
        MTWithOnlyNeutronProducts.append(MT)

frames = {1: xDataEnumsModule.Frame.lab, 2: xDataEnumsModule.Frame.centerOfMass}
FUDGE_EPS = endfFileToGNDSMiscModule.FUDGE_EPS
productNameToZA = {
    IDsPoPsModule.neutron: 1,
    'H1': 1001,
    'H2': 1002,
    'H3': 1003,
    'He3': 2003,
    'He4': 2004,
    IDsPoPsModule.photon: 0}
lightIsotopeNames = [IDsPoPsModule.neutron, 'H1', 'H2', 'H3', 'He3', 'He4']

crossSectionAxes = crossSectionModule.defaultAxes(energyUnit)
multiplicityAxes = multiplicityModule.defaultAxes(energyUnit)
averageProductEnergyAxes = averageProductEnergyModule.defaultAxes(energyUnit)
angularAxes = angularModule.defaultAxes(energyUnit)
energyAxes = energyModule.defaultAxes(energyUnit)
energyAngularAxes = energyAngularModule.defaultAxes(energyUnit)
angularEnergyAxes = angularEnergyModule.defaultAxes(energyUnit)
KalbachMann_f_Axes = KalbachMannModule.FSubform.defaultAxes(energyUnit)
KalbachMann_r_Axes = KalbachMannModule.RSubform.defaultAxes(energyUnit)
KalbachMann_a_Axes = KalbachMannModule.ASubform.defaultAxes(energyUnit)
fissionEnergyReleaseAxes = fissionEnergyReleaseModule.defaultAxes(energyUnit)


def particleZA(info, particleID):

    particle = info.reactionSuite.PoPs[particleID]
    return chemicalElementMiscPoPsModule.ZA(particle)


class MyIter:
    """Iterator that keeps track of line number."""

    def __init__(self, iterable):

        self.index = 0
        self.iterable = iter(iterable)
        self.length = len(iterable)

    def next(self):

        next_ = next(self.iterable)
        self.index += 1
        return next_


# Two useful shortcuts for reading ENDF data.
funkyF = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats


def funkyFI(a, logFile=sys.stderr) -> tuple[float, float, int, int, int, int]:   # read ENDF line with 2 floats and 4 ints

    return endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats(a, [2, 3, 4, 5], logFile=logFile)


# create some custom Exceptions:
class BadResonances(Exception): pass


class BadCovariance(Exception): pass


class DummyCrossSection:

    def __init__(self, domainMin, domainMax, unit):

        self.__domainMin = domainMin
        self.__domainMax = domainMax
        self.__unit = unit

    @property
    def domainMin(self) : return self.__domainMin

    @property
    def domainMax(self): return self.__domainMax

    @property
    def domainUnit(self): return self.__unit


def calculateZA(ZACompound, ZAOther, minus = True):
    """This function handles the removal (or addition) of ZAOther to ZACompound include natural compound (but not a natural other)."""

    if ZACompound % 1000 == 0: ZAOther = 1000 * (ZAOther // 1000)
    if minus: return ZACompound - ZAOther
    return ZACompound + ZAOther


def printAWR_mode(info, MT, MF, line, ZA, AWR, addToInfo=True, LIS=None):

    if addToInfo:
        info.ZA_massLineInfo.add(ZA, AWR, MT, MF, line, LIS=LIS)
    if info.AWR_mode is not None:
        info.AWR_mode.write("AWR_mode:: MT: %s: MF: %s: ZA: %s:  AWR: %s::\n" % (MT, MF, ZA, AWR))


def nudgeValue(value, sign):

    if value == 0: raise ValueError('FIXME')
    valueEpsilon = 10**(math.floor(math.log10(abs(value))) + epsilonExponent)
    return value + sign * valueEpsilon


def getCrossSectionForm(info, crossSectionRegions):

    axes = crossSectionAxes
    crossSectionRegions = [region for region in crossSectionRegions if len(region) > 1]
    if len(crossSectionRegions) == 1:         # Store as XYs1d.
        crossSectionForm = crossSectionModule.XYs1d(
            data=crossSectionRegions[0], label=info.style, axes=axes, interpolation=crossSectionRegions[0].interpolation)
    else:
        crossSectionForm = crossSectionModule.Regions1d(label=info.style, axes=axes)
        for region in crossSectionRegions:
            if len(region) > 1: crossSectionForm.append(region)
    return crossSectionForm


def getMultiplicity(multiplicity, EPrior, Ein):

    if isinstance(multiplicity, XYs1dModule.XYs1d):
        return multiplicity.evaluate(Ein)
    elif isinstance(multiplicity, regionsModule.Regions1d):
        for region in multiplicity:
            if Ein <= region.domainMax:
                if region.domainMax == EPrior == Ein: continue          # Next region is the one we want.
                return region.evaluate(Ein)
        return 0
    raise Exception('Unsupported multiplicity form "%s"' % multiplicity.moniker)


def uncorrelated(style, frame, angularSubform, energySubform):

    _angularSubform = uncorrelatedModule.AngularSubform(angularSubform)
    _energySubform = uncorrelatedModule.EnergySubform(energySubform)
    return uncorrelatedModule.Form(style, frame, _angularSubform, _energySubform)


def getMultiplicityPointwiseOrPieceWise(info, data, warningList):

    regions = []
    for region in data :
        if len(region) > 1: regions.append(region)
    if len(regions) == 1:
        multiplicity = multiplicityModule.XYs1d(data=regions[0], label=info.style, axes=multiplicityAxes,
                                                interpolation=regions[0].interpolation)
# BRB : fix me.
#    elif( ( len( regions ) == 2 ) and ( regions[0][-1] == regions[1][0] ) and
#            ( regions[0].axes[0].interpolation == regions[1].axes[0].interpolation ) ) :
#        xys = regions[1].copyDataToXYs( )      # This is a XYs1d data masquerading as Regions1d, convert to XYs1d.
#        xys[0][0] *= ( 1 + FUDGE_EPS )
#        xys = regions[0].copyDataToXYs( ) + xys
#        multiplicity = multiplicityModule.XYs1d( regions = xys, axes = axea )
    else:
        multiplicity = multiplicityModule.Regions1d(label=info.style, axes=multiplicityAxes)
        for region in regions:
            _region = region.copy()
            _region.axes = multiplicityAxes
            multiplicity.append(_region)
    return multiplicity


def getTotalOrPromptFission(info, MT, MTDatas, totalOrPrompt, warningList):

    MT456Data = MTDatas[MT][1]
    ZA, AWR, dummy, LNU, dummy, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats(
        MT456Data[0], logFile=info.logs)
    ZA = int(ZA)
    info.ZA_massLineInfo.add(ZA, AWR, MT, 1, 0)
    info.addMassAWR(ZA, AWR)
    LNU = int(LNU)
    info.logs.write('     %s fission neutron data: LNU = %d\n' % (totalOrPrompt, LNU))
    if LNU == 1:
        dataLine, polynomial = endfFileToGNDSMiscModule.getList(1, MT456Data, logFile=info.logs)
        domainMin, domainMax = 1e-5, 20e6        # BRB, these need to be set from data
        fissionMultiplicity = multiplicityModule.Polynomial1d(
            coefficients=polynomial['data'], label=info.style, axes=multiplicityAxes,
            domainMin=domainMin, domainMax=domainMax)
    else:
        dataLine, TAB1, multiplicityRegions = endfFileToGNDSMiscModule.getTAB1Regions(
            1, MT456Data, axes=multiplicityAxes, logFile=info.logs)
        fissionMultiplicity = getMultiplicityPointwiseOrPieceWise(info, multiplicityRegions, warningList)
    return fissionMultiplicity


def getDelayedFission(info, MT, MTDatas, warningList):

    info.logs.write('     Delayed fission neutron data (MT=455)')
    MT455Data = MTDatas[MT]
    MT455DataMF1 = MT455Data[1]
    ZA, AWR, LDG, LNU, dummy, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats(
        MT455DataMF1[0], logFile=info.logs)
    ZA = int(ZA)
    info.addMassAWR(ZA, AWR)
    info.ZA_massLineInfo.add(ZA, AWR, MT, 1, 0)
    LDG, LNU = int(LDG), int(LNU)
    info.logs.write(' LDG=%s LNU=%s' % (LDG, LNU))
    if LDG != 0: raise Exception("Only energy-independent delayed fission neutrons are supported")
    if LNU != 2: raise Exception("Only tables of delayed fission neutrons are supported")

    dataLine, decayRateData = endfFileToGNDSMiscModule.getList(1, MT455DataMF1, logFile=info.logs)
    NNF = int(decayRateData['NPL'])
    decayRates = decayRateData['data']

    dataLine, TAB1, multiplicityRegions = endfFileToGNDSMiscModule.getTAB1Regions(
        dataLine, MT455DataMF1, logFile=info.logs, axes=multiplicityAxes)
    info.totalDelayedMultiplicity = getMultiplicityPointwiseOrPieceWise(info, multiplicityRegions, warningList)

    interps = [region.interpolation for region in multiplicityRegions]
    if len(set(interps)) > 1:
        raise Exception("Currently only one interpolation flag is supported")
    nubarInterpolation = multiplicityRegions[0].interpolation

    if 5 in MT455Data:
        delayedNeutronEnergies, weights = readMF5(info, 455, MT455Data[5], warningList, delayNeutrons=True)
        weightsSum = XYs1dModule.XYs1d([], axes=weights[0].axes, interpolation=weights[0].interpolation)
        # Renormalize weights to sum to 1 (unless they are already normalized to ENDF precision)
        for weight in weights : weightsSum = weightsSum + weight
        if weightsSum.rangeMin < 0.999999 or weightsSum.rangeMax > 1.000001:
            for weight in weights:
                for i1, xy in enumerate(weight):
                    norm = weightsSum.evaluate(xy[0])
                    if norm == 0.0:
                        norm = 1.0
                    weight[i1] = [xy[0], xy[1] / norm]
    else:
        delayedNeutronEnergies = len(decayRates) * [None]
        if len(decayRates) > 1:
            warningList.append('More than one delayed fission neutron decay time but no MF = 5 data')

    if len(decayRates) != len(delayedNeutronEnergies):
        warningList.append("MF1 MT455 claims %d delayed neutron groups, but MF5 MT455 claims %d groups" %
                           (len(decayRates), len(delayedNeutronEnergies)))
        info.doRaise.append(warningList[-1])
        decayRates = []

    delayedNeutrons = []
    for i1, decayRate in enumerate(decayRates):
        energySubform = delayedNeutronEnergies[i1]
        if energySubform is not None:
            weight = weights[i1]
            weightsInterpolation = weight.interpolation
            if weightsInterpolation == xDataEnumsModule.Interpolation.flat:
                if len(weight) == 2 and nubarInterpolation == xDataEnumsModule.Interpolation.linlin:
                    weightsInterpolation = xDataEnumsModule.Interpolation.linlin
                    weight = XYs1dModule.XYs1d(data=[xy for xy in weight], axes=weight.axes)
            if weightsInterpolation != nubarInterpolation:
                raise Exception('For total nubar and energy interpolation differ which is currently not supported')
            if weightsInterpolation != xDataEnumsModule.Interpolation.linlin:
                raise Exception('For energy only "lin-lin" interpolation is supported: %s' % weightsInterpolation)
            totalDelayedM = info.totalDelayedMultiplicity
            if info.totalDelayedMultiplicity.domainMax > weight.domainMax:
                totalDelayedM = info.totalDelayedMultiplicity.domainSlice(domainMax=weight.domainMax)
            if weight.domainMax > totalDelayedM.domainMax:
                weight = weight.domainSlice(domainMax=totalDelayedM.domainMax)
            if weight.domainMin < totalDelayedM.domainMin:
                if len(weight) == 2:
                    if weight[0][1] == weight[1][1]: weight = weight.domainSlice(domainMin=totalDelayedM.domainMin)
            if totalDelayedM.domainMin > weight.domainMin:
                weight = weight.domainSlice(domainMin=totalDelayedM.domainMin)
            multiplicity = totalDelayedM * weight
            if isinstance(multiplicity, regionsModule.Regions1d):
                multiplicity = [region for region in multiplicity]
            else:
                multiplicity = [multiplicity]
            multiplicity = getMultiplicityPointwiseOrPieceWise(info, multiplicity, warningList)
            nuBar = multiplicity
        else :
            multiplicity = multiplicityModule.Unspecified(info.style)

        product = toGNDSMiscModule.newGNDSParticle(
            info, toGNDSMiscModule.getTypeNameGamma(info, 1), None, multiplicity=multiplicity)

        if energySubform is None:
            form = unspecifiedModule.Form(info.style, productFrame = xDataEnumsModule.Frame.lab)
        else:
            angularSubform = angularModule.Isotropic2d()
            form = uncorrelated(info.style, frames[1], angularSubform, energySubform)
        product.distribution.add(form)

        delayedNeutrons.append([decayRate, product])

    info.logs.write('\n')

    return delayedNeutrons


def getFissionEnergies(info, domainMin, domainMax, warningList):
    """
    For NPLY = 0 this data consists of pairs ( energy, standard deviation ).  The energies are:
    EFR: kinetic energy of all fission fragments
    ENP: energy of prompt fission neutrons
    END: energy of delayed fission neutrons
    EGP: energy of prompt fission gammas
    EGD: energy of delayed fission gammas
    EB: energy of delayed fission betas
    ENU: energy of fission neutrinos
    ER: EFR + ENP + END + EGP + EGD + EB, fission Q minus the neutrinos
    ET: ER + ENU, total fission Q
    For NPLY > 0 the fission energies are polynomials of degree NPLY in the incident energy,
    and the structure of the above table is repeated, once for each polynomial coefficient.
    """

    MF1Data = info.fissionEnergyReleaseData[1]
    dataLine = 0
    ZA, AWR, dummy, LFC, dummy, NFC = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats(
        MF1Data[dataLine], intIndices=[0, 5], logFile=info.logs)
    info.ZA_massLineInfo.add(ZA, AWR, 458, 1, 0)
    ZA = int(ZA)
    info.addMassAWR(ZA, AWR)

    dataLine += 1
    dummy, dummy, dummy, NPLY, N1, N2 = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats(
        MF1Data[dataLine], logFile=info.logs)
    if (N2 != (NPLY+1) * 9) or (N1 != N2 * 2): warningList.append("Inconsistent N1/N2/NPLY in section MF=1 MT=458!")
    nCoeffs = int(N2)  # total number of coefficients for all energy release components (each also has an uncertainty)

    dataLine += 1
    energies = endfFileToGNDSMiscModule.nFunkyFloatStringsToFloats(
        nCoeffs, dataLine, MF1Data, dimension=2, logFile=info.logs)
    dataLine += int(math.ceil(nCoeffs * 2 / 6.))

    FERterms = {}
    classes = (fissionEnergyReleaseModule.PromptProductKE, fissionEnergyReleaseModule.PromptNeutronKE,
               fissionEnergyReleaseModule.DelayedNeutronKE, fissionEnergyReleaseModule.PromptGammaEnergy,
               fissionEnergyReleaseModule.DelayedGammaEnergy, fissionEnergyReleaseModule.DelayedBetaEnergy,
               fissionEnergyReleaseModule.NeutrinoEnergy, fissionEnergyReleaseModule.NonNeutrinoEnergy,
               fissionEnergyReleaseModule.TotalEnergy)
    for idx, cls in enumerate(classes):
        coeffs, uncerts = zip(*energies[idx::9])
        poly1d = fissionEnergyReleaseModule.Polynomial1d(
            coefficients=coeffs, domainMin=domainMin, domainMax=domainMax, axes=fissionEnergyReleaseAxes,
            coefficientUncertainties=uncerts)

        FERterms[cls.moniker] = cls(data=poly1d)

    if LFC == 1:    # also have tabulated distributions for some or all energy release terms
        for NDCindex in range(NFC):
            dataLine, TAB1, regions = endfFileToGNDSMiscModule.getTAB1Regions(
                dataLine, MF1Data, allowInterpolation6=True, logFile=info.logs, axes=fissionEnergyReleaseAxes,
                cls=fissionEnergyReleaseModule.XYs1d)

            if len(regions) != 1:
                warningList.append("Tabular fission energy release with multiple regions not yet supported!")
                info.doRaise.append(warningList[-1])

            cls = classes[int(TAB1['L2'])-1]
            FERterms[cls.moniker] = cls(regions[0])

    return fissionEnergyReleaseModule.FissionEnergyRelease(label=info.style, **FERterms)


def angularLegendrePiecewiseToPointwiseIfPossible(piecewiseForm):

    if len(piecewiseForm) != 1: return piecewiseForm
    pointwise = angularModule.XYs2d(axes=angularAxes)
    for series in piecewiseForm[0]: pointwise.append(series)
    return pointwise


def angularLegendreToPointwiseOrPiecewiseLegendre(MT, angularData, warningList, MF, msg, subformPointwise=None):

    regions = []            # list of regions.
    i1 = 0
    priorEnergy = -1
    lists = angularData['Lists']
    try :
        for i2, interpolationFlag in angularData['interpolationInfo'] :
            if( interpolationFlag > 5 ) : raise Exception( 'Unsupported interpolation = %s for MF=%s, MT=%s' % ( interpolationFlag, MF, MT ) )
            interpolationQualifier, interpolationEin = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( interpolationFlag )
            region = [ interpolationQualifier, interpolationEin, [] ]
            if( i1 > 0 ) :
                if( lists[i1]['C2'] != lastRegion[0] ) : region[2].append( lastRegion )
            for i3 in range( i1, i2 ) :
                energy, coefficients = lists[i3]['C2'], [ 1.0 ] + lists[i3]['data']
                if( energy == priorEnergy ) :
                    if( coefficients == lastRegion[-1] ):
                        continue
                    if len(region[2]) > 1: regions.append( region )
                    region = [ interpolationQualifier, interpolationEin, [] ]
                lastRegion = ( energy, coefficients )
                region[2].append( lastRegion )
                priorEnergy = energy
            if len(region[2]) > 1: regions.append( region )
            i1 = i2
    except ValueError as err :
        warningList.append( 'ValueError occurred when constructing LegendreSeries: %s' % err )
        raise

    if( ( len( regions ) == 1 ) and ( subformPointwise is None ) ) :
        interpolationQualifier, interpolationEin, region = regions[0]
        subformLegendre = angularModule.XYs2d( axes = angularAxes, interpolation = interpolationEin,
                interpolationQualifier = interpolationQualifier )
        for i1, ( energy, coefficients ) in enumerate( region ) :
            subformLegendre.append( angularModule.Legendre( axes = angularAxes, coefficients = coefficients, outerDomainValue = energy ) )
    else :
        subformLegendre = angularModule.Regions2d( axes = angularAxes )
        for i1, regionInfo in enumerate( regions ) :
            interpolationQualifier, interpolationEin, region = regionInfo
            LegendreRegion = angularModule.XYs2d( axes = angularAxes, interpolation = interpolationEin,
                    interpolationQualifier = interpolationQualifier )
            for i2, ( energy, coefficients ) in enumerate( region ) :
                LegendreRegion.append( angularModule.Legendre( axes = angularAxes, coefficients = coefficients, outerDomainValue = energy ) )
            subformLegendre.append( LegendreRegion )
        if( subformPointwise is not None ) :
            if( isinstance( subformPointwise, angularModule.Regions2d ) ) :
                for region in subformPointwise:
                    subformLegendre.append( region )
            else :
                region = angularModule.XYs2d( interpolation = subformPointwise.interpolation,
                    interpolationQualifier = subformPointwise.interpolationQualifier, axes = angularAxes )
                for data in subformPointwise : region.append( data )
                subformLegendre.append( region )
    return( subformLegendre )


def convertNuclearPlusInterferenceDataToPiecewise(MT, angularData, warningList, MF, msg, identicalParticles):
    """
    Return three terms (nuclear + real/imaginary interference). These in turn contain
    Legendre expansions at various incident energies.
    """

    axes = angularAxes
    nuclear = angularModule.Regions2d( axes = axes )
    interferenceReal = angularModule.Regions2d( axes = axes )
    interferenceImaginary = angularModule.Regions2d( axes = axes )
    index, start, lastRegion = 0, 0, None
    lists = angularData['Lists']
    for end, interpolationFlag in angularData['interpolationInfo'] :
        interpolationQualifier, interpolationE_in = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( interpolationFlag )
        if( interpolationFlag > 6 ) : raise Exception( 'Unsupported interpolation = %s for MF=%s, MT=%s' % ( interpolationFlag, MF, MT ) )
        #interpolationE_in, interpolationCl, interpolationQualifier = endfFileToGNDSMiscModule.ENDFInterpolationToGNDSAxes3plusd( interpolationFlag )
        region_Nuc = angularModule.XYs2d( interpolation = interpolationE_in, axes = axes )
        region_IntReal = angularModule.XYs2d( interpolation = interpolationE_in, axes = axes )
        region_IntImaginary = angularModule.XYs2d( interpolation = interpolationE_in, axes = axes )
        priorEnergy = -1
        if( lastRegion is not None ) :
            if( lists[start]['C2'] != lastRegion ) : # ensure no gaps between regions
                region_Nuc.append( angularModule.Legendre( coefficients = nuclear_term, outerDomainValue = lastRegion, axes = axes ) )
                region_IntReal.append( angularModule.Legendre( coefficients = interference_termReal, outerDomainValue = lastRegion, axes = axes ) )
                region_IntImaginary.append( angularModule.Legendre( coefficients = interference_termImaginary, outerDomainValue = lastRegion, axes = axes ) )
        for idx in range( start, end ) :
            listNow = lists[idx]
            energy = listNow['C2']
            if( energy == priorEnergy ) :                           # This fixes a problem with some data having two same energy values.
                energy += FUDGE_EPS * energy
                warningList.append( 'same energies, second one being incremented for MT = %d, MF = %d, %s' % ( MT, MF, msg ) )
            priorEnergy = energy
            if( identicalParticles ) :
                splitPoint = listNow['N2'] + 1
            else :
                splitPoint = listNow['N2'] * 2 + 1
            nuclear_term = listNow['data'][:splitPoint]
            if( identicalParticles ) :  # nuclear_term only stores even-L coefficients. Add extra zeros for odd-L:
                newList = [nuclear_term[0]]
                for coef in nuclear_term[1:]: newList.extend([0,coef])
                nuclear_term = newList
            interference_term, interference_termReal, interference_termImaginary = listNow['data'][splitPoint:], [], []
            for jdx in range( 0, len( interference_term ), 2 ) :
                interference_termReal.append( interference_term[jdx] )
                interference_termImaginary.append( interference_term[jdx+1] )
            region_Nuc.append( angularModule.Legendre( coefficients = nuclear_term, outerDomainValue = energy, axes = axes ) )
            region_IntReal.append( angularModule.Legendre( coefficients = interference_termReal, outerDomainValue = energy, axes = axes ) )
            region_IntImaginary.append( angularModule.Legendre( coefficients = interference_termImaginary, outerDomainValue = energy, axes = axes ) )
        lastRegion = energy
        nuclear[index] = region_Nuc
        interferenceReal[index] = region_IntReal
        interferenceImaginary[index] = region_IntImaginary
        index += 1
        start = end
    nuclear = nuclearAmplitudeExpansionModule.NuclearTerm( angularLegendrePiecewiseToPointwiseIfPossible( nuclear ) )
    interferenceReal = nuclearAmplitudeExpansionModule.RealInterferenceTerm(
            angularLegendrePiecewiseToPointwiseIfPossible( interferenceReal ) )
    interferenceImaginary = nuclearAmplitudeExpansionModule.ImaginaryInterferenceTerm(
            angularLegendrePiecewiseToPointwiseIfPossible( interferenceImaginary ) )
    return( nuclear, interferenceReal, interferenceImaginary )


def convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1(MT, angularTAB1, warningList):

    angularData = angularTAB1['TAB2s']
    if( len( angularData ) == 1 ) :
        interplation, angularData = angularData[0]      # BRB: need to check interpolation?
        subform = angularModule.Regions2d(axes=angularAxes)
        xys2d = angularModule.XYs2d(axes=angularAxes)
        priorEnergy = -1.0
        for xys in angularData:
            if len(xys) > 1:
                raise NotImplementedError('help - need to support this')
            xys = xys[0]
            xys = angularModule.XYs1d(data=xys, axes=angularAxes, interpolation=xys.interpolation, outerDomainValue=xys.outerDomainValue)
            if xys.outerDomainValue == priorEnergy:
                if len(xys2d) > 0:
                    subform.append(xys2d)
                    xys2d = angularModule.XYs2d(axes=angularAxes)
            priorEnergy = xys.outerDomainValue
            xys2d.append(xys)
        if len(subform) > 0:
            if len(xys2d) > 0:
                subform.append(xys2d)
            if len(subform) == 1:
                subform = subform[0]
        else:
            subform = xys2d
    else :
        raise NotImplementedError( 'help - Regions2d tabulated angular not currently supported' )
    return subform

def convertAngularToPointwiseOrPiecewiseFromTAB2_List(MT, LANG, angularList, warningList):
    """
    Like convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1 except mu,P given as LISTs instead of TAB1s.
    """

    angularData = angularList['Lists']
    try :
        interpolation = { 12 : 2, 14 : 4 }[LANG]
    except :
        print( 'interpolation = LANG %d' % LANG )
        raise NotImplementedError( 'hell - what is this and how does it differ from 12' )
    interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS1d( interpolation )
    subform = angularModule.XYs2d( axes = angularAxes )
    e_in_Prior = -1
    for i1, listNow in enumerate( angularData ) :
        LANG_e_in = int( listNow['L1'] )
        if( LANG_e_in != LANG_e_in ) : raise NotImplementedError( 'hell - need to implement this' )
        e_in = listNow['C2']
        if( e_in == e_in_Prior ) : raise NotImplementedError( 'hell - need to implement this' )
        data = listNow['data']
        if( MT == 526 ) :        # Some MT 526 have data like '0.999998+0 2.046010+5 0.999998+0 2.870580+5' which this fixes.
            j1, j2 = None, None
            for j3 in range( 0, len( data ), 2 ) :
                if( j2 is not None ) :
                    if( data[j2] == data[j3] ) : data[j2] -= 0.5 * ( data[j2] - data[j1] )
                j1, j2 = j2, j3
        xys = angularModule.XYs1d( data = listNow['data'], dataForm = 'list', axes = angularAxes, interpolation = interpolation, outerDomainValue = e_in )
        subform.append( xys )
        e_in_Prior = e_in
    return( subform )


def toPointwiseOrPiecewiseEnergy(MT, TAB2):

    def getEpP( energyPrior, data ) :

        EpP = data[0]
        energy = float( EpP.outerDomainValue )
        if( len( data ) == 1 ) :
            EpP = energyModule.XYs1d( data = EpP, interpolation = EpP.interpolation, outerDomainValue = energy, axes = energyAxes )
        else :
            EpP = energyModule.Regions1d( outerDomainValue = energy, axes = energyAxes )
            for datum in data :
                EpP.append( energyModule.XYs1d( data = datum, interpolation = datum.interpolation, axes = energyAxes ) )
        return( energy, EpP )

    axes = energyAxes
    if( TAB2['NR'] == 1 ) :
        interpolation, data = TAB2['TAB2s'][0]
        interpolationQualifier, interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( interpolation )
        subform = energyModule.Regions2d( axes = axes )
        subformRegion = energyModule.XYs2d( axes = axes, interpolation = interpolation,
                interpolationQualifier = interpolationQualifier )
        energyPrior = -1
        for EpP in data :
            EIn = float( EpP[0].outerDomainValue )
            if( EIn == energyPrior ) :
                if( len( subformRegion ) > 1 ) : subform.append( subformRegion )
                subformRegion = energyModule.XYs2d( axes = axes, interpolation = interpolation, interpolationQualifier = interpolationQualifier )
                energyPrior = -1
            energyPrior, EpPp = getEpP( energyPrior, EpP )
            subformRegion.append( EpPp )
        if( len( subform ) == 0 ) :
            subform = subformRegion
        else :
            subform.append( subformRegion )
    else :
        subform = energyModule.Regions2d( axes = axes )
        TAB2s = TAB2['TAB2s']
        for i1, ( interpolation, TAB1s ) in enumerate( TAB2s ) :
            interpolationQualifier, interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( interpolation )
            region = energyModule.XYs2d( axes = axes, interpolation = interpolation, interpolationQualifier = interpolationQualifier )
            energyPrior = -1
            for EpP in TAB1s :
                energyPrior, EpPp = getEpP( energyPrior, EpP )
                region.append( EpPp )
            subform.append( region )

    return( subform )


def discreteOrPrimaryGamma(cls, energy, domainMin, domainMax):

    energySubform = cls(energy, domainMin, domainMax, axes = energyAxes)
    return energySubform


def translateENDFJpi(I, P):
    """
    endf uses weird convention for Jpi. Translate to simpler version
    """
    spin = abs(I)
    if I: parity = abs(I) / I
    else: parity = P or 1  # if (I,P)=(0,0) treat as '0+'
    return commonResonanceModule.Spin(spin), commonResonanceModule.Parity(parity)


def readMF2(info, MF2, warningList):
    """
    parse MF2 into resonances class (and sub-classes)
    """
    from xData import table as tableModule

    # store MT #s for all reactions that need to include resonance data:
    resonanceMTs = set()

    scatteringRadiusAxes = axesModule.Axes(2, labelsUnits={1: ('energy_in', 'eV'), 0: ('radius', 'fm')})

    def readResonanceSection(LRU, LRF, NRO, NAPS):
        """Helper function, read in resonance info for one energy range."""

        scatRadius = None
        if NRO != 0:  # energy-dependent scattering radius
            if NAPS == 2: raise BadResonances("NAPS=2 option not yet supported!")
            line1 = mf2.next()
            dum, dum, dum, dum, NR, NP = funkyFI(line1, logFile=info.logs)
            nLines = NR//3 + bool(NR % 3) + NP//3 + bool(NP % 3)
            data = [line1] + [mf2.next() for i in range(nLines)]
            dataLine, TAB1, regions = endfFileToGNDSMiscModule.getTAB1Regions(
                0, data,
                axes=axesModule.Axes(2,labelsUnits={1: ('energy_in', 'eV'), 0: ('radius', '10*fm')}),
                logFile=info.logs)
            if TAB1['NR'] != 1:
                raise NotImplementedError("multi-region scattering radius")
            data = regions[0]
            data = data.convertAxisToUnit(0, 'fm')
            scatRadius = scatteringRadiusModule.ScatteringRadius(data)

        if LRU == 0:  # scattering radius only. Note AP given in 10*fm
            SPI, AP, dum, dum, NLS, dum = funkyFI(mf2.next(), logFile=info.logs)
            info.particleSpins[info.target] = (commonResonanceModule.Spin(SPI), 0) # no parity information
            return AP*10, None

        elif LRU == 1 and LRF in (1, 2):   #SLBW or MLBW form
            SPI, AP, dum, dum, NLS, dum = funkyFI(mf2.next(), logFile=info.logs)
            info.particleSpins[info.target] = (commonResonanceModule.Spin(SPI), 0)
            resList = []
            negativeJs = False
            for lidx in range(NLS):
                AWRI_lineNumber = mf2.index
                AWRI, QX, L, LRX, tmp, NRS = funkyFI(mf2.next(), logFile=info.logs)
                info.ZA_massLineInfo.add(-1, AWRI, 151, 2, AWRI_lineNumber, column=0)
                if tmp != 6*NRS:
                    raise BadResonances("incorrect values in resonance section line %d" % mf2.index)
                for line in range(NRS):
                    e, j, gtot, gn, gg, gf = funkyF(mf2.next(), logFile=info.logs)
                    if j < 0: negativeJs = True
                    resList.append([e, L, j, gtot, gn, gg, gf])
            if negativeJs:
                raise BadResonances("Encountered negative J-values for SLBW/MLBW")

            table = tableModule.Table(
                columns=[
                    tableModule.ColumnHeader(0, name="energy", unit="eV"),
                    tableModule.ColumnHeader(1, name="L", unit=""),
                    tableModule.ColumnHeader(2, name="J", unit=""),
                    tableModule.ColumnHeader(3, name="totalWidth", unit="eV"),
                    tableModule.ColumnHeader(4, name="neutronWidth", unit="eV"),
                    tableModule.ColumnHeader(5, name="captureWidth", unit="eV"),
                    tableModule.ColumnHeader(6, name="fissionWidth", unit="eV")],
                data=sorted(resList, key=lambda res: res[0]))   # sort by energy only
            for column in ("totalWidth", "fissionWidth"):
                if not any(table.getColumn(column)):
                    table.removeColumn(column)

            if LRF == 1:
                approximation = resolvedModule.BreitWigner.Approximation.SingleLevel
            else:
                approximation = resolvedModule.BreitWigner.Approximation.MultiLevel

            BWresonances = resolvedResonanceModule.BreitWigner(
                info.style, approximation, commonResonanceModule.ResonanceParameters(table),
                scatteringRadius=scatRadius, calculateChannelRadius=not(NAPS))

            info.PoPsOverrides[BWresonances] = (AWRI, None)    # may need to override PoPs in resonance section

            return AP*10, BWresonances

        elif LRU == 1 and LRF == 3:     # Reich-Moore form, convert it to look like LRF=7
            from fudge.processing.resonances.reconstructResonances import getAllowedTotalSpins
            ENDFconversionFlags = ['LRF3']
            SPI, AP, LAD, dum, NLS, NLSC = funkyFI(mf2.next(), logFile=info.logs)
            if NLSC:
                ENDFconversionFlags.append('LvaluesNeededForConvergence=%d' % NLSC)
            info.particleSpins[info.target] = (commonResonanceModule.Spin(SPI), 0)  # store spin in GNDS particle list
            assert NRO == 0
            resDict = {}
            LdependentAP = {}
            haveFission = False
            for lidx in range(NLS):
                AWRI_lineNumber = mf2.index
                AWRI, APL, L, dum, tmp, NRS = funkyFI(mf2.next(), logFile=info.logs)
                info.ZA_massLineInfo.add(-1, AWRI, 151, 2, AWRI_lineNumber, column=0)
                if tmp != 6*NRS:
                    raise BadResonances("incorrect values in resonance section line %d" % mf2.index)
                if APL:
                    LdependentAP[L] = APL
                for line in range(NRS):
                    e, j, gn, gg, gf1, gf2 = funkyF(mf2.next(), logFile=info.logs)
                    if gf1 or gf2:
                        haveFission = True
                    channelSpin = abs(SPI + math.copysign(0.5, j))
                    j = abs(j)
                    if j == 0:
                        if SPI-L == 0.5:
                            channelSpin = SPI-0.5
                        elif SPI-L == -0.5:
                            channelSpin = SPI+0.5
                        else:
                            raise ValueError("Can't couple L=%s and S=%s to J=%s!" % (L,channelSpin,j))
                    resDict.setdefault(L, {}).setdefault(abs(j), {}).setdefault(channelSpin, []).append(
                        [e, gg, gn, gf1, gf2])

                for J in resDict[L]:
                    newChannelSpins = []
                    for channelSpin in resDict[L][J]:
                        if J not in getAllowedTotalSpins(L, channelSpin, useFactor2Trick=False):
                            # some evaluations don't indicate anything about channel spin, need to determine manually:
                            if J in getAllowedTotalSpins(L, SPI-0.5, useFactor2Trick=False):
                                newChannelSpins.append([SPI-0.5, channelSpin])
                                if 'ignoreChannelSpin' not in ENDFconversionFlags:
                                    ENDFconversionFlags.append('ignoreChannelSpin')
                            else:
                                raise ValueError("Can't couple L=%s and S=%s to J=%s! Allowed J values: %s" %
                                                 (L, channelSpin, J, getAllowedTotalSpins(L, channelSpin, False)))

                    for newChannelSpin, channelSpin in newChannelSpins:
                        resDict[L][J][newChannelSpin] = resDict[L][J][channelSpin]
                        del resDict[L][J][channelSpin]
            for L in range(max(NLS,NLSC)):
                if L not in resDict: resDict[L] = {}
                gsum = 0
                targetGSum = (2*L+1) * (2 * (2*SPI+1))
                for J in resDict[L]:
                    for S in resDict[L][J]:
                        gsum += (2*J+1)
                if gsum < targetGSum:  # add extra spin groups (with no resonances) for potential scattering:
                    for S in (SPI-0.5, SPI+0.5):
                        if S < 0: continue
                        J = abs(L - S)
                        jmax = L + S
                        while True:
                            jdict = resDict[L].setdefault(J,{})
                            if S not in jdict:
                                jdict[S] = []
                                gsum += (2*J+1)
                            J += 1
                            if J > jmax: break
                    if gsum != targetGSum:
                        raise ValueError("Method to fix missing gfactor failed! gsum should be %f, is %f" %
                                         (2*L+1,gsum / (2 * (2*SPI+1))))

            if LdependentAP:
                ENDFconversionFlags.append('explicitAPL')
            if AP == 0:
                ENDFconversionFlags.append('AP=0')
                AP = LdependentAP[0]

            resonanceReactions = commonResonanceModule.ResonanceReactions()
            mts = [102, 2]
            reactionLabels = {}
            if haveFission: mts.append(18)
            for MT in mts:
                gndsChannel, = [chan for chan in info.reactionSuite.reactions if chan.ENDF_MT == MT]
                reactionLabels[MT] = gndsChannel.label
                eliminated = (MT == 102)
                if MT == 2: ejectile = IDsPoPsModule.neutron
                elif MT == 102: ejectile = IDsPoPsModule.photon
                else: ejectile = None

                reactionLink = linkModule.Link(gndsChannel)
                resonanceReactions.add(
                    commonResonanceModule.ResonanceReaction(
                        label=gndsChannel.label, link=reactionLink, ejectile=ejectile,
                        eliminated=eliminated, scatteringRadius=None)
                )

            targetParity = 1    # parity is not usually stored in ENDF-6 (exception: LRF=7). For now assume it's positive
            if len(info.PoPs[info.target].nucleus.parity) > 0:
                targetParity = info.PoPs[info.target].nucleus.parity[0].value
            jdx = 0
            spinGroups = resolvedResonanceModule.SpinGroups()
            for L in sorted(resDict.keys()):
                for J in sorted(resDict[L].keys(), key=lambda val: abs(val)):
                    for channelSpin in sorted(resDict[L][J].keys()):

                        columnHeaders = [
                            tableModule.ColumnHeader(0, name="energy", unit="eV"),
                            tableModule.ColumnHeader(1, name=reactionLabels[102] + " width", unit="eV"),
                            tableModule.ColumnHeader(2, name=reactionLabels[2] + " width", unit="eV"),
                            tableModule.ColumnHeader(3, name="fission width_1", unit="eV"),
                            tableModule.ColumnHeader(4, name="fission width_2", unit="eV"),
                        ]

                        table = tableModule.Table(columns=columnHeaders[:], data=resDict[L][J][channelSpin])

                        channels = resolvedResonanceModule.Channels()
                        channels.add(resolvedResonanceModule.Channel(
                            "0", resonanceReactions[0].label, columnIndex=1, L=int(L),
                            channelSpin=commonResonanceModule.Spin(1)))    # capture
                        channels.add(resolvedResonanceModule.Channel(
                            "1", resonanceReactions[1].label, columnIndex=2, L=int(L),
                            channelSpin=commonResonanceModule.Spin(channelSpin)))
                        if L in LdependentAP and LdependentAP[L] != AP:
                            APL = constantModule.Constant1d(LdependentAP[L] * 10, domainMin=EL, domainMax=EH,
                                                            axes=scatteringRadiusAxes)
                            channels[-1].scatteringRadius = scatteringRadiusModule.ScatteringRadius(APL)
                            channels[-1].hardSphereRadius = scatteringRadiusModule.HardSphereRadius(APL)
                        if haveFission:
                            channels.add(resolvedResonanceModule.Channel(
                                "2", resonanceReactions[2].label, columnIndex=3, L=0,
                                channelSpin=commonResonanceModule.Spin(0)))
                            if len(table) > 0:
                                channels.add(resolvedResonanceModule.Channel(
                                    "3", resonanceReactions[2].label, columnIndex=4, L=0,
                                    channelSpin=commonResonanceModule.Spin(0)))
                                if not any(table.getColumn('fission width_2')):
                                    info.extraFissionWidths.append(table)
                        else:
                            table.removeColumn('fission width_2')
                            table.removeColumn('fission width_1')

                        parity = targetParity * (-1)**L
                        spinGroups.add(resolvedResonanceModule.SpinGroup(
                            str(jdx), commonResonanceModule.Spin(abs(J)), commonResonanceModule.Parity(parity),
                            channels, commonResonanceModule.ResonanceParameters(table)))
                        jdx += 1

            rmatrix = resolvedResonanceModule.RMatrix(
                info.style, resolvedResonanceModule.RMatrix.Approximation.ReichMoore, resonanceReactions, spinGroups,
                boundaryCondition=resolvedResonanceModule.BoundaryCondition.EliminateShiftFunction,
                calculateChannelRadius=not(NAPS), supportsAngularReconstruction=bool(LAD),
                relativisticKinematics=False, reducedWidthAmplitudes=False
            )

            info.PoPsOverrides[rmatrix] = (AWRI, None)

            if ENDFconversionFlags:
                info.ENDFconversionFlags.add(rmatrix, ",".join(ENDFconversionFlags))
            return AP*10, rmatrix

        elif LRU == 1 and LRF == 4:     # Adler-Adler, not currently supported
            raise BadResonances("Adler-Adler resonance formalism not yet supported!")

        elif LRU == 1 and LRF == 7:     # R-Matrix Limited
            dum, dum, IFG, KRM, NJS, KRL = funkyFI(mf2.next(), logFile=info.logs)
            assert IFG in (0, 1), f"IFG={IFG} not supported"
            if KRM == 3:
                approximation = resolvedResonanceModule.RMatrix.Approximation.ReichMoore
            elif KRM == 4:
                approximation = resolvedResonanceModule.RMatrix.Approximation.RMatrix
            else:
                raise BadResonances("R-Matrix with KRM=%d not yet implemented!\n" % KRM)

            dum, dum, NPP, dum, tmp1, tmp2 = funkyFI(mf2.next(), logFile=info.logs)
            if tmp1 != 12*NPP or tmp2 != 2*NPP:
                raise BadResonances("incorrect LRF7 header!")

            # some helper functions:
            def getOutgoingParticles(MT, targZA, projZA):
                reacStr = endf_endlModule.ENDF_MTZAEquation(projZA, targZA, MT)[1]
                outgoing = reacStr.split('->')[1].strip()
                pA, pB = outgoing.split()[::2]
                return pA, pB

            # ENDF R-Matrix starts by listing outgoing particle pairs
            # these are referred back to later on:
            resonanceReactions = commonResonanceModule.ResonanceReactions()
            SHFs = []
            for idx in range(NPP):
                MA_MB_lineNumber = mf2.index
                MA, MB, ZA, ZB, IA, IB = funkyF(mf2.next(), logFile=info.logs)
                Q, PNT, SHF, MT, PA, PB = funkyF(mf2.next(), logFile=info.logs)
                MT = int(MT)
                resonanceMTs.add(MT)

                if SHF == -1: SHF = 0   # format changed between ENDF-5 and -6, some evaluations still have old version
                SHFs.append(SHF)

                # identify the channel using ZA and MT:
                if MT in (18, 19):
                    pA = None
                else:
                    pA, pB = getOutgoingParticles(MT, int(ZAM), info.projectileZA)
                    # get target spin. In future, this should already be present in particle list
                    info.particleSpins[pA] = translateENDFJpi(IA, PA)
                    if MT != 102:   # spin/parity of 2nd particle are always 0 for capture in ENDF
                        info.particleSpins[pB] = translateENDFJpi(IB, PB)

                    # note: ZA and ZB in ENDF are charges, not ZA numbers. Compute ZA and add to particle mass dictionary:
                    ZA_A, ZA_B = endf_endlModule.ENDF_MTZAEquation(int(ZAM), info.projectileZA, MT)[0]
                    # ZA_A and ZA_B may not be in the same order as MA and MB, need to figure which ones go together
                    if ZA != ZB:
                        if (ZA_A // 1000) == ZB: ZA_A, ZA_B = ZA_B, ZA_A
                    else:   # two isotopes, so higher A goes with larger mass
                        if ZA_A > ZA_B: ZA_A, ZA_B = ZA_B, ZA_A
                        if MA > MB: MA, MB = MB, MA
                    info.addMassAWR(ZA_A, MA)
                    info.ZA_massLineInfo.add(ZA_A, MA, MT, 2, MA_MB_lineNumber)
                    info.addMassAWR(ZA_B, MB)
                    info.ZA_massLineInfo.add(ZA_B, MB, MT, 2, MA_MB_lineNumber)

                gndsChannel, = [chan for chan in info.reactionSuite.reactions if chan.ENDF_MT == MT]
                channelName = gndsChannel.label

                # Sanity check on PNT and SHF
                if PNT == 1 and MT in (18, 19, 102):
                    info.doRaise.append("Can't compute penetrability for MT%d!" % MT)
                elif PNT == -1 and MT not in (18, 19, 102):
                    info.doRaise.append("Unsupported: not computing penetrability for MT%d!" % MT)

                computeShift = {
                    -1: False,  # old (pre-2012) ENDF format
                    0: False,   # current ENDF format
                    1: True,    # For boundaryCondition = Given or NegativeOrbitalMomentum
                    # 2: False,   # SHF=2 now indicates Brune transform, but this is not yet handled
                }.get(SHF)
                if computeShift is None:
                    info.doRaise.append("Unexpected value for SHF: %d" % SHF)

                Qval = None
                eliminated = (KRM==3 and MT==102)
                if Q == 0 and gndsChannel.outputChannel.Q.evaluated.value >= 0:
                    pass
                elif gndsChannel.outputChannel.Q.evaluated.value != Q:
                    warningList.append("Resonance region Q-value doesn't match the rest of the evaluation")
                    originalQ = gndsChannel.outputChannel.Q.evaluated
                    newQ = QModule.Constant1d(Q, domainMin=originalQ.domainMin, domainMax=originalQ.domainMax,
                                              axes=originalQ.axes.copy(), label=info.style)
                    Qval = QModule.Component()
                    Qval.add(newQ)

                reactionLink = linkModule.Link(gndsChannel)
                resonanceReactions.add(
                    commonResonanceModule.ResonanceReaction(
                        label=channelName, link=reactionLink, ejectile=pA, Q=Qval, eliminated=eliminated))

            if 51 in resonanceMTs:
                resonanceMTs.add(4)
            if 600 in resonanceMTs:
                resonanceMTs.add(103)

            # next we have NJS spin groups, each containing channels and resonances
            spinGroups = resolvedResonanceModule.SpinGroups()
            radii = {}
            boundaryConditions = []
            for spinGroupIndex in range(NJS):
                # read scattering radius, binding, etc:
                AJ, PJ, KBK, KPS, tmp, NCH = funkyFI(mf2.next(), logFile=info.logs)
                if KPS != 0:
                    raise NotImplementedError("KPS != 0 in resonances LRF=7")
                if tmp != 6*NCH:
                    raise BadResonances("incorrect LRF7 header, line %d" % mf2.index)
                channels = resolvedResonanceModule.Channels()
                columnHeaders = [tableModule.ColumnHeader(0, name="energy", unit="eV")]
                channelNames = []
                widthUnit = "eV**(1/2)" if IFG==1 else "eV"
                for idx in range(NCH):
                    IPP, L, SCH, BND, APE, APT = funkyF(mf2.next(), logFile=info.logs)
                    thisChannel = resonanceReactions[int(IPP)-1]
                    BC = None
                    if BND not in (0, -L):
                        BC = BND
                    channels.add(resolvedResonanceModule.Channel(
                        str(idx), thisChannel.label, columnIndex=idx+1, L=int(L),
                        channelSpin=commonResonanceModule.Spin(SCH), boundaryConditionValue=BC))

                    channelName = "%s width" % thisChannel.label
                    jdx = 2
                    while True:
                        if channelName not in channelNames:
                            channelNames.append(channelName); break
                        channelName = '%s width_%d' % (thisChannel.label, jdx)
                        jdx += 1

                    radii.setdefault(thisChannel, []).append((channels[-1], APT*10, APE*10))
                    boundaryConditions.append((L, BND))
                    columnHeaders.append(tableModule.ColumnHeader(idx+1, name=channelName, unit=widthUnit))

                # resonances for this J:
                dum, dum, dum, NRS, tmp, NX = funkyFI(mf2.next(), logFile=info.logs)
                if tmp != 6*NX:
                    raise BadResonances("incorrect LRF7 header, line %d" % mf2.index)
                if NRS == 0: mf2.next()   # skip empty line
                resonances = []
                for i in range(NRS):
                    nlines = int(math.ceil((NCH+1)/6.0))  # Extra "1" is for the Eres column
                    vals = []
                    for j in range(nlines):
                        vals += funkyF(mf2.next(), logFile=info.logs)
                    resonances.append(vals[:NCH+1])
                table = tableModule.Table(columns=columnHeaders, data=resonances)
                # done with this spin group:
                J, pi = translateENDFJpi(AJ,PJ)
                spinGroups.add(resolvedResonanceModule.SpinGroup(
                    str(spinGroupIndex), J, pi, channels, commonResonanceModule.ResonanceParameters(table)))

                for kbk_idx in range(KBK):
                    from fudge.resonances import externalRMatrix as externalRMatrixModule
                    dum, dum, LCH, LBK, dum, dum = funkyFI(mf2.next(), logFile=info.logs)
                    if LBK == 0:
                        info.ENDFconversionFlags.add(channels[LCH-1], "LBK=0")
                        continue
                    if LBK == 2:
                        ED, EU, dum, dum, five, dum = funkyFI(mf2.next(), logFile=info.logs)
                        assert five == 5, "Malformed LBK=2 section in MF=2 line %d" % mf2.index
                        R0, R1, R2, S0, S1, dum = funkyF(mf2.next(), logFile=info.logs)
                        params = {}
                        for label, val, unit in (('constantExternalR', R0, ''),
                                                 ('linearExternalR', R1, '1/eV'),
                                                 ('quadraticExternalR', R2, '1/eV**2'),
                                                 ('constantLogarithmicCoefficient', S0, ''),
                                                 ('linearLogarithmicCoefficient', S1, '1/eV'),
                                                 ('singularityEnergyBelow', ED, 'eV'),
                                                 ('singularityEnergyAbove', EU, 'eV')):
                            params[label] = quantityModule.Double(label, val, unit)
                        external = externalRMatrixModule.SAMMY(**params)
                    elif LBK == 3:
                        ED, EU, dum, dum, three, dum = funkyFI(mf2.next(), logFile=info.logs)
                        assert three == 3, "Malformed LBK=3 section in MF=2 line %d" % mf2.index
                        R0, S0, GA, dum, dum, dum = funkyF(mf2.next(), logFile=info.logs)
                        params = {}
                        for label, val, unit in (('averageRadiationWidth', GA, 'eV'),
                                                 ('constantExternalR', R0, ''),
                                                 ('poleStrength', S0, ''),
                                                 ('singularityEnergyBelow', ED, 'eV'),
                                                 ('singularityEnergyAbove', EU, 'eV')):
                            if val != 0:
                                params[label] = quantityModule.Double(label, val, unit)
                        external = externalRMatrixModule.Froehner(**params)
                    else:
                        raise NotImplementedError("External R-Matrix with LBK=%d at MF=2 line %d" % (LBK, mf2.index))

                    channels[LCH-1].externalRMatrix = external

            # Use elastic scattering radius as 'default' value for entire resonance section
            from collections import Counter
            elastic, = [reac for reac in resonanceReactions if reac.link.link.ENDF_MT == 2]
            dum, trueRad, dum = zip(*radii[elastic])
            AP = Counter(trueRad).most_common(1)[0][0]

            # for each resonanceReaction, store default value for true radius and effective radius.
            # Only include them in specific channels if we need to override the default value.
            for resonanceReac in resonanceReactions:
                channel, trueRad, effRad = zip(*radii[resonanceReac])
                if not any(trueRad) and not any(effRad): continue     # ignore radii for (n,gamma)

                APTmostCommon = Counter(trueRad).most_common(1)[0][0]
                if APTmostCommon != AP:
                    resonanceReac.scatteringRadius = scatteringRadiusModule.ScatteringRadius(
                        constantModule.Constant1d(APTmostCommon, domainMin=EL, domainMax=EH, axes=scatteringRadiusAxes)
                    )

                APEmostCommon = Counter(effRad).most_common(1)[0][0]
                if APEmostCommon != APTmostCommon:
                    resonanceReac.hardSphereRadius = scatteringRadiusModule.HardSphereRadius(
                        constantModule.Constant1d(APEmostCommon, domainMin=EL, domainMax=EH, axes=scatteringRadiusAxes)

                    )

                for (channel,val1,val2) in radii[resonanceReac]:
                    if val1 != APTmostCommon:
                        channel.scatteringRadius = scatteringRadiusModule.ScatteringRadius(
                            constantModule.Constant1d(val1, domainMin=EL, domainMax=EH, axes=scatteringRadiusAxes)
                        )
                    if val2 != APEmostCommon:
                        channel.hardSphereRadius = scatteringRadiusModule.HardSphereRadius(
                            constantModule.Constant1d(val2, domainMin=EL, domainMax=EH, axes=scatteringRadiusAxes)
                        )

            # determine boundary condition. Most common: all == 0 or -L
            irregular = [bc for bc in boundaryConditions if bc[1] not in (0,-bc[0])]
            for tmp in irregular: boundaryConditions.remove(tmp)
            Ls, BCs = zip(*boundaryConditions)
            BCset = set(BCs)
            BCset_irregular = set()
            boundaryConditionValue = None
            if len(irregular) > 0:
                Ls_irregular, BCs_irregular = zip(*irregular)
                BCset_irregular = set(BCs_irregular)
            if BCset == {0}:    # BND = SHF, SAMMY convention
                boundaryCondition = resolvedResonanceModule.BoundaryCondition.EliminateShiftFunction
            elif all([shf == 2 for shf in SHFs]):
                boundaryCondition = resolvedResonanceModule.BoundaryCondition.Brune
                raise NotImplementedError("Brune transform (SHF=2) proposed but not yet implemented")
            elif all( [L == -BC for L, BC in boundaryConditions] ):
                boundaryCondition = resolvedResonanceModule.BoundaryCondition.NegativeOrbitalMomentum
            elif len(BCset) == 0 and len(BCset_irregular) == 1:  # non-zero constant boundary condition for all channels
                boundaryCondition = resolvedResonanceModule.BoundaryCondition.Given
                boundaryConditionValue = BCset_irregular.pop()
                for sg in spinGroups:
                    for chan in sg.channels:
                        chan.boundaryConditionValue = None
            else:
                raise Exception("Can't decipher boundary condition!")

            # end of spin groups. write RMatrix class:
            return AP, resolvedResonanceModule.RMatrix(
                info.style, approximation, resonanceReactions, spinGroups,
                relativisticKinematics=bool(KRL), reducedWidthAmplitudes=bool(IFG),
                calculateChannelRadius=False, boundaryCondition=boundaryCondition,
                boundaryConditionValue=boundaryConditionValue,
                supportsAngularReconstruction=True)

        elif LRU == 2:  # unresolved
            levelSpacingAxes = axesModule.Axes(2, labelsUnits={1: ('energy_in', 'eV'), 0: ('levelSpacing', 'eV')})
            widthAxes = axesModule.Axes(2, labelsUnits={1: ('energy_in', 'eV'), 0: ('average width', 'eV')})

            resonanceReactions = commonResonanceModule.ResonanceReactions()
            reactionLabels = {}
            for MT in (2, 102, 18):
                # ENDF-6 also supports a 'competitive' URR channel, but it is not associated with any reaction
                gndsChannel = [reac for reac in info.reactionSuite.reactions if reac.ENDF_MT == MT]
                if not gndsChannel: continue
                gndsChannel, = gndsChannel
                reactionLabels[MT] = gndsChannel.label
                reactionLink = linkModule.Link(gndsChannel)
                resonanceReactions.add(
                    commonResonanceModule.ResonanceReaction(label=gndsChannel.label, link=reactionLink, ejectile=None))

            L_list = unresolvedResonanceModule.Lsections()
            flags = []
            if LRF != 2:
                flags.append('LRF%d' % LRF)
            if LFW != int(info.reactionSuite.hasFission()):
                flags.append('LFW%d' % LFW)

            SPI, AP, LSSF, dum, NE, NLS = funkyFI(mf2.next(), logFile=info.logs)
            info.LSSF = bool(LSSF)
            if info.target not in info.particleSpins:
                info.particleSpins[info.target] = (commonResonanceModule.Spin(SPI), 0)
            if NRO == 0:
                scatRadius = scatteringRadiusModule.ScatteringRadius(
                    constantModule.Constant1d(AP*10, domainMin=EL, domainMax=EH, axes=scatteringRadiusAxes))

            if LFW == 0 and LRF == 1:   # 'Case A', see ENDF 2017 manual page 76
                NLS = NE
                for lidx in range(NLS):
                    J_list = unresolvedResonanceModule.Jsections()
                    AWRI_lineNumber = mf2.index
                    AWRI, dum, L, dum, tmp, NJS = funkyFI(mf2.next(), logFile=info.logs)
                    info.ZA_massLineInfo.add(-1, AWRI, MT, 2, AWRI_lineNumber, column=0)
                    if tmp != 6*NJS:
                        raise BadResonances("bad unresolved flag, line %d" % mf2.index)
                    for jidx in range(NJS):
                        D, AJ, AMUN, GNO, GG, dum = funkyF(mf2.next(), logFile=info.logs)

                        levelSpacing = unresolvedResonanceModule.LevelSpacing(
                            constantModule.Constant1d(D, domainMin=EL, domainMax=EH, axes=levelSpacingAxes))
                        widths = unresolvedResonanceModule.Widths()
                        widths.add(
                            unresolvedResonanceModule.Width(
                                '0',
                                reactionLabels.get(2),
                                degreesOfFreedom=AMUN,
                                data=constantModule.Constant1d(GNO, domainMin=EL, domainMax=EH, axes=widthAxes))
                        )
                        widths.add(
                            unresolvedResonanceModule.Width(
                                '1',
                                reactionLabels.get(102),
                                degreesOfFreedom=0,
                                data=constantModule.Constant1d(GG, domainMin=EL, domainMax=EH, axes=widthAxes))
                        )
                        J_list.add(unresolvedResonanceModule.Jsection(
                            str(jidx), commonResonanceModule.Spin(AJ), levelSpacing, widths))
                    L_list.add(unresolvedResonanceModule.Lsection(str(lidx), L, J_list))

            elif LFW == 1 and LRF == 1:  # 'Case B', only fission width is energy-dependent
                nlines = int(math.ceil(NE/6.0))
                energyList = []
                for i in range(nlines):
                    energyList += funkyF(mf2.next(), logFile = info.logs)
                energyList = energyList[:NE]

                for lidx in range(NLS):
                    J_list = unresolvedResonanceModule.Jsections()
                    AWRI_lineNumber = mf2.index
                    AWRI, dum, L, dum, NJS, dum = funkyFI(mf2.next(), logFile=info.logs)
                    info.ZA_massLineInfo.add(-1, AWRI, MT, 2, AWRI_lineNumber, column=0)
                    for jidx in range(NJS):
                        dum, dum, L, MUF, tmp, dum = funkyFI(mf2.next(), logFile=info.logs)
                        if tmp != NE+6:
                            raise BadResonances("Bad unresolved flag, line %d" % mf2.index)
                        D, AJ, AMUN, GNO, GG, dum = funkyF(mf2.next(), logFile=info.logs)
                        fissionWidth = []
                        for i in range(nlines):
                            fissionWidth += funkyF(mf2.next(), logFile=info.logs)
                        fissionWidth = fissionWidth[:NE]

                        levelSpacing = unresolvedResonanceModule.LevelSpacing(
                            constantModule.Constant1d(D, domainMin=EL, domainMax=EH, axes=levelSpacingAxes))
                        widths = unresolvedResonanceModule.Widths()
                        widths.add(
                            unresolvedResonanceModule.Width(
                                '0',
                                reactionLabels.get(2),
                                degreesOfFreedom=AMUN,
                                data=constantModule.Constant1d(GNO, domainMin=EL, domainMax=EH, axes=widthAxes))
                        )
                        widths.add(
                            unresolvedResonanceModule.Width(
                                '1',
                                reactionLabels.get(102),
                                degreesOfFreedom=0,
                                data=constantModule.Constant1d(GG, domainMin=EL, domainMax=EH, axes=widthAxes))
                        )

                        widths.add(
                            unresolvedResonanceModule.Width(
                                '2',
                                reactionLabels.get(18),
                                degreesOfFreedom=MUF,
                                data=XYs1dModule.XYs1d(list(zip(energyList, fissionWidth)), axes=widthAxes))
                        )
                        J_list.add(unresolvedResonanceModule.Jsection(
                            str(jidx), commonResonanceModule.Spin(AJ), levelSpacing, widths))
                    L_list.add(unresolvedResonanceModule.Lsection(str(lidx), L, J_list))

            elif LRF == 2:            # 'Case C', most common in ENDF-VII.1
                NLS = NE
                for Lidx in range(NLS):
                    J_list = unresolvedResonanceModule.Jsections()
                    AWRI_lineNumber = mf2.index
                    AWRI, dum, L, dum, NJS, dum = funkyFI(mf2.next(), logFile=info.logs)
                    info.ZA_massLineInfo.add(-1, AWRI, MT, 2, AWRI_lineNumber, column=0)
                    for jidx in range(NJS):
                        resList = []
                        AJ, dum, INT, dum, tmp, NE = funkyFI(mf2.next(), logFile=info.logs)
                        interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS1d(INT)
                        if tmp != 6*NE+6:
                            raise BadResonances("bad unresolved flag, line %d" % mf2.index)

                        dum, dum, AMUX, AMUN, AMUG, AMUF = funkyF(mf2.next(), logFile=info.logs)
                        for i in range(NE):
                            resList.append(funkyF(mf2.next(), logFile=info.logs))

                        # temporarily store in a table for convenience:
                        table = tableModule.Table(columns= [
                            tableModule.ColumnHeader(0, name="energy", unit="eV"),
                            tableModule.ColumnHeader(1, name="levelSpacing", unit="eV"),
                            tableModule.ColumnHeader(2, name="competitive", unit="eV"),
                            tableModule.ColumnHeader(3, name="elastic", unit="eV"),
                            tableModule.ColumnHeader(4, name="capture", unit="eV"),
                            tableModule.ColumnHeader(5, name="fission", unit="eV")],
                            data=resList)

                        elist = table.getColumn('energy')
                        if len(elist) != len(set(elist)):   # duplicates detected!
                            import collections
                            for energy, count in collections.Counter(elist).most_common():
                                if count == 1: break
                                rows = [row for row in table.data if row[0]==energy]
                                firstMatch = table.data.index(rows[0])
                                for row in rows[1:]:
                                    if row != rows[0]:
                                        raise NotImplementedError("Discontinuity in unresolved widths")
                                    index = table.data.index(rows[0], firstMatch+1)
                                    table.data.pop(index)
                                    warningList.append("removing duplicate energy from unresolved parameters")
                            elist = table.getColumn('energy')

                        levelSpacing = unresolvedResonanceModule.LevelSpacing(
                            XYs1dModule.XYs1d(list(zip(elist,table.getColumn('levelSpacing'))),
                                axes=levelSpacingAxes, interpolation=interpolation))

                        widths = unresolvedResonanceModule.Widths()
                        label = 0
                        for MT, width, DOF in ((2,'elastic',AMUN),(102,'capture',AMUG),
                                             (18,'fission',AMUF),(-1,'competitive',AMUX)):
                            column = table.getColumn(width)
                            if any(column) or DOF:  # keep all channels with DOF > 0 (even if widths are all 0)
                                if 0 in column and interpolation in (
                                        xDataEnumsModule.Interpolation.loglin, xDataEnumsModule.Interpolation.loglog):
                                    zidx = len(column) - column[::-1].index(0) + 1  # index of first non-zero value
                                    if zidx -1 != column.count(0):
                                        raise BadResonances("In URR: L=%d, J=%s, channel=%s widths fluctuate 0 to nonzero and back"
                                                            % (L, AJ, width))
                                    if zidx == len(column) or not any(column):
                                        data = XYs1dModule.XYs1d(list(zip(elist, column)), axes=widthAxes,
                                                                 interpolation=xDataEnumsModule.Interpolation.linlin)
                                    else:   # break up into regions
                                        data = regionsModule.Regions1d(axes=widthAxes)
                                        data.append(
                                            XYs1dModule.XYs1d(list(zip(elist[:zidx],column[:zidx])),
                                                              interpolation=xDataEnumsModule.Interpolation.linlin,
                                                              axes=widthAxes ) )
                                        data.append(
                                            XYs1dModule.XYs1d(list(zip(elist[zidx-1:],column[zidx-1:])),
                                                             interpolation=interpolation, axes=widthAxes))
                                else:
                                    data = XYs1dModule.XYs1d(list(zip(elist, column)),
                                                             axes=widthAxes, interpolation=interpolation)
                                reactionLabel = reactionLabels.get(MT, width)
                                widths.add(unresolvedResonanceModule.Width(str(label), reactionLabel, DOF, data))
                                label += 1

                        J_list.add(unresolvedResonanceModule.Jsection(
                            str(jidx), commonResonanceModule.Spin(AJ), levelSpacing, widths))
                    L_list.add(unresolvedResonanceModule.Lsection(str(Lidx), L, J_list))

            urr = unresolvedResonanceModule.TabulatedWidths(
                info.style, 'SingleLevelBreitWigner', resonanceReactions, L_list,
                scatteringRadius=scatRadius, useForSelfShieldingOnly=info.LSSF)

            info.PoPsOverrides[urr] = (AWRI, (commonResonanceModule.Spin(SPI), None))

            if flags:
                info.ENDFconversionFlags.add(urr, ",".join(flags))
            return AP * 10, urr

        else:
            info.logs.write("Unexpected LRU=%d, LRF=%d encountered\n" % (LRU, LRF))

    # end of helper functions.
    # now read MF2 data:
    mf2 = MyIter(MF2) # mf2.next() to get each line

    scatteringRadii = {}
    energyBounds = set()
    resolvedList = []
    unresolvedList = []

    # read MF2 header:
    AWR_lineNumber = mf2.index
    ZAM, AWR, dum, dum, NIS, dum = funkyFI(mf2.next(), logFile=info.logs)
    ZAM = int(ZAM)
    printAWR_mode(info, 151, 2, AWR_lineNumber, ZAM, AWR)
    info.addMassAWR(ZAM, AWR)
    if NIS != 1: info.logs.write("careful, more than one isotope in MF2!")
    ZAI, ABN, dum, LFW, NER, dum = funkyFI(mf2.next(), logFile=info.logs)

    for erange in range(NER):
        # each energy range
        EL, EH, LRU, LRF, NRO, NAPS = funkyFI(mf2.next(), logFile=info.logs)
        radius, resonanceSection = readResonanceSection(LRU, LRF, NRO, NAPS)
        scatteringRadii[LRU] = radius
        energyBounds.update([EL, EH])
        if resonanceSection is not None:
            resonanceMTs.update([1, 2, 3, 18, 19, 102])
        if LRU == 1:
            resolvedList.append((resonanceSection, EL, EH))
        elif LRU == 2:
            unresolvedList.append((resonanceSection, EL, EH))
    if not resolvedList: resolved = None
    elif len(resolvedList) == 1:
        form, domainMin, domainMax = resolvedList[0]
        resolved = resolvedResonanceModule.Resolved(domainMin, domainMax, domainUnit='eV')
        resolved.add(form)
    else:
        warningList.append("multiple resolved energy intervals are deprecated!")
        form = commonResonanceModule.EnergyIntervals(info.style)
        idx = 0
        for resonanceSection, EL, EH in resolvedList:
            interval = commonResonanceModule.EnergyInterval(idx, resonanceSection, EL, EH, domainUnit='eV')
            form.append(interval)
            idx += 1
        resolved = resolvedResonanceModule.Resolved(
            domainMin=resolvedList[0][1], domainMax=resolvedList[-1][2], domainUnit='eV')
        resolved.add(form)
    if not unresolvedList: unresolved = None
    elif len(unresolvedList) == 1:
        form, domainMin, domainMax = unresolvedList[0]
        unresolved = unresolvedResonanceModule.Unresolved(domainMin, domainMax, domainUnit='eV')
        unresolved.add(form)
        # reconstructCrossSection=not info.LSSF
    else:
        raise BadResonances("multiple unresolved regions not supported")

    if mf2.index != mf2.length:
        warningList.append("Not all resonance data converted!")
        info.doRaise.append(warningList[-1])

    for LRU in range(3):
        if LRU in scatteringRadii:
            AP = scatteringRadii[LRU]
            break
    scatteringRadius = scatteringRadiusModule.ScatteringRadius(
        constantModule.Constant1d(
            AP, domainMin=min(energyBounds), domainMax=max(energyBounds),
            axes=axesModule.Axes(2, labelsUnits = {1: ('energy_in', 'eV'), 0: ('radius', 'fm')}),
            label=info.style))
    if 2 in scatteringRadii:    # unresolved radius may be redundant
        if (1 not in scatteringRadii) or (scatteringRadii[2] == scatteringRadii[1]):
            if not isinstance(unresolved.evaluated.scatteringRadius.evaluated, XYs1dModule.XYs1d):
                unresolved.evaluated.scatteringRadius = None

    resonances = resonancesModule.Resonances(scatteringRadius, resolved=resolved, unresolved=unresolved)
    return resonances, sorted(resonanceMTs)


def readMF3( info, MT, MF3Data, warningList ) :

    ZA, AWR, dum, dum, dum, dum = funkyFI( MF3Data[0], info.logs )
    ZA = int( ZA )
    info.ZA_massLineInfo.add(ZA, AWR, MT, 3, 0)
    info.addMassAWR( ZA, AWR )
    dataLine, TAB1, crossSectionRegions = endfFileToGNDSMiscModule.getTAB1Regions( 1, MF3Data, allowInterpolation6 = True,
            logFile = info.logs, axes = crossSectionAxes, cls = crossSectionModule.XYs1d )
    QM, QI, LR = TAB1['C1'], TAB1['C2'], int( TAB1['L2'] )
    breakupProducts = None
    if(   LR == 0 ) :
        pass
    elif( LR in [ 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36 ] ) :
        info.logs.write( ' : MF=3, LR=%s' % LR )
        breakupProducts, productCounts = {}, endf_endlModule.endfMTtoC_ProductLists[LR].productCounts
        for product in productCounts :
            if( productCounts[product] != 0 ) : breakupProducts[product] = productCounts[product]
        breakupProducts[info.projectile] -= 1
        if( breakupProducts[info.projectile] == 0 ) : del breakupProducts[info.projectile]
    elif( LR == 31 ) :
        warningList.append( 'Invalid LR = %s for MT = %s is being ignored' % ( LR, MT ) )
    elif( LR == 1 ) :
        info.logs.write( ' : MF=3, LR=1' )
        pass
    elif( LR in [ 39, 40 ] ) :
        if( LR == 40 ) :
            warningList.append( 'LR = %s for MT = %s is being ignored' % ( LR, MT ) )
        else :
            if( MT != 5 ) :
                warningList.append( "Breakup LR = %s is not supported: MT = %s" % ( LR, MT ) )
                raise NotImplementedError( "Breakup LR = %s is not supported: MT = %s" % ( LR, MT ) )
    else :
        raise Exception( "Invalid breakup flag LR %s: MT = %d" % ( LR, MT ) )

    crossSection = getCrossSectionForm( info, crossSectionRegions )

    return( QM, QI, crossSection, LR, breakupProducts )


def readMF4(info, product, MT, MF4Data, formClass, warningList):

    if MT not in MTWithOnlyNeutronProducts: info.MF4ForNonNeutrons.append(MT)
    ZA, AWR, LVT, LTT, dummy, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF4Data[0], logFile = info.logs )
    ZA = int( ZA )
    printAWR_mode(info, MT, 4, 0, ZA, AWR)
    info.addMassAWR( ZA, AWR )
    LVT = int( LVT )                # 1: transformation matrix given. Must be 0 for endf/b6 format but not older formats.
    LTT = int( LTT )                # 0: isotropic, 1: Legendre, 2: table, 3: Legendre for low E and table for high E.

    dummy, AWR_, LI, LCT, NK, NM = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF4Data[1], logFile = info.logs )
    if AWR != AWR_:
        printAWR_mode(info, MT, 4, 1, ZA, AWR)
    else:
        info.ZA_massLineInfo.add(ZA, AWR_, MT, 4, 1)
    LI = int( LI )                  # if 1, gammas isotropic
    LCT = int( LCT )                # 1 for lab frame, 2 for center of mass
    NK = int( NK )                  # number of entries in transformation matrix
    NM = int( NM )                  # maximum Legendre order
    if( ( LCT != 2 ) and ( formClass == angularModule.TwoBody ) ):
        raise ValueError( "Discrete two-body must be in the center-of-mass frame: LCT = %d MT = %d." % ( LCT, MT ) )

    firstDataLine = 2
    if( LVT != 0 ) :
        warningList.append( 'MF = 4, MT = 2 contains obsolete matrix used to transform Legendre coefficients between frames.' )
        firstDataLine += ( NK + 5 ) // 6

    info.logs.write( ' : MF=4, LTT = %s' % LTT )
    if( LTT == 0 ) :                # Purely isotropic angular distribution.
        subform = angularModule.Isotropic2d( )
    elif( LTT == 1 ) :              # Legendre polynomial coefficient
        nextDataLine, angularData = endfFileToGNDSMiscModule.getTAB2_Lists( firstDataLine, MF4Data, logFile = info.logs )
        subform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 4, 'LTT = 1' )
    elif( LTT == 2 ) :              # Tabulated probability distribution
        nextDataLine, angularTable = endfFileToGNDSMiscModule.getTAB2_TAB1s( firstDataLine, MF4Data, logFile = info.logs )
        subform = convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1( MT, angularTable, warningList )
    elif( LTT == 3 ) :              # Mixed Legendre and Tabulated probability distribution
        nextDataLine, angularData = endfFileToGNDSMiscModule.getTAB2_Lists( firstDataLine, MF4Data, logFile = info.logs )
        nextDataLine, angularTable = endfFileToGNDSMiscModule.getTAB2_TAB1s( nextDataLine, MF4Data, logFile = info.logs )
        subformPointwise = convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1( MT, angularTable, warningList )
        subform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 4, 'LTT = 3', subformPointwise )
    else:
        raise ValueError("Encountered unknown LTT=%d in MF4" % LTT)

    if( formClass is None ) : return( subform )
    form = formClass( info.style, frames[LCT], subform )
    product.distribution.add( form )
    return( form )


def readMF5(info, MT, MF5Data, warningList, delayNeutrons=False, product=None):

    ZA, AWR, dummy, dummy, NK, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF5Data[0], logFile = info.logs )
    ZA = int( ZA )
    printAWR_mode( info, MT, 5, 0, ZA, AWR )
    info.addMassAWR( ZA, AWR )
    NK = int( NK )                 # number of partial energy distributions.
    dataLine = 1
    energySubforms, weights = [], []
    for k in range( NK ) :
        dataLine, productData = endfFileToGNDSMiscModule.getTAB1( dataLine, MF5Data, logFile = info.logs )
        if( productData['NR'] > 1 ) :
            oldInterpolations = productData['interpolationInfo']
            if( oldInterpolations == [ [ 3, 1 ], [ 4, 2 ], [ 5, 1 ] ] ) :       # This is a kludge for about 5 data sets, but does a good job.
                productData['NR'] = 1
                productData['data'].insert( 1, [ ( 1. - FUDGE_EPS ) * productData['data'][1][0], productData['data'][0][1] ] )
                productData['interpolationInfo'] = [ [ len( productData['data'] ), 2 ] ]
            else :
                raise ValueError( "Currently only one interpolation flag is supported" )
        LF = int( productData['L2'] )           # Energy distribution law
        xPrior, addWarning = None, True
        if( productData['data'][0][0] == productData['data'][1][0] ) : del productData['data'][0]
        weight = None
        if( ( NK > 1 ) and not( delayNeutrons ) ) :
            if( productData['NR'] != 1 ) : raise ValueError( "Currently only one interpolation flag is supported for weight for MF=5, MT=%s" % MT )
            x1, y1 = productData['data'][0]
            discontinuities = []
            if( productData['interpolationInfo'][0][1] != 1 ) :
                for i1, xy in enumerate( productData['data'][1:] ) :                # Check if data can be treated as flat interpolation.
                    x2, y2 = xy
                    if( x1 != x2 ) :
                        if( y1 != y2 ) : raise ValueError( 'Weight data must be convertible to flat interpolation' )
                    else :
                        discontinuities.insert( 0, i1 )
                    x1, y1 = x2, y2
                for discontinuity in discontinuities : del productData['data'][discontinuity]
            interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS1d( 1 )    # Flat interpolation.
            axes = axesModule.Axes(2, labelsUnits = { 1 : ( 'energy_in' , 'eV' ), 0 : ( 'weight' , '' ) } )
            weight = XYs1dModule.XYs1d( data = productData['data'], axes = axes, interpolation = interpolation )
        else :
            for xy in productData['data'] :
                if( xPrior is not None ) :
                    if( xy[0] < xPrior ) : raise ValueError( 'xy[0] = %s < xPrior = %s for MT=%d, MF=5' % ( xy[0], xPrior, MT ) )
                    if( xy[0] == xPrior ) :
                        xy[0] *= ( 1 + FUDGE_EPS )
                        if( addWarning ) : warningList.append( 'weights have same energies, second one being incremented for MT=%d, MF=5' % MT )
                        addWarning = False
                xPrior = xy[0]
        interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS1d( productData['interpolationInfo'][0][1] )
        axes = multiplicityAxes.copy( )
        axes[0].label = ''
        weights.append( XYs1dModule.XYs1d( data = productData['data'], axes = axes, interpolation = interpolation ) )   # weights is only used for delayed nu_bar data

        info.logs.write( ' : MF=5, LF=%s' % LF )
        U = physicalQuantityModule.U( PQUModule.PQU_float.surmiseSignificantDigits( productData['C1'] ), 'eV' ) # Upper energy limit.
        if( LF == 1 ) :
            dataLine, EEpETable = endfFileToGNDSMiscModule.getTAB2_TAB1s( dataLine, MF5Data, logFile = info.logs, axes = energyAxes )
            subform = toPointwiseOrPiecewiseEnergy( MT, EEpETable )
        elif( LF == 5 ) :
            dataLine, TAB1, thetas = endfFileToGNDSMiscModule.toEnergyFunctionalData( info, dataLine, MF5Data, 5, 'theta', 'eV' )
            dataLine, TAB1, gs = endfFileToGNDSMiscModule.toEnergyFunctionalData( info, dataLine, MF5Data, 5, 'g', '',
                    xLabel = 'energy_out / theta( energy_in )', xUnit = '' )
            subform = energyModule.GeneralEvaporation( U, thetas, gs )
        elif( LF == 7 ) :
            dataLine, TAB1, thetas = endfFileToGNDSMiscModule.toEnergyFunctionalData( info, dataLine, MF5Data, 7, 'theta', 'eV' )
            subform = energyModule.SimpleMaxwellianFission( U, thetas )
        elif( LF == 9 ) :
            dataLine, TAB1, thetas = endfFileToGNDSMiscModule.toEnergyFunctionalData( info, dataLine, MF5Data, 9, 'theta', 'eV' )
            subform = energyModule.Evaporation( U, thetas )
        elif( LF == 11 ) :
            dataLine, TAB1, a = endfFileToGNDSMiscModule.toEnergyFunctionalData( info, dataLine, MF5Data, 11, 'a', 'eV' )
            dataLine, TAB1, b = endfFileToGNDSMiscModule.toEnergyFunctionalData( info, dataLine, MF5Data, 11, 'b', '1/eV' )
            subform = energyModule.Watt( U, a, b )
        elif( LF == 12 ) :
            dataLine, TAB1, Ts = endfFileToGNDSMiscModule.toEnergyFunctionalData( info, dataLine, MF5Data, 12, 'T_M', 'eV' )
            EFL = physicalQuantityModule.EFL( PQUModule.PQU_float.surmiseSignificantDigits( TAB1['C1'] ), 'eV' )
            EFH = physicalQuantityModule.EFH( PQUModule.PQU_float.surmiseSignificantDigits( TAB1['C2'] ), 'eV' )
            subform = energyModule.MadlandNix( EFL, EFH, Ts )
        else :
            raise ValueError( "Unsupported LF = %d" % LF )
        if( not( delayNeutrons ) and ( NK > 1 ) ) : subform.weight = weight

        energySubforms.append( subform )

    if( not( delayNeutrons ) ) :
        if( NK > 1 ) :
            info.logs.write( ' -- using energy.WeightedFunctionals subform --' )
            weightedSubform = energyModule.WeightedFunctionals( )
            for functional in energySubforms :
                weight = energyModule.Weighted( functional.weight, functional )
                weightedSubform.addWeight( weight )
                del functional.weight
            energySubforms = weightedSubform
        else :
            energySubforms = energySubforms[0]
    return( energySubforms, weights )


def readMF6(MT, info, MF6Data, productList, warningList, undefinedLevelInfo, isTwoBody, crossSection, LR, compoundZA=None):

    twoBodyIndex = 0
    ZA, AWR, JP, LCT, NK, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF6Data[0], logFile = info.logs )
    ZA = int( ZA )
    doInfo = not (ZA == 0 and AWR != 0)
    printAWR_mode(info, MT, 6, 0, ZA, AWR, doInfo)
    info.addMassAWR( ZA, AWR )
    JPP=int(JP)//10
    JPN=int(JP)%10
    if JP > 0:
        warningList.append(warningModule.NotImplemented("Multiplicity dependent fission data, P(nu)"))

    LCT = int( LCT )
    if( LCT == 0 ) :        # Happens for electro-atomic data.
        LCT = 1
        if( MT in (525, 526) ) : LCT = 2       # Not sure this is right. Maybe it should be 1.
    LCTLight, LCTWeight = LCT, LCT
    if( LCT in (3,4) ) : LCTLight, LCTWeight = 2, 1
    NK = int( NK )                  # number of outgoing particle data sets

    dataLine, discreteGammas, discretePrimaryGammas = 1, {}, []
    info.logs.write( ' : MF=6' )
    missingMetaStable = []
    for outGoing in range( NK ) :
        isLegendreConvertedToEnergy = False
        AWP_lineNumber = dataLine
        dataLine, productData, multiplicityRegions = endfFileToGNDSMiscModule.getTAB1Regions( dataLine, MF6Data, logFile = info.logs,
                axes = multiplicityAxes )
            # ZAP is ZA for outgoing particle; AWP is its atomic mass, LIP: 0 for residual in ground state, 1 for first excited state, etc
        ZAP, AWP, LIP, LAW, NP = int( productData['C1'] ), productData['C2'], productData['L1'], productData['L2'], productData['NR']
        if( ZAP < 2005 ) : LIP = 0 # LIP has multiple meanings. For light products, signifies different products with the same ZAP.
        if( ZAP == 11 ) :          # 11 is the ZAP for an electron, however, ENDF/B-VII.1 mislabels the photo as 11 also.
            ZAP = 0
            if( ( LAW == 8 ) or ( MT in [ 525, 526, 528 ] ) or ( MT >= 534 ) ) : ZAP = 9
            if( ZAP == 0 ) : warningList.append( 'photon most likely labelled as an electron (ZAP = 11), converting to ZAP = 0' )
        printAWR_mode(info, MT, 6, AWP_lineNumber, ZAP, AWP, LIS=LIP)
        info.addMassAWR( ZAP, AWP, asTarget=False )
        LCTp = LCTLight
        if( ZAP % 1000 > 4 ) : LCTp = LCTWeight
        if LCT==4:
            if outGoing == 0:
                LCTp = 2
                ZA_recoilPartner = info.projectileZA + info.targetZA - ZAP
            elif outGoing == 1:     # check whether this is the recoil partner or a break-up product
                if ZAP == ZA_recoilPartner:
                    LCTp = 2
                else:
                    LCTp = 1
            else:
                LCTp = 1 # lab frame for all remaining break-up products.

        LAW = int( LAW )
        frame = frames[LCTp]

        info.logs.write( ' : ZAP=%s, LAW=%s' % ( ZAP, LAW ) )
        form = None
        averageProductEnergy = None
        if( LAW == 0 ) :
            form = unspecifiedModule.Form( info.style, frame )
        elif( LAW == 1 ) :              # Continuum Energy-Angle distributions
            dummy, dummy, LANG, LEP, NR, NE = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )
            LANG = int( LANG )          # identifies the type of data
            info.logs.write( ', LANG=%s' % LANG )
            LEP = int( LEP )            # interpolation type for outgoing energy
            interpolationLEP = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS1d( LEP )
            if( LANG == 1 ) :               # Legendre expansion
                NR = int( NR )              # number of interpolation regions for incident energy
                if( NR != 1 ) : raise Exception( "Currently only one interpolation region is supported for MF = 6, LAW = 1, LANG = 2; MT = %s" % MT )
                NE = int( NE )              # number of incident energies

                dataLine += 1
                EinInterpolationTypes = endfFileToGNDSMiscModule.nStringsToInts( NR, dataLine, MF6Data, dimension = 2 )
                interpolationQualifier, EinInterpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( EinInterpolationTypes[0][1] )
                dataLine += 1 + ( NR - 1 ) // 3    # the next data is energy-angle distributions

                massRatio = AWR / ( 1. + AWR )
                maxLegendre = 0
                EEpClsData = []
                for EinCount in range( NE ) :
                    dummy, Ein, ND, NA, NW, NEP = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF6Data[dataLine], logFile = info.logs )
                    ND = int( ND )          # number of discrete gammas (nonzero only for gammas)
                    NA = int( NA )          # number of angular parameters (i.e., lMax).
                    NW = int( NW )          # number of data values for this incident energy
                    NEP = int( NEP )        # number of outgoing energy values
                    maxLegendre = max( maxLegendre, NA )
                    dataLine += 1

                    if ND != 0:
                        if NA > 0: raise Exception("Logic here currently only supports isotropic scattering.")
                        discreteGammasAtE, EoutData = endfFileToGNDSMiscModule.readDiscreteAndLegendre(
                            ND, NEP - ND, dataLine, MF6Data, dimension=2 + NA, logFile=info.logs)

                        discreteNorm = sum([discrete[-1][0] for discrete in discreteGammasAtE])
                        discretePrimaryGammasAtE = []
                        for Eg, duplicateCounter, P in discreteGammasAtE:               # Divide gammas into primary and discrete.
                            if Eg < 0:                                                  # Primary gamma.
                                Epg = -Eg
                                if not info.convertJENDL_stylePrimarygammas:
                                    if( Ein > 1e-2 ) : Epg -= massRatio * Ein           # Only do projectile correction for Ein > 1e-2 as typically Eg ~ 1e6.
                                discretePrimaryGammasAtE.append( [ Epg, [ Ein, P ] ] )
                            else :                                                      # Discrete gamma.
                                if( (Eg, duplicateCounter) not in discreteGammas ) : discreteGammas[(Eg,duplicateCounter)] = []
                                discreteGammas[(Eg,duplicateCounter)].append( [ Ein, P ] )
                        if( len( discretePrimaryGammasAtE ) > 0 ) :
                            if( len( discretePrimaryGammas ) == 0 ) :
                                discretePrimaryGammas = discretePrimaryGammasAtE
                            else :      # Now we need to match level energies (aka binding energies) from prior with Ein's with current Ein.
                                        # This is needed since for different Ein's, the calculation of Epg will vary slightly (hopefully less than 1e-4).
                                if( len( discretePrimaryGammas ) != len( discretePrimaryGammasAtE ) ) :
                                    raise Exception( 'number of primary gammas at Ein = %s is different than first incident energy' % Ein )
                                for index, Eg in enumerate( discretePrimaryGammas ) :
                                    if( abs( Eg[0] - discretePrimaryGammasAtE[index][0] ) > 1e-4 * Eg[0] ) :
                                        raise Exception( 'For Ein = %s, primary energy of %s is not close to primary energy %s'
                                                         % ( Ein, Eg[0], discretePrimaryGammasAtE[index][0] ) )
                                    Eg.append( discretePrimaryGammasAtE[index][1] )
                    else :
                        EoutData = endfFileToGNDSMiscModule.nFunkyFloatStringsToFloats(
                            NEP, dataLine, MF6Data, dimension=2 + NA, logFile=info.logs)
                        discreteNorm = 0

                    if len(EoutData) > 0:
                        # Test for some common problems.
                        # 1: all outgoing energies == 0 for continuum gammas.
                        # This usually only affects one incident energy, so just replace that energy with empty outgoing distribution:
                        energy_out_list = [ e[0] for e in EoutData ]
                        if( sum( energy_out_list ) == 0 ) :
                            warningList.append("At Ein=%s eV, continuum gamma outgoing energies are all 0.0 (MT=%d, ZAP=%d)!"
                                    % ( Ein, MT, ZAP ) )
                            EoutData = [ [ 0.0, 0.0 ], [ 1.0, 0.0 ] ]
                        # 2: trouble with duplicate outgoing energies:
                        elif( ( max( [ energy_out_list.count( a ) for a in energy_out_list ] ) > 2 ) or ( energy_out_list.count( 0.0 ) > 1 ) ) :
                            warningList.append( "Too many duplicate outgoing energies for Ein=%s eV (MT=%d, ZAP=%d)!"
                                    % ( Ein, MT, ZAP ) )
                            tmp = []
                            i1 = 0
                            if( energy_out_list.count( 0.0 ) > 1 ) :
                                finalZeroIndex = energy_out_list.count( 0.0 ) - 1
                                energy_out_list = energy_out_list[finalZeroIndex:]
                                EoutData = EoutData[finalZeroIndex:]
                            while( i1 < len( energy_out_list ) ) :
                                eout = energy_out_list[i1]
                                eoutCount = energy_out_list.count( eout )
                                if( eoutCount > 2 ) :
                                    tmp.extend( [ EoutData[i1], EoutData[i1+eoutCount-1] ] )
                                    i1 += eoutCount
                                else :
                                    tmp.append( EoutData[i1] )
                                    i1 += 1
                            EoutData = tmp                          # End of TENDL-specific tests.
                        EpClsDatas, EpPrior, ClsPrior = [], -1, []
                        if( LEP == 2 ) :                            # Remove some data that are not needed and complicates things.
                            if( len( EoutData ) > 2 ) :
                                if( EoutData[-1][0] == EoutData[-2][0] ) :
                                    if( EoutData[-1][1] == EoutData[-2][1] ) :
                                        EoutData.pop()
                                    elif (EoutData[-1][1] == 0):   # Discontinuity used to make final point == 0. Bump previous energy back slightly
                                         EoutData[-2][0] = nudgeValue( EoutData[-2][0], -1 )
                        regions = []
                        skippedLast = False
                        for i1, EpCs in enumerate( EoutData ) :
                            e_out, Cls = EpCs[0], EpCs[1:]
                            if( e_out == EpPrior ) :
                                if( Cls == ClsPrior ) :
                                    if( skippedLast ) : raise ValueError( 'hell - need to fix' )
                                    warningList.append( 'skipping duplicate e_out = %s, i1 = %s %s %s' % ( e_out, i1, EinCount, Ein ) )
                                    skippedLast = True
                                    continue
                                else:   # create new region:
                                    regions.append( EpClsDatas )
                                    EpClsDatas = []
                            skippedLast = False
                            EpClsData = [ e_out, Cls ]
                            EpClsDatas.append( EpClsData )
                            EpPrior, ClsPrior = e_out, Cls
                        if regions:
                            regions.append( EpClsDatas )
                            EEpClsData.append( [ Ein, regions, True ] )
                        else:
                            EEpClsData.append( [ Ein, EpClsDatas, False ] )

                        # integrate continuum probability:
                        Ein, regionsData, isRegions = EEpClsData[-1]
                        if not isRegions:
                            regionsData = [regionsData]
                        continuumNorm = 0
                        for region in regionsData:
                            Eouts = [r[0] for r in region]
                            probs = [r[1][0] for r in region]
                            continuumRegion = XYs1dModule.XYs1d((Eouts, probs), dataForm="XsAndYs", interpolation=interpolationLEP)
                            continuumNorm += continuumRegion.integrate()
                    else:
                        continuumNorm = 0

                    normalization = discreteNorm + continuumNorm
                    if normalization < 0.9999 or normalization > 1.0001:
                        warningList.append(
                            "MT=%d unnormalized gammas at Ein=%s eV: sum = %.6f (discrete = %.6f, continuum = %.6f)"
                            % (MT, Ein, normalization, discreteNorm, continuumNorm))
                        # FIXME also add to info.doRaise list?

                    dataLine += 1 + ( NW -  1 ) // 6    # Offset for the next incident energy

                if( maxLegendre == 0 ) :    # Only have l=0 for each outgoing energy. Convert this to
                                            # uncorrelated with P(E_out|E_in) and isotropic angular distribution.
                    if( len( EEpClsData ) != 0 ) :  # length == 0 happens when there are only discrete gammas and no continuum gamma data.
                        isLegendreConvertedToEnergy = True
                        angularSubform = angularModule.Isotropic2d( )

                        EPrior = -1
                        energySubform = energyModule.Regions2d( axes = energyAxes )
                        energySubformRegion = energyModule.XYs2d( axes = energyAxes, interpolation = EinInterpolation,
                                interpolationQualifier = interpolationQualifier )
                        for index, ( Ein, EpCls, multiRegion ) in enumerate( EEpClsData ) :
                            if multiRegion:
                                xData = energyModule.Regions1d( outerDomainValue = Ein, axes = energyAxes )
                                for region in EpCls:
                                    EpP = [ [ Ep, P[0] ] for Ep, P in region ]
                                    if( len( EpP ) < 2 ) : continue
                                    xData.append( energyModule.XYs1d( data = EpP, axes = energyAxes, interpolation = interpolationLEP ) )
                            else:
                                EpP = [ [ Ep, P[0] ] for Ep, P in EpCls ]
                                xData = energyModule.XYs1d( data = EpP, interpolation = interpolationLEP, outerDomainValue = Ein, axes = energyAxes )
                            if( Ein == EPrior ) :
                                energySubform.append( energySubformRegion )
                                energySubformRegion = energyModule.XYs2d( axes = energyAxes, interpolation = EinInterpolation,
                                interpolationQualifier = interpolationQualifier )
                            energySubformRegion.append( xData )
                            EPrior = Ein
                        if( len( energySubform ) > 0 ) :
                            if( len( energySubformRegion ) > 1 ) : energySubform.append( energySubformRegion )
                        else :
                            energySubform = energySubformRegion

                        form = uncorrelated( info.style, frame, angularSubform, energySubform )
                else:
                    EPrior = -1
                    energyAngularSubform = energyAngularModule.Regions3d(axes = energyAngularAxes)
                    energyAngularSubformRegion = energyAngularModule.XYs3d(axes = energyAngularAxes, interpolation = EinInterpolation,
                            interpolationQualifier = interpolationQualifier)
                    for i1, (e_in, EpCls, multiRegion) in enumerate(EEpClsData):
                        if multiRegion:
                            raise NotImplemented
                        multiD_2d = energyAngularModule.XYs2d(outerDomainValue = e_in, interpolation = interpolationLEP, axes = energyAngularAxes)
                        for i2, (e_out, Cls) in enumerate(EpCls):
                            multiD_2d.append(energyAngularModule.Legendre(coefficients = Cls, outerDomainValue = e_out, axes = energyAngularAxes))
                        if e_in == EPrior:
                            if len(energyAngularSubformRegion) > 1:
                                energyAngularSubform.append(energyAngularSubformRegion)
                            energyAngularSubformRegion = energyAngularModule.XYs3d(axes=energyAngularAxes, interpolation=EinInterpolation,
                            interpolationQualifier = interpolationQualifier)
                        energyAngularSubformRegion.append(multiD_2d)
                        EPrior = e_in
                    if len(energyAngularSubform) > 1:
                        energyAngularSubform.append(energyAngularSubformRegion)
                    else:
                        energyAngularSubform = energyAngularSubformRegion
                    form = energyAngularModule.Form( info.style, frame, energyAngularSubform)

            elif LANG == 2:                 # Kalbach-Mann data
                if LCTp != 2:
                    raise Exception('LCT = %s != 2 as required for Kalbach-Mann data for MF=6, MT=%s' % (LCTp, MT))
                dataLine, KalbachMannData = endfFileToGNDSMiscModule.getTAB2_Lists(dataLine, MF6Data, logFile=info.logs)
                if KalbachMannData['NR'] != 1:
                    raise Exception("Currently only one interpolation flag is supported for MF = 6, LAW = 1, LANG = 2; MT = %s" % MT)
                interpolationQualifier, interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd(KalbachMannData['interpolationInfo'][0][1])
                NA = int(KalbachMannData['Lists'][0]['L2'])
                dataPerPoint = NA + 2

                fData = [energyModule.XYs2d(axes=KalbachMann_f_Axes, interpolation=interpolation,
                                           interpolationQualifier=interpolationQualifier)]

                ra_interpolationQualifier = {
                    xDataEnumsModule.InterpolationQualifier.none: xDataEnumsModule.InterpolationQualifier.none,
                    xDataEnumsModule.InterpolationQualifier.unitBase: xDataEnumsModule.InterpolationQualifier.unitBaseUnscaled,
                    xDataEnumsModule.InterpolationQualifier.correspondingPoints: xDataEnumsModule.InterpolationQualifier.correspondingPoints
                }[interpolationQualifier]

                rData = [KalbachMannModule.XYs2d(axes=KalbachMann_r_Axes, interpolation=interpolation,
                                                 interpolationQualifier=ra_interpolationQualifier)]

                aData = None
                if NA == 2:
                    aData = [KalbachMannModule.XYs2d(axes=KalbachMann_a_Axes, interpolation=interpolation,
                                                     interpolationQualifier=ra_interpolationQualifier)]

                priorE, priorEp_f_r_ = -1, [-1, -1, -1]
                for i1, data in enumerate(KalbachMannData['Lists']):
                    value = data['C2']
                    Ep_ = data['data'][::dataPerPoint]
                    f_ = data['data'][1::dataPerPoint]
                    r_ = data['data'][2::dataPerPoint]
                    priorEp_f_r__ = [Ep_, f_, r_]
                    if value == priorE:
                        if i1 == 1:
                            fData[-1].pop(0)
                            rData[-1].pop(0)
                            if aData is not None:
                                aData[-1].pop(0)
                        else:
                            if priorEp_f_r_ == priorEp_f_r__:
                                continue        # For TENDL files with duplicate data.

                            # Discontinuous data, need to split into regions
                            fData.append(KalbachMannModule.XYs2d(axes=KalbachMann_f_Axes, interpolation=interpolation,
                                                                 interpolationQualifier=interpolationQualifier))
                            rData.append(KalbachMannModule.XYs2d(axes=KalbachMann_r_Axes, interpolation=interpolation,
                                                                 interpolationQualifier=ra_interpolationQualifier))
                            if aData is not None:
                                aData.append(
                                    KalbachMannModule.XYs2d(axes=KalbachMann_a_Axes, interpolation=interpolation,
                                                            interpolationQualifier=ra_interpolationQualifier))

                    priorE, priorEp_f_r_ = value, priorEp_f_r__
                    if Ep_[-2] == Ep_[-1]:
                        Ep_[-1] *= 1.00000001
                    fData[-1].append(energyModule.XYs1d(data=(Ep_, f_), dataForm='xsandys', axes=KalbachMann_f_Axes,
                                                 outerDomainValue=value, interpolation=interpolationLEP))
                    rData[-1].append(XYs1dModule.XYs1d(data=(Ep_, r_), dataForm= 'xsandys', axes=KalbachMann_r_Axes,
                                                 outerDomainValue=value, interpolation=interpolationLEP))
                    if aData is not None:
                        a_ = data['data'][3::dataPerPoint]
                        aData[-1].append(XYs1dModule.XYs1d(data=(Ep_, a_), dataForm='xsandys', axes=KalbachMann_a_Axes,
                                                     outerDomainValue=value, interpolation=interpolationLEP))

                def toRegionsIfNeeded(section):
                    if section is None: return
                    if len(section) == 1: return section[0]  # XYs2d
                    regions_ = KalbachMannModule.Regions2d()
                    for subsection in section:
                        regions_.append(subsection)
                    regions_.axes = regions_[0].axes
                    return regions_

                fSubform = KalbachMannModule.FSubform(toRegionsIfNeeded(fData))
                rSubform = KalbachMannModule.RSubform(toRegionsIfNeeded(rData))
                aSubform = KalbachMannModule.ASubform(toRegionsIfNeeded(aData))
                form = KalbachMannModule.Form(info.style, frame, fSubform, rSubform, aSubform)

            elif( LANG in [ 11, 12, 13, 14, 15 ] ) :    # P(E',mu|E)
                NR = int( NR )                      # number of interpolation regions for incident energy
                if( NR != 1 ) : raise Exception( "Currently only one interpolation region is supported for MF = 6, LAW = 1, LANG = %s; MT = %s" % (LANG, MT) )
                NE = int( NE )                      # number of incident energies

                dataLine += 1
                EinInterpolationTypes = endfFileToGNDSMiscModule.nStringsToInts( NR, dataLine, MF6Data, dimension = 2 )
                interpolationQualifier, EinInterpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( EinInterpolationTypes[0][1] )
                dataLine += 1 + ( NR - 1 ) // 3      # the next data is energy-angle distributions

                EPrior = -1
                energyAngularSubformRegion = energyAngularModule.XYs3d(axes=energyAngularAxes, interpolation=EinInterpolation,
                            interpolationQualifier=interpolationQualifier)

                muInterpolationQualifier, muInterpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( LANG - 10 )
                for EinCount in range( NE ) :
                    dummy, Ein, ND, NA, NW, NEP = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF6Data[dataLine], logFile = info.logs )
                    ND = int( ND )          # number of discrete gammas (nonzero only for gammas)
                    NA = int( NA )          # number of angular parameters (i.e., the number of mu values).
                    NW = int( NW )          # number of data values for this incident energy
                    NEP = int( NEP )        # number of outgoing energy values
                    dataLine += 1

                    if( ND != 0 ) : raise Exception( 'Discrete gammas not support for MF=6, LAW=%d, LANG=%d; MT=%d' %
                            ( LAW, LANG, MT ) )

                    if( Ein == EPrior ) : raise Exception( 'regions not implemented for MF=6, LAW=%d, LANG=%d; MT=%d; Ein = %s' %
                            ( LAW, LANG, MT, Ein ) )
                    EPrior = Ein

                    EpF0_fOfMu = endfFileToGNDSMiscModule.nFunkyFloatStringsToFloats( NW // ( NA + 2 ), dataLine, MF6Data,
                            logFile = info.logs, dimension = NA + 2 )
                    dataLine += 1 + ( NW -  1 ) // 6    # Offset for the next incident energy

                    fOfMu_givenEp = energyAngularModule.XYs2d(outerDomainValue=Ein, interpolation=interpolationLEP, axes = energyAngularAxes)
                    for data in EpF0_fOfMu :
                        Ep = data.pop( 0 )
                        f0 = data.pop( 0 )
                        fOfMu = f0 * XYs1dModule.XYs1d( data, dataForm = 'list' )
                        fOfMu = energyAngularModule.XYs1d(fOfMu, outerDomainValue=Ep, interpolation=muInterpolation, axes = energyAngularAxes)
                        fOfMu_givenEp.append( fOfMu )

                    energyAngularSubformRegion.append( fOfMu_givenEp )

                form = energyAngularModule.Form( info.style, frame, energyAngularSubformRegion )

            else :
                raise Exception( "Unsupported LANG = %d for continuum energy-angle distribution MF = 6: ZAP = %d, LAW = %d: MT = %d" % \
                        ( LANG, ZAP, LAW, MT ) )
        elif( LAW == 2 ) :
            if( LCT not in (2,4) ): raise Exception( "Discrete two-body must be in the center-of-mass frame: LCT = %d MT = %d." % ( LCT, MT ) )
            isTwoBody = True
            twoBodyIndex = len(productList)
            dataLine, angularData = endfFileToGNDSMiscModule.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
            LANG = int( angularData['Lists'][0]['L1'] )
            info.logs.write( ', LANG=%s' % LANG )
            if( angularData['NR'] != 1 ) :
                raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 2, LANG = %s; MT = %s" % ( LANG, MT ) )
            interpolationQualifier, interpolationE_in = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( angularData['interpolationInfo'][0][1] )

            if( LANG == 0 ) :
                angularSubform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 6, 'LAW = 2, LANG = 0' )
                form = angularModule.TwoBody( info.style, frame, angularSubform )
            elif( LANG in [ 12, 14 ] ) :
                angularSubform = convertAngularToPointwiseOrPiecewiseFromTAB2_List( MT, LANG, angularData, warningList )
                form = angularModule.TwoBody( info.style, frame, angularSubform )
            else :
                raise Exception( "LANG = %d for LAW = %d not supported: MT = %d" % ( LANG, LAW, MT ) )
            if( ( ZAP == 0 ) and ( AWP != 0 ) ) :
                form = angularModule.TwoBody( info.style, frame, angularSubform )
        elif( LAW == 3 ) :
            subform = angularModule.Isotropic2d( )
            form = angularModule.TwoBody( info.style, frame, subform )
        elif( LAW == 4 ) :
            subform = angularModule.Recoil( product.distribution[info.style], relative=True )
            form = angularModule.TwoBody( info.style, frame, subform )
        elif( LAW == 5 ) :  # charged-particle elastic scattering
            assert LCT == 2, "Charged-particle elastic must be in the center-of-mass frame: LCT = %d MT = %d." % ( LCT, MT )
            dataLine, angularData = endfFileToGNDSMiscModule.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
            SPI = angularData['C1']
            LIDP = angularData['L1']    # identical particles?
            info.particleSpins[info.projectile] = ( commonResonanceModule.Spin( SPI ), 0 )                         # no parity information
            LTP = int( angularData['Lists'][0]['L1'] )
            info.logs.write( ', LTP=%s' % LTP )
            interpolationQualifier, interpolationE_in = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( angularData['interpolationInfo'][0][1] )

            nuclearPlusInterference = None
            nuclearAmplitudeExpansion = None
            # LTP flag changes interpretation of the data:
            if LTP == 1:
                (nuclear, real, imaginary) = convertNuclearPlusInterferenceDataToPiecewise(MT, angularData, warningList, 6, 'LAW = 5, LTP = %d'%LTP, LIDP)
                nuclearAmplitudeExpansion = nuclearAmplitudeExpansionModule.NuclearAmplitudeExpansion(
                        nuclearTerm=nuclear, realInterference=real, imaginaryInterference=imaginary)
            elif LTP == 2:
                raise NotImplementedError("MF=6 LAW=5 LTP=2 not yet implemented (MT%d)!" % MT)
            elif LTP in (12, 14, 15):
                subform = convertAngularToPointwiseOrPiecewiseFromTAB2_List(MT, LTP, angularData, warningList)
                assert len(set([tmp.domainMax for tmp in subform])) == 1, "mu cutoff should not depend on energy!"

                muCutoff = subform[0].domainMax
                crossSection = crossSectionModule.XYs1d(data = 2*math.pi*crossSection, axes = crossSection.axes)
                nuclearPlusInterference = nuclearPlusInterferenceModule.NuclearPlusInterference(
                    muCutoff = muCutoff,
                    crossSection = nuclearPlusInterferenceModule.CrossSection(crossSection),
                    distribution = nuclearPlusInterferenceModule.Distribution(subform)
                )
            else:
                raise Exception("unknown LTP encountered for MF=6, LAW=5, MT=%s" % MT)

            dSigma_form = CoulombPlusNuclearElasticModule.Form( info.projectile, info.style, nuclearPlusInterference = nuclearPlusInterference,
                    nuclearAmplitudeExpansion = nuclearAmplitudeExpansion, identicalParticles = ( LIDP == 1 ) )
            info.dSigma_form = dSigma_form
            # also make a link from 'normal' distribution to differential part:
            form = referenceModule.CoulombPlusNuclearElastic( link=dSigma_form, label=info.style, relative=True )
        elif( LAW == 6 ) :
            APSX, dummy, dummy, dummy, dummy, NPSX = endfFileToGNDSMiscModule.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )
            dataLine += 1
            APSX *= info.massTracker.neutronMass
            angularSubform = angularModule.Isotropic2d( )
            mass = physicalQuantityModule.Mass( PQUModule.PQU_float.surmiseSignificantDigits( APSX ), 'amu' )
            energySubform = energyModule.NBodyPhaseSpace( int( NPSX ), mass )
            # Some data has the wrong frame, should always be com hence frames[2].
            form = uncorrelated( info.style, frames[2], angularSubform, energySubform )

        elif( LAW == 7 ) :
            dataLine, NEData = endfFileToGNDSMiscModule.getTAB2Header( dataLine, MF6Data, logFile = info.logs )
            NR = int( NEData['NR'] )                    # number of interpolation regions for this incident energy
            if( NR != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 7; MT = %s" % MT )
            interpolationQualifier, interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( NEData['interpolationInfo'][0][1] )
            angularEnergySubform = angularEnergyModule.XYs3d( axes = angularEnergyAxes, interpolation = interpolation,
                    interpolationQualifier = interpolationQualifier )
            for i1 in range( int( NEData['NZ'] ) ) :   # Loop over incident energies
                dataLine, muEpPTable = endfFileToGNDSMiscModule.getTAB2_TAB1s( dataLine, MF6Data, logFile = info.logs, axes = angularEnergyAxes )
                muEpP = muEpPTable['TAB2s']
                if( len( muEpP ) != 1 ) : raise ValueError( 'hell - need to fix' )
                interpolation, muEpP = muEpP[0]
                interpolationQualifier, interpolation = endfFileToGNDSMiscModule.ENDFInterpolationToGNDS2plusd( muEpPTable['interpolationInfo'][0][1] )
                muEpP_multiD = angularEnergyModule.XYs2d( outerDomainValue = muEpPTable['C2'],
                        interpolation = interpolation, interpolationQualifier = interpolationQualifier, axes = angularEnergyAxes )
                for i2, EpP in enumerate( muEpP ) :
                    if( len( EpP ) > 1 ) : raise ValueError( 'hell - need to fix' )
                    muEpP_multiD.append( angularEnergyModule.XYs1d.returnAsClass( EpP[0], EpP[0] ) )
                angularEnergySubform.append( muEpP_multiD )
            form = angularEnergyModule.Form( info.style, frame, angularEnergySubform )
        elif( LAW == 8 ) :
            dataLine, TAB1, energyLoss = endfFileToGNDSMiscModule.getTAB1Regions(dataLine, MF6Data, logFile=info.logs, axes=averageProductEnergyAxes)
            if( len( energyLoss ) != 1 ) : raise ValueError( 'hell - fix me' )
            energyLoss = energyLoss[0]
            data = [ [ energyLoss.domainMin, energyLoss.domainMin ], [ energyLoss.domainMax, energyLoss.domainMax ] ]
            energyLoss = energyLoss.__class__( data, axes = energyLoss.axes ) - energyLoss
            averageProductEnergy = averageProductEnergyModule.XYs1d( label = info.style, axes = averageProductEnergyAxes,
                    data = energyLoss, interpolation = energyLoss.interpolation )
        else :
            raise Exception( "LAW = %d not implemented: MT = %d." %  ( LAW, MT ) )

        multiplicityData = productData['data']
        rangeMin, integerMultiplicity = multiplicityData[0][1], False
        if (ZAP != 0) and (int(rangeMin) == rangeMin):
            integerMultiplicity = True
            for x, y in multiplicityData : integerMultiplicity &= (rangeMin == y)
        if integerMultiplicity:
            multiplicity = int(rangeMin)
        else:
            multiplicity = getMultiplicityPointwiseOrPieceWise(info, multiplicityRegions, warningList)
        if isinstance(multiplicity, multiplicityModule.XYs1d) and (multiplicity.rangeMin < 0):
            warningList.append("Negative multiplicity encountered for MF6, MT%d %s" %
                               (MT, toGNDSMiscModule.getTypeNameENDF(info, ZAP, undefinedLevelInfo)))

        if (ZAP == 0) and (AWP == 0):  # Gammas. Appear to always be stored using LAW = 1 and LANG = 1.

            def addGammaProduct(form, multiplicity):

                if isinstance(multiplicity, multiplicityModule.Regions1d):
                    if len(multiplicity) == 1:
                        label = multiplicity.label
                        multiplicity = multiplicity[0]
                        multiplicity.label = label
                        multiplicity.index = None
                if isinstance(multiplicity, multiplicityModule.XYs1d):
                    if len(multiplicity) < 2: return None
                product = toGNDSMiscModule.newGNDSParticle(
                    info, toGNDSMiscModule.getTypeNameENDF(info, ZAP, None), crossSection,
                    multiplicity=multiplicity)
                product.distribution.add(form)
                productList.append(product)
                return product      # May be required for LAW = 4 logic.

            def addGammaAdjustWeight(angularSubform, energySubform, totalGammaMultiplicity):
                """
                Have Legendre section with total multiplicity and only l=0 for discrete gammas.
                Convert to isotropic angular distributions, plus adjust weights.
                """

                def getPointwiseMultiplicity(angularSubform, totalGammaMultiplicity, EPrior, axes):

                    data = []
                    for energy_in in angularSubform:
                        multiplicity = getMultiplicity(totalGammaMultiplicity, EPrior, energy_in.outerDomainValue)
                        data.append([energy_in.outerDomainValue, energy_in.coefficients[0] * multiplicity])
                        EPrior = energy_in.outerDomainValue
                    multiplicity = multiplicityModule.XYs1d(data=data, axes=axes)
                    return EPrior, multiplicity

                if len(angularSubform) < 2: return None
                axes = multiplicityAxes
                if isinstance(angularSubform, angularModule.XYs2d):
                    EPrior, multiplicity = getPointwiseMultiplicity(angularSubform, totalGammaMultiplicity, angularSubform[0].outerDomainValue, axes)
                    multiplicity.label = info.style
                elif isinstance(angularSubform, angularModule.Regions2d):
                    multiplicity = multiplicityModule.Regions1d(label=info.style, axes=axes)
                    EPrior = -1
                    for region in angularSubform:
                        EPrior, multiplicityRegion = getPointwiseMultiplicity(region, totalGammaMultiplicity, EPrior, axes)
                        multiplicity.append(multiplicityRegion)
                else:
                    raise Exception('Unsupported angular subform = "%s"' % angularSubform.moniker)

                product = toGNDSMiscModule.newGNDSParticle(
                    info, toGNDSMiscModule.getTypeNameENDF(info, ZAP, None), multiplicity,
                    multiplicity=multiplicity)

                angularSubform = angularModule.Isotropic2d()
                form = uncorrelated(info.style, frame, angularSubform, energySubform)
                product.distribution.add(form)
                productList.append(product)
                return product       # May be required for LAW = 4 logic.

            def addPrimaryOrDiscreteGamma(distinctType, Eg, ELegendres, totalGammaMultiplicity, frame):

                domainMin = ELegendres[0][0]
                domainMax = ELegendres[-1][0]
                if distinctType == 'discrete':
                    energySubform = discreteOrPrimaryGamma(energyModule.DiscreteGamma, Eg, domainMin, domainMax)
                else:
                    energySubform = discreteOrPrimaryGamma(energyModule.PrimaryGamma, Eg, domainMin, domainMax)

                angularSubform = angularModule.XYs2d(axes=angularAxes)
                EPrior = -1
                maxLegendre = 0
                angularRegion = angularModule.Regions2d(axes=angularAxes)
                for i1, (energy, coefficients) in enumerate(ELegendres):
                    maxLegendre = max(maxLegendre, len(coefficients))
                    if energy == EPrior:
                        angularRegion.append(angularSubform)
                        angularSubform = angularModule.XYs2d(axes=angularAxes)
                    angularSubform.append(
                        angularModule.Legendre(coefficients=coefficients, outerDomainValue=energy, axes=angularAxes))
                    EPrior = energy
                if len(angularRegion) > 0:
                    angularRegion.append(angularSubform)
                    angularSubform = angularRegion

                if maxLegendre == 1:    # Only have l=0 for each incident energy. Distribution needs to be normalized.
                    return addGammaAdjustWeight(angularSubform, energySubform, totalGammaMultiplicity)

                form = uncorrelated(info.style, frame, angularSubform, energySubform)

                return addGammaProduct(form, totalGammaMultiplicity)

            totalGammaMultiplicity = multiplicity               # Makes the name more explicit.
            if len(discreteGammas) + len(discretePrimaryGammas) > 0:
                info.totalMF6_12_13Gammas[MT] = [6, multiplicity]

            if form is not None:                       # Continuum gamma present

                if discreteGammas or discretePrimaryGammas:

                    def fixContinuumGammaSpectrum(energyDist):

                        EpPrior = 2e-5
                        energyDist[0] = (0, 1e-16)            # Make distribution with small probability near 0.
                        energyDist[1] = (EpPrior, 0)
                        for i1 in range(2, len(energyDist)):
                            Ep, P1 = energyDist[i1]
                            if Ep < EpPrior: Ep = 2 * EpPrior
                            EpPrior = Ep
                            energyDist[i1] = (Ep, 0)

                    if isinstance(totalGammaMultiplicity, XYs1dModule.XYs1d):
                        newMultiplicity = multiplicityModule.XYs1d(
                            data=[], axes=totalGammaMultiplicity.axes,
                            interpolation=totalGammaMultiplicity.interpolation, label=info.style)
                    elif isinstance(totalGammaMultiplicity, regionsModule.Regions):
                        newMultiplicity = multiplicityModule.Regions1d(axes=totalGammaMultiplicity.axes, label=info.style)
                    else:
                        raise Exception('Unsupported multiplicity form "%s"' % totalGammaMultiplicity.moniker)
                    if isinstance(energySubform, energyModule.XYs2d): energySubform = [energySubform]
                    EPrior = energySubform[0][0].outerDomainValue
                    for i1, region in enumerate(energySubform):
                        for i2, energyDist in enumerate(region):
                            Ein = energyDist.outerDomainValue
                            integral = energyDist.integrate()
                            if integral == 0:
                                if not info.continuumSpectraFix:
                                    raise Exception(
                                        "Zero norm continuum gamma spectrum energies = %s (MT=%d). Try option 'continuumSpectraFix'"
                                        % (energyDist.outerDomainValue, MT))
                                fixContinuumGammaSpectrum(energyDist)
                                integral = energyDist.integrate()
                            if isinstance(newMultiplicity, XYs1dModule.XYs1d):
                                newMultiplicity.setValue(Ein, getMultiplicity(totalGammaMultiplicity, EPrior, Ein) * integral)
                            elif isinstance(newMultiplicity, regionsModule.Regions):
                                if Ein == EPrior:
                                    if i1 + i2 > 0: newMultiplicity.append(newRegionMultiplicity)
                                    newRegionMultiplicity = multiplicityModule.XYs1d(
                                        axes=totalGammaMultiplicity.axes, data=[])
                                newRegionMultiplicity.setValue(
                                    Ein,getMultiplicity(totalGammaMultiplicity, EPrior, Ein) * integral)
                            region[i2] = region[i2].normalize()
                            EPrior = Ein
                    if isinstance(energySubform, list): energySubform = energySubform[0]
                    form = uncorrelated(info.style, form.productFrame, form.angularSubform.data, energySubform)
                    if isinstance(newMultiplicity, multiplicityModule.Regions1d):
                        newMultiplicity.append(newRegionMultiplicity)
                    product = addGammaProduct(form, newMultiplicity)
                else:
                    product = addGammaProduct(form, totalGammaMultiplicity)
            for Eg_key in sorted(discreteGammas):
                product = addPrimaryOrDiscreteGamma(
                    'discrete', Eg_key[0], discreteGammas[Eg_key], totalGammaMultiplicity, frame)
            for EgEinPs in discretePrimaryGammas:
                product = addPrimaryOrDiscreteGamma(
                    'primary', EgEinPs[0], EgEinPs[1:], totalGammaMultiplicity, frame)
        else:                     # Non gamma particle.
            if (info.targetLevel > 0) and (MT in range(51, 91)) and (ZAP == info.targetZA):
                # for isomeric targets, need to adjust the product level: MT51 goes to ground state, etc.
                if MT - 50 <= info.targetLevel: undefinedLevelInfo['levelIndex'] -= 1
            thisParticle = toGNDSMiscModule.getTypeNameENDF(info, ZAP, undefinedLevelInfo)
            if LIP > 0:
                thisParticleBaseName = thisParticle.id
                metaStableName = PoPsAliasModule.MetaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex(
                    thisParticleBaseName, int(LIP))
                if metaStableName in info.PoPs:
                    levelName = info.PoPs[metaStableName].pid
                    thisParticle = info.reactionSuite.PoPs[levelName]
                else:
                    try:
                        levelIndex, level = metaStableDataModule.metaStableData[metaStableName]
                    except:
                        levelIndex, level = 1, metaStableDataModule.unknownMetaStableEnergy
                        missingMetaStable.append(metaStableName)
                    try:
                        thisParticle = toGNDSMiscModule.getTypeNameGamma(info, ZAP, level=level, levelIndex=levelIndex)
                        info.PoPs.add(PoPsAliasModule.MetaStable(metaStableName, thisParticle.id, int(LIP)))
                    except:
                        raise KeyError('Meta stable data not available for %s' % metaStableName)
            product = toGNDSMiscModule.newGNDSParticle(info, thisParticle, crossSection, multiplicity = multiplicity)
            if form is not None: product.distribution.add(form)
            if averageProductEnergy is not None: product.averageProductEnergy.add(averageProductEnergy)

            productList.append( product )

        if (len(discreteGammas) > 0) and (productData['interpolationInfo'][0][1] != 2):
            info.logs.write('interpolation %s is not linear for gamma multiplicity'
                            % productData['interpolationInfo'][0][1])

    if len(missingMetaStable) > 0:
        print()
        for k in missingMetaStable:
            print('    "%s": (1, unknownMetaStableEnergy),' % k)
        raise KeyError('Meta stable data not available for %s' % missingMetaStable)

    if twoBodyIndex != 0:
        product = productList.pop(twoBodyIndex)
        productList.insert(0, product)

    for product in productList:
        if product.pid == IDsPoPsModule.photon:
            info.ENDFconversionFlags.add(product, 'MF6')
        elif info.reactionSuite.projectile in (IDsPoPsModule.neutron, IDsPoPsModule.photon):
            if isTwoBody or (product.pid == IDsPoPsModule.neutron):
                info.ENDFconversionFlags.add(product, 'MF6')

    if LR == 1 and MT != 5:
        residualZA = compoundZA
        for product in productList:
            if product.pid == IDsPoPsModule.photon:
                continue
            residualZA -= int(particleZA(info, product.pid) * product.multiplicity.getConstant())
        if residualZA != 0:
            if residualZA < 0:
                raise Exception('Negative residual ZA calculated as %s.' % residualZA)
            thisParticle = toGNDSMiscModule.getTypeNameGamma(info, residualZA)
            product = toGNDSMiscModule.newGNDSParticle( info, thisParticle, crossSection)
            productList.append(product)

    return( isTwoBody )


def readMF8(info, MT, MTData, warningList):
    """Regular decay data."""

    firstLMF, radioactiveDatas = None, []
    if 8 in MTData:
        dataLine, MF8Data = 1, MTData[8]
        ZA, AWR, LIS, LISO, NS, NO = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2, 3, 4, 5 ], logFile = info.logs )
        info.addMassAWR( ZA, AWR )
        printAWR_mode(info, MT, 8, 0, ZA, AWR, LIS=LIS)
        MF9Data = readMF9or10( info, MT, MTData, 9, LIS, warningList )
        MF10Data = readMF9or10( info, MT, MTData, 10, LIS, warningList )
        metastables = {}
        for idx in range( NS ) :
            ZAP, ELFS, LMF, LFS, ND6, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats(
                MF8Data[dataLine], intIndices = [ 0, 2, 3, 4 ], logFile = info.logs )

            if LMF in (9, 10):
                if LMF == 9:
                    IZAP = MF9Data[idx][0]['L1']
                else:
                    if MF10Data is None:                        # Happens for some JENDL-5 data.
                        IZAP = ZAP
                    else:
                        if idx > len(MF10Data)+1:
                            warningList.append("For MT=%d, MF8 lists %s residuals but MF10 only contains %s" %
                                               (MT, NS, len(MF10Data)))
                            info.doRaise.append(warningList[-1])
                            return None, []
                        IZAP = MF10Data[idx][0]['L1']
                if ZAP != IZAP:
                    warningList.append("For MT=%d, inconsistent ZAP in MF8 and MF%d: %d vs %d" % (MT, LMF, ZAP, IZAP))
                    info.doRaise.append(warningList[-1])
                    return None, []

            if MT == 18:
                if MF9Data is None and MF10Data is None:    # meaningless section, appears in many JENDL-5 evaluations
                    return None, []
                if ZAP != -1:
                    warningList.append("For sub-actinide fission (MT=18 MF=8/10) expected ZAP=-1, got %d" % ZAP)
                    if not info.acceptBadMF10FissionZAP:
                        info.doRaise.append(warningList[-1])
                radioactiveDatas.append([0, 0, 0, None, MF10Data[0][-1], 0, ELFS])
                return 10, radioactiveDatas

            if ZAP not in metastables: metastables[ZAP] = 0
            if LFS == 0 and ELFS != 0:
                warningList.append("MF8 claims non-zero ELFS = %s for the ground state, MT = %s. Converting ELFS to 0"
                                   % (ELFS, MT))
                info.doRaise.append(warningList[-1])
                ELFS = 0
            if ELFS != 0: metastables[ZAP] += 1

            if firstLMF is None: firstLMF = LMF
            if LMF != firstLMF: raise Exception('LMF changing from %s to %s is not allowed' % (firstLMF, LMF))
            dataLine += 1

            crossSection, multiplicity = None, 1
            QM, QI = None, None
            if LMF == 3:
                pass
            elif LMF == 6:  # multiplicity in MF6, which hasn't been read. Set multiplicity to 'unspecified', to be overridden later.
                multiplicity = multiplicityModule.Unspecified( info.style )
            elif( LMF in [ 9, 10 ] ) :
                if( ( ( LMF == 9 ) and not( MF9Data ) ) or ( ( LMF == 10 ) and not( MF10Data ) ) ) :    # BRB ????? Why are we checking for bad data.
                    LMF1 = ( 10 if( LMF == 9 ) else 9 )
                    warningList.append( 'LMF = %d, but no MF%d found for MT%d. Trying LMF=%d instead' % ( LMF, LMF, MT, LMF1 ) )
                    print('     ========', warningList[-1])
                    if( ( ( LMF1 == 9 ) and not( MF9Data ) ) or ( ( LMF1 == 10 ) and not( MF10Data ) ) ) :
                        # neither MF9 or MF10 exist.
                        continue
                    LMF = LMF1
                if LMF == 9:
                    TAB1, multiplicity = MF9Data[idx]
                    multiplicity = getMultiplicityPointwiseOrPieceWise( info, multiplicity, warningList )
                else:
                    TAB1, crossSection = MF10Data[idx]
                QM, QI, IZAP, LFS9or10 = TAB1['C1'], TAB1['C2'], int( TAB1['L1'] ), int( TAB1['L2'] )
                ELFS9or10 = QM - QI
                if abs(ELFS - ELFS9or10) > 2e-4 * abs(ELFS) and abs(ELFS - ELFS9or10) > 1e-6 * abs(QM):
                    if not info.convertJENDL_stylePrimarygammas or abs(ELFS - ELFS9or10) > 1e-3 * abs(ELFS):
                        warningList.append(
                            f"MF8 residual level energy = {ELFS} for level {LIS} of ZA = {ZAP} not close "
                            f"to MF{LMF}'s value = {ELFS9or10} for MT = {MT}")
                        info.doRaise.append(warningList[-1])
                if LFS != LFS9or10:
                    warningList.append("For MT%d, MF8 claims level index = %d for ZA = %d but MF9/10 claim level index = %d"
                                       % (MT, LFS, ZAP, LFS9or10))
                    info.doRaise.append( warningList[-1] )

            radioactiveDatas.append( [ ZAP, ELFS, LFS, multiplicity, crossSection, LFS, QI ] )

            if metastables[ZAP] and ZAP != 0:
                residual = toGNDSMiscModule.getTypeNameGamma(info, ZAP, level=ELFS, levelIndex=LFS)
                if ELFS != 0:
                    if abs(residual.energy.float('eV') - ELFS) / ELFS > 0.005:
                        warningList.append(
                            f"MF8 residual level energy = {ELFS} for level {LFS} of ZA = {ZAP} not close "
                            f"to value in PoPs (likely computed from MF=3 Q-values) for MT = {MT}")

                        # is another energy level a better match?
                        levels = [nuclide.energy.float('eV') for nuclide in residual.ancestor]
                        matches = [idx for (idx, level) in enumerate(levels) if abs(level-ELFS) / ELFS < 0.005]
                        if matches:
                            warningList[-1] += f". Possible match: LFS = {matches[0]}"
                        info.doRaise.append(warningList[-1])

                isotopeName = residual.isotope.key
                aliasName = PoPsAliasModule.MetaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex(isotopeName, metastables[ZAP])
                if aliasName not in info.PoPs:
                    info.PoPs.add(PoPsAliasModule.MetaStable(aliasName, residual.id, metastables[ZAP]))

            if NO == 0: dataLine += ( ND6 + 5 ) // 6
    return firstLMF, radioactiveDatas


def readMF9or10(info, MT, MTData, MF, targetLIS, warningList):

    if MF not in MTData.keys(): return None
    dataLine, MFData, MF9or10 = 1, MTData[MF], []
    ZA, AWR, LIS, dummy, NS, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats(
        MFData[0], intIndices=[0, 2, 4], logFile=info.logs)
    ZA = int(ZA)
    printAWR_mode(info, MT, MF, 0, ZA, AWR, LIS=LIS)
    info.addMassAWR(ZA, AWR)
    if LIS != targetLIS:
        warningList.append("residual's LIS = %s does not match target's LIS = %s: MT=%d, MF=%d"
                           % (LIS, targetLIS, MT, MF))
        info.doRaise.append(warningList[-1])
    if MF == 10:
        axes = crossSectionAxes
        XYsclass = crossSectionModule.XYs1d
    else:
        axes = axesModule.Axes(2, labelsUnits={1: ('energy_in', 'eV'), 0: ('weight', '')})
        XYsclass = multiplicityModule.XYs1d
    for i1 in range(NS):
        dataLine, TAB1, regions = endfFileToGNDSMiscModule.getTAB1Regions(
            dataLine, MFData, axes=axes, logFile=info.logs, cls=XYsclass)
        if MF == 10:
            XSec = getCrossSectionForm(info, regions)
        else :
            XSec = regions
        MF9or10.append([TAB1, XSec])
    return MF9or10


def readMF12_13(info, MT, MTData, productList, warningList, crossSection, _dummyCrossSection, gammaBRTolerance=1e-6):

    def addMF12_13GammaToList(gList, EGk, ESk, LP, LF, regions):
        """EGk is the gamma's energy and ESk is the gamma's origination level energy."""

        if EGk in gList: raise ValueError('Gammas with the same energy (%s) are not supported: MT=%s' % (EGk, MT))
        gList.append({'EGk': EGk, 'ESk': ESk, 'LP': LP, 'LF': LF, 'yield': regions, 'angularSubform': None, 'energySubform': None})

    def checkAngularData(gamma):

        angularSubform = gamma['angularSubform']
        if angularSubform is None: angularSubform = angularModule.Isotropic2d()
        gamma['angularSubform'] = angularSubform

    def addGammaProduct(info, MF, gamma, productList, warningList, ESk):

        yields = gamma['yield']
        conversionFlags = []
        multiplicity = getMultiplicityPointwiseOrPieceWise(info, yields, warningList)
        if MF == 13:
            conversionFlags.append('MF13')
            multiplicity._temp_divideByCrossSection = True
        if gamma['ESk'] != 0:
            conversionFlags.append('ESk=%s' % ESk)

        pid = toGNDSMiscModule.getTypeNameENDF(info, 0, None)
        product = toGNDSMiscModule.newGNDSParticle(info, pid, multiplicity, multiplicity=multiplicity)

        angularSubform = gamma['angularSubform']
        energySubform = gamma['energySubform']
        form = uncorrelated(info.style, frames[1], angularSubform, energySubform)
        product.distribution.add(form)
        if conversionFlags: info.ENDFconversionFlags.add(product, ','.join(conversionFlags))
        productList.append(product)

    if 12 in MTData:
        if 13 in MTData:
            raise Exception('MF = 12 and 13 present for MT=%s, this is not supported' % MT)
        MF = 12
    elif 13 in MTData:
        MF = 13
    elif (14 in MTData) or (15 in MTData):
        warningList.append('MF 14 and/or 15 data and no MF 12 or 13 data: MT=%s MFs=%s' % (MT, MTData.keys()))
        info.doRaise.append(warningList[-1])
        return
    else:
        return

    MF12_13Data = MTData[MF]
    MF6gammas = [product for product in productList if product.pid == IDsPoPsModule.photon]
    if MF6gammas:
        warningList.append(f'MT{MT} photons are likely double-counted, in MF6 and MF{MF}')
        info.doRaise.append(warningList[-1])
    ZA, AWR, LO, LG, NK, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats(
        MF12_13Data[0], intIndices=[0, 2, 3, 4], logFile=info.logs)
    printAWR_mode(info, MT, MF, 0, ZA, AWR)
    info.addMassAWR(ZA, AWR)
    info.logs.write(' : MF=%s LO=%s : ZAP=0 ' % (MF, LO))
    dataLine, continuousGamma, discreteGammas, primaryGammas, branchingGammas, regions = 1, [], [], [], [], None
    if ((MF == 12) and (LO == 1)) or ((MF == 13) and (LO == 0)):
        if MF == 12:
            axes = multiplicityAxes
        else:
            axes = crossSectionAxes
        if NK > 1:
            dataLine, TAB1, regions = endfFileToGNDSMiscModule.getTAB1Regions(
                dataLine, MF12_13Data, axes=axes, logFile=info.logs)
            info.totalMF6_12_13Gammas[MT] = [MF, getMultiplicityPointwiseOrPieceWise(info, regions, warningList)]
        for i in range(NK):
            dataLine, TAB1, regions = endfFileToGNDSMiscModule.getTAB1Regions(
                dataLine, MF12_13Data, axes=axes, logFile=info.logs)
            EGk, ESk, LP, LF = TAB1['C1'], TAB1['C2'], int(TAB1['L1']), int(TAB1['L2'])
            if EGk == 0.:
                if LP not in [0, 1]: raise Exception('LP = %s != 0 for continuous gamma for MT = %s' % (LP, MT))
                if len(continuousGamma) == 0:
                    addMF12_13GammaToList(continuousGamma, EGk, ESk, LP, LF, regions)
                else:
                    raise Exception('continuous gamma information for MF=%s, MT = %s already exist' % (MF, MT))
            elif LP == 2:
                addMF12_13GammaToList(primaryGammas, EGk, ESk, LP, LF, regions)
            else:
                addMF12_13GammaToList(discreteGammas, EGk, ESk, LP, LF, regions)
    elif (MF == 12) and (LO == 2):
        dataLine, LO2 = endfFileToGNDSMiscModule.getList(dataLine, MF12_13Data, logFile = info.logs)
        LP = int(LO2['L1'])
        if LP == 2 and MT not in (91, 649, 699, 749, 799, 849, 999):
            warningList.append("Incorrect 'primary gamma' flag for MF12 MT%d" % MT)
        NT, LGp = LO2['N2'], LG + 1
        NK = NT
        for idx in range(NT):
            parentEnergy, finalEnergy = LO2['C1'], LO2['data'][idx*LGp]
            if parentEnergy == 0:
                raise Exception("Gamma decay from ground state in MF12 MT%d" % MT)
            if abs((parentEnergy-finalEnergy)/parentEnergy) < 0.0001:
                raise Exception("Zero-energy gamma from %f eV to %f eV in MF12 MT%d" % (parentEnergy, finalEnergy, MT))
            branchingGammas.append({'ES': LO2['C1'], 'EGk': 0, 'ESk': LO2['data'][idx*LGp], 'angularSubform': None,
                                     'LG': LG, 'branching': LO2['data'][idx*LGp+1:idx*LGp+LGp]})
            if LG==2:
                GammaP = branchingGammas[-1]['branching'][-1]
                if not (0 <= GammaP <= 1):
                    warningList.append(f"Unphysical gamma emission probability {GammaP} encountered in MF=12 MT={MT}")
                    info.doRaise.append(warningList[-1])
        gammaBRList = [g['branching'][0] for g in branchingGammas]
        sumGammaBRList = sum(gammaBRList)
        if abs(sumGammaBRList - 1.0) > gammaBRTolerance:
            warningList.append("sum of gamma BR's for MT="+str(MT)+" MF=12 is " + str(sumGammaBRList)+' != 1.0')
        info.MF12_LO2[MT] = branchingGammas
    else:
        raise Exception('LO=%s is not valid for MF=%s, MT=%s' % (LO, MF, MT))

    readMF14(info, MT, MTData, MF, NK, warningList, discreteGammas, primaryGammas, continuousGamma, branchingGammas)
    readMF15(info, MT, MTData, continuousGamma, warningList)

    for gamma in branchingGammas:
        if gamma['angularSubform'] is not None:
            info.logs.write('NON-isotropic gamma')
            break
    if 14 in MTData: info.logs.write(': MF=14 ')
    if 15 in MTData: info.logs.write(': MF=15 ')

    if len(continuousGamma):
        gamma = continuousGamma[0]
        checkAngularData(gamma)
        addGammaProduct(info, MF, gamma, productList, warningList, gamma['ESk'])

    if crossSection is None and regions is not None:
        crossSection = DummyCrossSection(regions[0].domainMin, regions[-1].domainMax, 'eV')
        _dummyCrossSection.append(crossSection)

    if crossSection is not None:
        for gamma in discreteGammas:
            checkAngularData(gamma)
            gamma['energySubform'] = discreteOrPrimaryGamma(energyModule.DiscreteGamma, gamma['EGk'],
                                                            crossSection.domainMin, crossSection.domainMax)
            addGammaProduct(info, MF, gamma, productList, warningList, gamma['ESk'])
        for gamma in primaryGammas:
            checkAngularData(gamma)
            gamma['energySubform'] = discreteOrPrimaryGamma(energyModule.PrimaryGamma, gamma['EGk'],
                                                             crossSection.domainMin, crossSection.domainMax)
            addGammaProduct(info, MF, gamma, productList, warningList, gamma['ESk'])
        for gamma in branchingGammas: checkAngularData(gamma)


def readMF14(info, MT, MTData, MF, NK, warningList, discreteGammas, primaryGammas, continuousGamma, branchingGammas):

    allGammas = discreteGammas + primaryGammas + continuousGamma + branchingGammas

    def addAngularData( EGk, ESk, subform ) :

        possibleGammas = []
        for gamma in allGammas :
            if( ( abs( EGk - gamma['EGk'] ) < 1e-4 ) and ( abs( ESk - gamma['ESk'] ) < 1e-4 ) ) : possibleGammas.append( gamma )
        if( len( possibleGammas ) == 0 ) :
            sys.stderr.write( "MF = 14 gamma data: EGk = %s, ESk = %s\n" % ( EGk, ESk ) )
            sys.stderr.write( 'MF = %s gamma list is,\n' % MF )
            for gamma in allGammas : sys.stderr.write( "    EGk = %s, ESk = %s\n" % ( gamma['EGk'], gamma['ESk'] ) )
            raise Exception( 'No matching gamma from MF = %s found for MF = 14 gamma: MT=%s: see above' % ( MF, MT ) )
        if( len( possibleGammas ) != 1 ) : raise Exception( 'Multiple possible gammas found for MF=14 data: EGk=%s, ESk=%s, MT=%s' % ( EGk, ESk, MT ) )
        gamma = possibleGammas[0]
        if( gamma['angularSubform'] is None ) :
            gamma['angularSubform'] = subform
        else :
            raise Exception( 'Gamma already has MF=14 angular data: EGk=%s, ESk=%s, MT=%s' % ( EGk, ESk, MT ) )

    if( not( 14 in MTData ) ) : return
    MF14Data = MTData[14]
    ZA, AWR, LI, LTT, NK14, NI = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats( MF14Data[0], intIndices = [ 0, 2, 3, 4, 5 ], logFile = info.logs )
    printAWR_mode( info, MT, 14, 0, ZA, AWR )
    info.addMassAWR( ZA, AWR )
    if( ( NK14 != NK ) and info.printBadNK14 ) :
        warningList.append( 'MF14 NK = %s != MF12/13 NK = %s for MT = %s' % ( NK14, NK, MT ) )
    dataLine, frame = NI + 1, xDataEnumsModule.Frame.lab
    if( LTT == 0 ) :                                     # All distributions are isotropic
        pass
    elif( LTT == 1 ) :
        for i in range( NK14 - NI ) :
            dataLine, angularData = endfFileToGNDSMiscModule.getTAB2_Lists( dataLine, MF14Data, logFile = info.logs )
            if( angularData['NR'] != 1 ) : raise Exception( 'Currently only one interpolation flag is supported: NR=%s, MT=%s' % ( angularData['NR'], MT ) )
            EGk, ESk = angularData['C1'], angularData['C2']
            subform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 4, 'LTT = 1' )
            addAngularData( EGk, ESk, subform )
    else :
        raise Exception( 'MF=14, LI=%s not implemented' % LI )


def readMF15(info, MT, MTData, continuousGamma, warningList):

    if 15 not in MTData:
        if len(continuousGamma):
            warningList.append('Continous gamma with no MF=15 data: MT=%s' % MT)
            info.doRaise.append(warningList[-1])
        return
    if len(continuousGamma) == 0:
        warningList.append('MF=15 data and no continous gamma MF=12,13 data: MT=%s' % MT)
        info.doRaise.append(warningList[-1])
    MF15Data = MTData[15]
    ZA, AWR, dummy, dummy, NC, dummy = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats(
        MF15Data[0], intIndices=[0, 4], logFile=info.logs)
    printAWR_mode(info, MT, 15, 0, ZA, AWR)
    info.addMassAWR(ZA, AWR)
    if NC != 1: raise Exception('NC = %s > 1 is not supported: MT=%s' % (NC, MT))

    dataLine, weights = endfFileToGNDSMiscModule.getTAB1(1, MF15Data, logFile = info.logs)
    for Ein, weight in weights['data']:
        if weight != 1: raise Exception('For MF15 data weight other than 1 is not currently supportd: MT%d' % MT)

    LF = int(weights['L2'])
    if LF != 1: raise Exception('For MF=15 data, only LF=1 is currently supported, not LF=%s : MT%s' % (LF, MT))
    dataLine, EEpETable = endfFileToGNDSMiscModule.getTAB2_TAB1s(dataLine, MF15Data, logFile=info.logs, axes=energyAxes)
    continuousGamma[0]['energySubform'] = toPointwiseOrPiecewiseEnergy(MT, EEpETable)


def genID(cov_info, MT, MF, MT2=None, MF2=None, MAT2=None, QI=None, QI2=None, linkType='rowColumn'):
    """
    For covariances we need a unique id for each section, and also need link info.
    This is messy: lots of special cases.
    """

    evalStyle = cov_info['style']
    MTdict = cov_info['MTdict']

    def getReaction( MT, MF, QI=None ):
        # find the <reaction> corresponding to this covariance info, if any
        if MT in (452,455):
            for multiplicitySum in cov_info['multiplicitySums']:
                if multiplicitySum.ENDF_MT == MT:
                    return multiplicitySum
        if MT == 456:
            MT = 18 # redirect prompt nubar to fission
        if MT in MTdict:
            chThisMT = MTdict[MT]
            if len(chThisMT)==1:
                return chThisMT[0]
            elif MF==40:   # MF40 section may contain covariances for multiple excited states of residual
                thisChannel = [ch for ch in chThisMT if isinstance(ch, productionModule.Production) and ch.getQ('eV')==QI ]
                if len(thisChannel)==1: return thisChannel[0]
                elif MT in cov_info:
                    # Q-values disagree between MF10/MF40. Assume 1st covariance is for lowest-energy residual, etc
                    index = cov_info[MT]
                    cov_info[MT] += 1
                else:
                    index = 0; cov_info[MT] = 1
                if index >= len(thisChannel):
                    raise BadCovariance( "Can't find production reaction corresponding to MF40 MT%d QI=%f"
                            % (MT, QI) )
                return thisChannel[ index ]
            elif MF==33:
                # residual must be in ground state unless MT in 51-90, etc
                thisChannel = [ch for ch in chThisMT if (isinstance(ch, reactionModule.Reaction) and
                        sum( [p.getLevelAsFloat('eV') for p in ch.outputChannel] )==0) or
                        isinstance(ch, sumsModule.CrossSectionSum)]
                if len(thisChannel) != 1: raise BadCovariance("MF33 MT%d covariance doesn't correspond to any channel!"
                        % MT )
                return thisChannel[0]
            else :
                raise BadCovariance( "Can't determine which reaction this covariance (MF%d MT%d) corresponds to!" % ( MF, MT ) )
        return

    rowReaction = getReaction( MT, MF, QI )

    def makeID(MT, reaction):
        if reaction is not None:
            Id = str(reaction)
        elif MT in (452,455,456):
            Id = { 452 : totalToken, 455 : delayedToken, 456 : promptToken }[MT]
        elif MT in (1,3,4,103,104,105,106,107):
            Id = {1:'total', 3:'nonelastic', 4:'sum(z,n)', 103:'sum(z,p)', 104:'sum(z,d)',
                    105:'sum(z,t)', 106:'sum(z,He3)', 107:'sum(z,a)'}[MT]
        elif 850 < MT < 871:
            Id = "lump%d" % (MT-851)
        else:
            Id = "MF%d_MT%d" % (MF,MT)

        return Id

    rowId = makeID( MT, rowReaction )

    otherTarget = None
    if MAT2:
        # cross-material covariance
        try:
            otherTarget = endf_endlModule.getParticleNameFromMAT(MAT2)
        except KeyError:
            raise BadCovariance("Encountered illegal MAT number %d in covariances!" % MAT2)

        colReaction = None  # must create an 'externalReaction' and link to it below
        versus = ' vs. '
        colId = "%s MF%d MT%d" % (otherTarget,MF2,MT2)
    elif (MT2 and MT2!=MT):
        MF2 = MF2 or MF
        colReaction = getReaction( MT2, MF2, QI2 )
        versus = ' vs. '
        colId = makeID( MT2, colReaction )
    else:
        colId = versus = ''

    qualifier = {31:' [nubar]', 33:'', 34:' [angular distribution]', 35:' [spectrum]', 40:''}[MF]
    ID = rowId + versus + colId + qualifier

    # also create links from covariances to data:

    def makeLink( MT, reaction, Id, MAT2, linkClass ):
        link_ = linkClass(ENDF_MFMT="%d,%d" % (MF,MT))
        link_.root = "$reactions"
        if reaction is not None:
            if MF in (33,40):
                link_.link = reaction.crossSection[evalStyle]
            elif MF==31:
                if isinstance(reaction, sumsModule.MultiplicitySum):
                    link_.link = reaction.multiplicity[evalStyle]
                else:
                    link_.link = reaction.outputChannel[0].multiplicity[evalStyle]
            elif MF==34:
                link_.link = reaction.outputChannel[0].distribution[evalStyle].subforms[0]
            elif MF==35:
                distribution = reaction.outputChannel[0].distribution[evalStyle]
                if( isinstance( distribution, uncorrelatedModule.Form ) ) :
                    link_.link = distribution.energySubform
                else :
                    link_.link = distribution
        elif MAT2:
            # cross-material covariance: link needs to point to an external file
            cov_info['externalFiles'].append( otherTarget )
            link_.root = "$%s" % otherTarget
            link_.path = "/FIXME/path/to/MT%d" % MT2
        else:
            if 850 < MT < 871:
                pass
            elif MT == 4:
                cov_info['MTL_2'][(MT,MF)] =  [(smt, 33) for smt in range(50, 92)]
            elif MT == 103:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(600, 650)]
            elif MT == 104:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(650, 700)]
            elif MT == 105:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(700, 750)]
            elif MT == 106:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(750, 800)]
            elif MT == 107:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(800, 850)]
            elif MT == 102:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(900, 1000)]
            elif MT == 16:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(875, 892)]
            elif MT == 1:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(2, 1000)]
            elif MT == 3:
                cov_info['MTL_2'][(MT,MF)] = [(smt, 33) for smt in range(4, 1000)]
            else:
                message = f"Can't find reaction corresponding to MF{MF}/MT{MT} covariance. Treating as lumped sum."
                print("       WARNING:", message)
                cov_info['doRaise'].append(message)

            if (MT,MF) not in cov_info['lumpedChannels']:
                cov_info['lumpedChannels'][(MT,MF)] = sumsModule.CrossSectionSum( label=Id, ENDF_MT=MT )
            quant = cov_info['lumpedChannels'][(MT,MF)]

            link_.link = quant.crossSection

        return link_

    linkClass = covarianceSummedModule.Summand
    if linkType == 'rowColumn':
        linkClass = covarianceSectionModule.RowData
    rowData = makeLink( MT, rowReaction, rowId, MAT2, linkClass )
    colData = None
    if (MT2 and MT2!=MT) or MAT2:
        if linkType == 'rowColumn':
            linkClass = covarianceSectionModule.ColumnData
        colData = makeLink( MT2, colReaction, colId, MAT2, linkClass )
    return ID, [rowData, colData]


def readMatrix(info, MF, MT, LS,LB, NT,NP, dat, warningList):
    """ matrix format is very similar for MF31, MF33 and MF35, so we can
    generalize parsing the matrix """

    if NP in (0,1):   # matrix has at most one energy boundary, no data
        return None
    nlines, remainder = divmod(NT, 6)
    subsec = []
    for line in range(nlines): subsec += funkyF(dat.next(), logFile=info.logs)
    if remainder: subsec += funkyF(dat.next(), logFile=info.logs)[:remainder]
    # LB flag tells how to interpret this data:
    if LB in (0, 1, 2, 8, 9): # only diagonal is stored
        subsec = subsec[:NT]
        energyBounds = [subsec[::2]]
        data = subsec[1::2]
        if LB == 2:
            import numpy
            data = numpy.array(data)[:-1]
            negatives, = numpy.where(data < 0)
            if any(negatives):
                info.LB2_firstNegativeIndex = negatives[0]
            rawMatrix = numpy.outer(data, data)
            lowerDiagonal = rawMatrix[numpy.tril_indices(len(data))]
            matrix = arrayModule.Full(rawMatrix.shape, lowerDiagonal, symmetry=arrayModule.Symmetry.lower)
        else:
            if data[-1] != 0:
                warningList.append('Ignoring non-zero trailing element in LB=%d matrix for MF%d/MT%d' % (LB, MF, MT))
            data = data[:-1]
            matrix = arrayModule.Diagonal((len(data), len(data)), data)
    elif LB == 5:
        if LS == 1:  # symmetric upper-diagonal
            energyBounds = [subsec[:NP],]
            data = linearAlgebraModule.switchSymmetry(subsec[NP:NT])
            # ACD skip this, first row is first entry, second row is 2nd two, 3rd is next three
            # ACD parser will give list of numbers; maybe make numpy array that's 100 long and reshape to 10x10 matrix
            # ACD: arrayModule has containers that GNDS recognizes for covarianceSuite
            matrix = arrayModule.Full((NP-1, NP-1), data, symmetry=arrayModule.Symmetry.lower)
        else:  # asymmetric, but same energy grids for both axes:
            energyBounds = [subsec[:NP]]
            matrix = subsec[NP:]
            matrix = arrayModule.Full((NP-1,NP-1), matrix)
    elif LB == 6:  # asymmetric
        NER = NP
        NEC = (NT-1)//NER
        energyBounds = [subsec[:NER], subsec[NER:NER+NEC]]
        matrix = subsec[NER+NEC:NT]
        NEC -= 1; NER -= 1  # matrix dimensions
        matrix = arrayModule.Full((NER,NEC), matrix)
    else:
        return None

    for energies in energyBounds:
        if energies != sorted(energies):
            warningList.append('Energy boundaries out of order in MF%d MT%d' % (MF, MT))
            info.doRaise.append(warningList[-1])

    # FIXME: covariances need 'flat' interpolation along both independent axes, but gridded doesn't support that yet
    axes = axesModule.Axes(3, labelsUnits={
        0: ('matrix_elements', ''),
        1: ('column_energy_bounds', 'eV'),
        2: ('row_energy_bounds', 'eV')})
    axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit,
                              style=xDataEnumsModule.GridStyle.boundaries, values=valuesModule.Values(energyBounds[0]))
    if len(energyBounds) == 2:
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit,
                                  style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values(energyBounds[1]))
    else:
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit,
                                  style=xDataEnumsModule.GridStyle.boundaries,
                                  values=linkModule.Link(link=axes[2].values, relative=True))
# BRB FIX ME
    return griddedModule.Gridded2d(axes, matrix)


def readMF31_33(info, dat, mf, mt, cov_info, warningList):
    """ nubar and cross section covariances have basically the same form
    in ENDF, so we can treat them the same way: """

    dat = MyIter(dat)
    # dat contains one MFMT section, but in gnds each cross-correlation gets its own section:
    sectionList, linkData = [], []

    AWR_lineNumber = dat.index
    ZA, AWR, dum, MTL, dum, NL = funkyFI(dat.next(), logFile=info.logs)
    ZA = int(ZA)
    info.ZA_massLineInfo.add(ZA, AWR, mt, mf, AWR_lineNumber)
    info.addMassAWR(ZA, AWR)
    if MTL != 0:
        # MTL!=0 implies section is a placeholder pointing to a lumped channel
        if (MTL, mf) in cov_info['MTL']:
            cov_info['MTL'][(MTL, mf)].append((mt, mf))
        else: cov_info['MTL'][(MTL, mf)] = [(mt, mf)]
    for subsection in range(NL):
        XMF1, XLFS1, MAT1, MT1, NC, NI = funkyFI(dat.next(), logFile=info.logs)
        if not 0 < MT1 < 1000:
            print("    WARNING: invalid MT1 = %d" % MT1)
            # FIXME should this be an error? Evaluations occasionally use MT1=0 to mean MT1==mt
        if MAT1 == info.MAT:
            # Some files explicitly give MAT1==MAT for internal covariance.
            # for simplicity, just set MAT1=0 unless it is actually a different material
            MAT1 = 0
        XMF1, XLFS1 = int(XMF1), int(XLFS1)

        covarsThisSection = []

        if XMF1 != 0:     # field is usually 0, but could be 1 for MF31 or 3 for MF33.
            if XMF1 == mf - 30:
                print("    WARNING: XMF1 = %d, should be 0." % XMF1)
                XMF1 = 0
            else:
                raise BadCovariance("XMF1=%d for MF=%d in covariances not currently supported!" % (XMF1, mf))

        if XLFS1 != 0:
            raise BadCovariance("non-zero XLFS1 in covariances not currently supported!")

        for NCdx in range(NC):
            dum, dum, dum, LTY, dum, dum = funkyFI(dat.next(), logFile=info.logs)
            if LTY == 0:
                E1, E2, dum, dum, NCI2, NCI = funkyFI(dat.next(), logFile=info.logs)
                subsec = []
                nlines = int(math.ceil(NCI2/6.0))
                for line in range(nlines): subsec += funkyF(dat.next(), logFile=info.logs)
                summands = []
                for coef, mtnum in zip(subsec[:NCI2][::2], subsec[:NCI2][1::2]):
                    Id, pointers = genID(cov_info, int(mtnum), mf, linkType='summand')
                    link = pointers[0]
                    link.coefficient = coef
                    summands.append(link)
                covarsThisSection.append(covarianceSummedModule.SummedCovariance(
                    label=info.style, domainMin=E1, domainMax=E2, domainUnit='eV', summands=summands))
                cov_info['NC_data'].append(covarsThisSection[-1])
            else:
                warningList.append('non-zero LTY in MF33')

        for NIdx in range(NI):
            dum, dum, LS, LB, NT, NP = funkyFI(dat.next(), logFile=info.logs)
            matrix = readMatrix(info, mf,mt,LS,LB,NT,NP, dat, warningList)
            if LB not in (0, 1, 2, 5, 6, 8):
                warningList.append('skipping LB%d section for MF%d MT%d' % (LB, mf, mt ))
                continue
            if matrix is None:
                warningList.append('skipping empty matrix for MF%d MT%d' % (mf, mt))
                info.doRaise.append(warningList[-1])
                continue
            Type = covarianceEnumsModule.Type.relative
            if LB in (0, 8, 9):
                Type = covarianceEnumsModule.Type.absolute
                matrix.axes[0].unit = 'b**2'

            ENDFconversionFlag = None
            if LB in (2, 3, 4): ENDFconversionFlag = "LB=%d" % LB

            if LB == 2 and hasattr(info, 'LB2_firstNegativeIndex'):
                ENDFconversionFlag += ",firstNegativeIndex=%d" % info.LB2_firstNegativeIndex
                del info.LB2_firstNegativeIndex

            if LB in (8, 9):
                varianceDependence = {8: shortRangeSelfScalingVarianceModule.DependenceOnProcessedGroupWidth.inverse,
                                      9: shortRangeSelfScalingVarianceModule.DependenceOnProcessedGroupWidth.direct}[LB]
                covmatrix = shortRangeSelfScalingVarianceModule.ShortRangeSelfScalingVariance(
                    label=info.style, type=Type, dependenceOnProcessedGroupWidth=varianceDependence, matrix=matrix)
            else:
                covmatrix = covarianceMatrixModule.CovarianceMatrix(label=info.style, type=Type, matrix=matrix)
            if ENDFconversionFlag:
                info.ENDFconversionFlags.add(covmatrix, ENDFconversionFlag)
            covarsThisSection.append(covmatrix)

        # create unique id for each section:
        idNow, pointers = genID(cov_info, mt, mf, MT2=MT1, MF2=(XMF1 or mf), MAT2=MAT1)
        rowdat, coldat = pointers
        # ACD: instantiating a section that will go into covarianceSuite
        section = covarianceSectionModule.CovarianceSection(label=idNow, rowData=rowdat, columnData=coldat)

        if len(covarsThisSection) > 1:
            form = covarianceMixedModule.MixedForm(label=info.style, components=covarsThisSection)
            for idx in range(len(form)):
                form[idx].label = str(idx)
        elif len(covarsThisSection) == 1:
            form = covarsThisSection[0]
        else:
            #raise Exception("Encountered empty covariance section!!!")
            info.logs.write("Missing covariance data from section!")
            continue
        section.add(form)

        sectionList.append(section)
        linkData.append((mt, mf, MT1, XMF1, idNow))
        # end loop over NL subsections

    if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%d MT%d" % (mf, mt))
    return sectionList, linkData


def readMF32(info, dat, mf, mt, cov_info, warningList):
    # MF=32 resonance parameter covariances. Must be synchronized with MF=2 resonance parameters

    import numpy
    resonances = cov_info['resonances']

    def swaprows(matrix, i1, i2, nrows):
        # matrix rows may be out-of-order and need resorting
        rows = matrix[i1:i1+nrows].copy()
        matrix[i1:i1+nrows] = matrix[i2:i2+nrows]
        matrix[i2:i2+nrows] = rows
        cols = matrix[:, i1:i1+nrows].copy()
        matrix[:, i1:i1+nrows] = matrix[:, i2:i2+nrows]
        matrix[:, i2:i2+nrows] = cols

    def read_LCOMP1(dim, matrixSize, dat):
        data = []
        nLines, rem = divmod(matrixSize, 6)
        for i in range(nLines): data += funkyF(dat.next(), logFile=info.logs)
        if rem: data += funkyF(dat.next(), logFile=info.logs)[:rem]
        matrix = numpy.zeros((dim, dim))
        start, length = 0, dim
        for i in range(dim):
            # data stores upper-diagonal matrix. Symmetrize:
            matrix[i, i:] = matrix[i:, i] = data[start:start + length]
            start = start + length
            length = length - 1
        return matrix

    def read_LCOMP2_correlation(NNN, NM, NDIGIT, dat):
        """
        :param NNN: matrix dimension
        :param NM: number of lines to read from dat
        :param NDIGIT: Number of significant digits stored
        :param dat: open ENDF MF=32 section
        :return: numpy.array with shape (NNN,NNN) storing correlations
        """
        matrix = numpy.zeros((NNN, NNN))
        for idx in range(NM):
            row, col, vals = endfFileToGNDSMiscModule.readEndfINTG(dat.next(), NDIGIT)
            vals = vals[:row-col]
            # go to 0-based index:
            row -= 1
            col -= 1
            if row >= NNN or col >= row:
                raise BadCovariance("Matrix indices out of range for MF32 LCOMP=2 matrix: (%d,%d) vs. %d" % (row,col,NNN))
            matrix[row, col:col+len(vals)] = vals
        for idx in range(NNN):    # symmetrize
            matrix[idx, idx:] = matrix[idx:, idx]
        matrix[matrix > 0] += 0.5
        matrix[matrix < 0] -= 0.5
        matrix += numpy.eye(NNN, dtype="float") * 10**NDIGIT
        matrix /= float(10**NDIGIT)
        return matrix

    dat = MyIter(dat)
    AWR_lineNumber = dat.index
    ZA, AWR, dum, dum, NIS, dum = funkyFI(dat.next(), logFile=info.logs)
    ZA = int(ZA)
    info.ZA_massLineInfo.add(ZA, AWR, mt, mf, AWR_lineNumber)
    info.addMassAWR(ZA, AWR)
    if NIS!=1: raise BadCovariance("Can't handle multi-isotope file 32!")
    ZAI, ABN, dum, LFW, NER, dum = funkyFI(dat.next(), logFile=info.logs)

    sections = []
    for subsection in range(NER):
        EL, EH, LRU, LRF, NRO, NAPS = funkyFI(dat.next(), logFile=info.logs)
        if NRO!=0: raise BadCovariance("Can't handle non-zero NRO in MF32!")
        # format is determined mainly by combination of LCOMP and LRU/LRF
        if LRU == 1:  # resolved resonance covariance section
            ENDFconversionFlags = []
            if LRF in (1, 2, 3):  # Breit-Wigner and simplified Reich-Moore formats are similar
                if not hasattr(resonances.resolved, 'evaluated'):
                    warningList.append("Resonance covariance data for non-existent resonance region")
                    break
                if LRF in (1, 2):
                    mf2_elist = list(zip(
                        resonances.resolved.evaluated.resonanceParameters.table.getColumn('energy'),
                        resonances.resolved.evaluated.resonanceParameters.table.getColumn('neutronWidth'),
                        resonances.resolved.evaluated.resonanceParameters.table.getColumn('captureWidth')))
                else:
                    ENDFconversionFlags.append('LRF3')
                    elasticLabel, = [r.label for r in resonances.resolved.evaluated.resonanceReactions
                                     if r.ejectile == IDsPoPsModule.neutron]
                    captureLabel, = [r.label for r in resonances.resolved.evaluated.resonanceReactions
                                     if r.ejectile == IDsPoPsModule.photon]
                    mf2_elist = [[], [], []]
                    for spinGroup in resonances.resolved.evaluated.spinGroups:
                        mf2_elist[0].extend( spinGroup.resonanceParameters.table.getColumn('energy'))
                        mf2_elist[1].extend(spinGroup.resonanceParameters.table.getColumn(elasticLabel + ' width'))
                        mf2_elist[2].extend(spinGroup.resonanceParameters.table.getColumn(captureLabel + ' width'))
                    mf2_elist = list(zip(*mf2_elist))
                SPI, AP, dum, LCOMP, NLS, ISR = funkyFI(dat.next(), logFile=info.logs)
                DAP = []
                if ISR > 0:   # scattering radius uncertainty
                    if LRF in (1, 2):
                        dum, dap, dum, dum, dum, dum = funkyFI(dat.next(), logFile=info.logs)
                        DAP = [10*dap]
                    else:  # LRF==3
                        dum, dum, dum, dum, MLS, one = funkyFI(dat.next(), logFile=info.logs)
                        for idx in range(int(math.ceil(MLS/6.0))):
                            DAP.extend(funkyF(dat.next(), logFile=info.logs) )
                        DAP = [10*dap for dap in DAP[:MLS]]  # convert from 10*fm to fm
                if LCOMP == 0:        # internal correlations given for each resonance, no cross-resonance terms
                    ENDFconversionFlags.append('LCOMP=0')
                    mf32_resonances, mf32_covars = [], []
                    NRS = 0
                    for Lval in range(NLS):
                        AWRI_lineNumber = dat.index
                        AWRI, dum, L, dum, tmp, nrs_ = funkyFI(dat.next(), logFile=info.logs)
                        info.ZA_massLineInfo.add(ZA, AWRI, mt, mf, AWRI_lineNumber, column=0)
                        NRS += nrs_
                        for i in range(nrs_):
                            mf32_resonances.append(funkyF(dat.next(), logFile=info.logs))
                            mf32_covars.append(
                                funkyF(dat.next(), logFile=info.logs) + funkyF(dat.next(), logFile=info.logs))

                    dEsq, dNsq, dNdG, dGsq, dNdF, dGdF, dFsq, dJdN, dJdG, dJdF, dJsq, dum = zip(*mf32_covars)
                    MPAR = 3
                    if any(dFsq): MPAR = 4
                    if any(dJsq): raise BadCovariance("Encountered uncertainty on J in MF32!")
                    matrix = numpy.zeros((MPAR*NRS, MPAR*NRS))
                    for ridx in range(NRS):
                        matrix[ridx*MPAR, ridx*MPAR] = dEsq[ridx]
                        matrix[ridx*MPAR+1, ridx*MPAR+1] = dNsq[ridx]
                        matrix[ridx*MPAR+2, ridx*MPAR+1:ridx*MPAR+3] = [dNdG[ridx], dGsq[ridx]]
                        if MPAR == 4:
                            matrix[ridx*MPAR+3, ridx*MPAR+1:ridx*MPAR+4] = [dNdF[ridx], dGdF[ridx], dFsq[ridx]]
                    # symmetrize:
                    for ridx in range(MPAR*NRS):
                        matrix[ridx, ridx:] = matrix[ridx:, ridx]

                    mf32_elist = [(lis[0], lis[3], lis[4]) for lis in mf32_resonances]
                    nResonances = len(mf32_elist)
                    _type = covarianceEnumsModule.Type.absolute
                    matrixClass = arrayModule.Flattened

                elif LCOMP == 1:
                    AWRI_lineNumber = dat.index
                    AWRI, dum, dum, dum, NSRS, NLRS = funkyFI(dat.next(), logFile=info.logs)
                    info.ZA_massLineInfo.add(ZA, AWRI, mt, mf, AWRI_lineNumber, column=0)
                    dum, dum, MPAR, dum, tmp, NRB = funkyFI(dat.next(), logFile=info.logs)
                    dim = NRB * MPAR  # num. of resonances * num. parameters per resonance
                    matrixSize = dim * (dim+1) // 2
                    if matrixSize + 6*NRB != tmp:
                        raise BadCovariance("Incorrect dimension for the matrix!")

                    # resonances are listed again (redundant!):
                    mf32_resonances = [funkyF(dat.next(), logFile=info.logs) for i in range(NRB)]
                    if LRF in (1, 2):
                        mf32_elist = [(lis[0], lis[3], lis[4]) for lis in mf32_resonances]
                    else:
                        mf32_elist = [(lis[0], lis[2], lis[3]) for lis in mf32_resonances]
                    matrix = read_LCOMP1(dim, matrixSize, dat)
                    nResonances = len(mf32_elist)
                    _type = covarianceEnumsModule.Type.absolute
                    matrixClass = arrayModule.Full

                elif LCOMP == 2:
                    ENDFconversionFlags.append('LCOMP=2')
                    AWRI_lineNumber = dat.index
                    AWRI, QX, dum, LRX, tmp, NRSA = funkyFI(dat.next(), logFile=info.logs)
                    info.ZA_massLineInfo.add(ZA, AWRI, mt, mf, AWRI_lineNumber, column=0)
                    # resonance parameters + uncertainties:
                    mf32_resonances = [funkyF(dat.next(), logFile=info.logs) for i in range(NRSA*2)]
                    if LRF in (1, 2):
                        mf32_elist = [(lis[0], lis[3], lis[4]) for lis in mf32_resonances[::2]]
                    else:
                        mf32_elist = [(lis[0], lis[2], lis[3]) for lis in mf32_resonances[::2]]
                    # for LCOMP==2, off-diagonal terms are given as correlation matrix:
                    dum, dum, NDIGIT, NNN, NM, dum = funkyFI(dat.next(), logFile=info.logs)
                    MPAR = NNN//NRSA
                    diagonal = []
                    for idx in range(NRSA):
                        if LRF in (1, 2):
                            dE, dum, dum, dGammaN, dGammaG, dGammaF = mf32_resonances[2*idx+1]
                            diagonal.extend([dE, dGammaN,dGammaG])
                            if MPAR == 4: diagonal.append(dGammaF)
                        elif LRF == 3:
                            dE, dum, dGammaN, dGammaG, dGammaF1, dGammaF2 = mf32_resonances[2*idx+1]
                            diagonal.extend([dE, dGammaN, dGammaG])
                            if MPAR == 4: diagonal.extend([dGammaF1])
                            elif MPAR == 5: diagonal.extend([dGammaF1,dGammaF2])
                    if len(diagonal) != NNN:
                        raise BadCovariance("Incorrect dimensions for LCOMP=2 matrix! Expected NNN=%d, got %d." %
                                            (len(diagonal), NNN))

                    # off-diagonal parts of matrix are stored as sparse correlation matrix:
                    ENDFconversionFlags.append('NDIGIT=%d' % NDIGIT)
                    matrix = read_LCOMP2_correlation(NNN, NM, NDIGIT, dat)

                    # convert correlation -> absolute covariance matrix
                    matrix = matrix * numpy.outer(diagonal, diagonal)
                    nResonances = NRSA
                    _type = covarianceEnumsModule.Type.absolute
                    matrixClass = arrayModule.Flattened

                if LRF in (1, 2):
                    start = 0
                    MPAR += 3    # expand matrix with zeros to account for L,J and totalWidth columns
                    index = []
                    for ridx in range(nResonances):
                        index.extend([start+1, start+2, start+3])
                        start += MPAR
                    n_b = matrix.shape[0] + len(index)
                    dim = nResonances * MPAR
                    assert n_b == dim
                    not_index = numpy.array([k for k in range(n_b) if k not in index])
                    matrix2 = numpy.zeros((dim, dim))
                    matrix2[not_index.reshape(-1, 1), not_index] = matrix
                    matrix = matrix2
                    matrixClass = arrayModule.Flattened  # even if originally LCOMP=1

                # mf32 may not contain all resonances from mf2:
                mf2_elist_sorted = sorted(mf2_elist, key=lambda res: res[0])
                mf32_elist_sorted = sorted(mf32_elist, key=lambda res: res[0])
                if mf32_elist != mf32_elist_sorted or LCOMP == 0:
                    ENDFconversionFlags.append('sortByL')

                if not set(mf32_elist).issubset(mf2_elist):
                    onlyInMF32 = set(mf32_elist).difference(mf2_elist)
                    ndiffs = len(onlyInMF32)
                    warningList.append("MF32 resonance parameters differ for %d resonances. For example:" % ndiffs)
                    for mf32res in sorted(onlyInMF32):
                        # find closest match (by resonance energy) in MF=2
                        eres = mf32res[0]
                        for idx, mf2res in enumerate(mf2_elist_sorted):
                            if mf2res[0] >= eres: break
                        if abs(mf2_elist_sorted[idx-1][0] - eres) < abs(mf2res[0] - eres):
                            idx -= 1
                            mf2res = mf2_elist_sorted[idx-1]
                        warningList.append("    resonance #%d: MF2 = %s, MF32 = %s" % (idx, mf2res, mf32res))
                        if not info.verboseWarnings: break

                    raise BadCovariance("MF32 resonance parameters don't match MF2 parameters!")

                if len(mf2_elist) != len(mf32_elist):
                    dim = len(mf2_elist) * MPAR
                    matrix2 = numpy.zeros((dim, dim))
                    matrix2[:len(matrix), :len(matrix)] = matrix
                    matrix = matrix2
                    for mf2res in mf2_elist:
                        if mf2res not in mf32_elist:
                            mf32_elist.append(mf2res)
                    matrixClass = arrayModule.Flattened  # since some rows will be all 0

                if mf32_elist != mf2_elist or LCOMP == 0:
                    # rearrange order of MF32 resonances to match GNDS storage order
                    mf32_elist_extended = mf32_elist + [v for v in mf2_elist if v not in mf32_elist]
                    for i1 in range(len(mf2_elist)):
                        i2 = mf32_elist_extended.index(mf2_elist[i1])
                        if i2 != i1:
                            swaprows(matrix, MPAR * i1, MPAR * i2, MPAR)
                            # also need to swap values in elist2:
                            val = mf32_elist_extended[i1]
                            mf32_elist_extended[i1] = mf32_elist_extended[i2]
                            mf32_elist_extended[i2] = val

                if LRF == 3:  # also swap elastic and capture widths to follow LRF=7 convention
                    for i1 in range(len(mf2_elist)):
                        swaprows(matrix, MPAR * i1 + 1, MPAR * i1 + 2, 1)

                parameters = covarianceModelParametersModule.Parameters()
                parameter_values = []
                startIndex = 0
                if DAP:     # scattering radius uncertainty was specified. Expand matrix to include it:
                    if LRF in (1, 2):
                        scatteringRadius = resonances.resolved.evaluated.getScatteringRadius()
                    else:
                        scatteringRadius = resonances.scatteringRadius
                    parameters.add(covarianceModelParametersModule.ParameterLink(
                        label="scatteringRadius", root="$reactions", link=scatteringRadius,
                        matrixStartIndex=startIndex, nParameters=1))
                    parameter_values.append(scatteringRadius.evaluated.value)
                    startIndex += 1
                    if len(DAP) == 1:
                        DAP_matrix = numpy.array([[DAP[0]**2]])
                    else:
                        # L-dependent scattering radius uncertainty, only supported for LRF=3 or 7
                        RMatrix = resonances.resolved.evaluated
                        radius = scatteringRadius.evaluated
                        DAPS_per_L = []
                        uniqueId = 0
                        for lidx, DAPnow in enumerate(DAP[1:]):
                            DAPS_per_L.append([])
                            if not DAPnow: continue
                            for spinGroup in RMatrix.spinGroups:
                                for channel in spinGroup.channels:
                                    if RMatrix.resonanceReactions[channel.resonanceReaction].eliminated:
                                        continue
                                    if channel.L == lidx:
                                        if channel.scatteringRadius is None:
                                            channel.scatteringRadius = scatteringRadiusModule.ScatteringRadius(radius)
                                            channel.hardSphereRadius = scatteringRadiusModule.HardSphereRadius(radius)
                                        parameters.add(covarianceModelParametersModule.ParameterLink(
                                            label=f"scatteringRadius_{uniqueId}", root="$reactions",
                                            link=channel.scatteringRadius, matrixStartIndex=startIndex,
                                            nParameters=1))
                                        parameter_values.append(channel.scatteringRadius.evaluated.value)
                                        parameters.add(covarianceModelParametersModule.ParameterLink(
                                            label=f"hardSphereRadius_{uniqueId}", root="$reactions",
                                            link=channel.hardSphereRadius, matrixStartIndex=startIndex+1,
                                            nParameters=1))
                                        parameter_values.append(channel.hardSphereRadius.evaluated.value)
                                        DAPS_per_L[-1] += [DAPnow, DAPnow]
                                        uniqueId += 1
                                        startIndex += 2

                        # DAP for each L forms a fully correlated block in the parameter matrix
                        N_DAPs = sum([len(tmp) for tmp in DAPS_per_L]) + 1
                        DAP_matrix = numpy.zeros((N_DAPs, N_DAPs))
                        DAP_matrix[0,0] = DAP[0]**2
                        matrix_start_index = 1
                        for tmp in DAPS_per_L:
                            if tmp:
                                block = numpy.outer(tmp, tmp)
                                slice_ = slice(matrix_start_index, matrix_start_index+len(tmp))
                                DAP_matrix[slice_, slice_] = block
                                matrix_start_index += len(tmp)
                        assert matrix_start_index == len(DAP_matrix)

                    dim = len(DAP_matrix) + len(matrix)
                    new_matrix = numpy.zeros((dim, dim))
                    new_matrix[:len(DAP_matrix), :len(DAP_matrix)] = DAP_matrix
                    new_matrix[len(DAP_matrix):, len(DAP_matrix):] = matrix
                    matrix = new_matrix

                # store into GNDS:
                if LRF in (1, 2):
                    resData = resonances.resolved.evaluated.resonanceParameters.table
                    nParams = resData.nColumns * resData.nRows
                    for row in resData.data:
                        parameter_values += row
                    parameters.add(covarianceModelParametersModule.ParameterLink(
                        label="resonanceParameters", root="$reactions", link=resData,
                        matrixStartIndex=startIndex, nParameters=nParams))
                else:
                    # for RMatrix need links to each spinGroup
                    excess_rows = []
                    for spinGroup in resonances.resolved.evaluated:
                        sgTable = spinGroup.resonanceParameters.table
                        if sgTable in info.extraFissionWidths:
                            # 2nd fission width is all 0 for this spin group.
                            # We can drop that column unless we have non-zero covariance data
                            cIndex = sgTable.nColumns - 1
                            assert sgTable.columns[cIndex].name == 'fission width_2'
                            width2_rows = numpy.array(range(startIndex, startIndex+len(sgTable), cIndex))
                            if MPAR < 5 or (MPAR == 5 and not numpy.any(matrix[width2_rows])):
                                sgTable.removeColumn('fission width_2')
                                spinGroup.channels.pop(spinGroup.channels.labels()[-1])
                                if MPAR == 5:
                                    # mark corresponding rows / columns of the matrix for deletion
                                    excess_rows += list(width2_rows)
                            info.extraFissionWidths.remove(sgTable)
                        nParams = sgTable.nColumns * sgTable.nRows
                        for row in sgTable.data:
                            parameter_values += row
                        if nParams == 0: continue
                        parameters.add(covarianceModelParametersModule.ParameterLink(
                            label=spinGroup.label, link=sgTable, root="$reactions",
                            matrixStartIndex=startIndex, nParameters=nParams
                        ))
                        startIndex += nParams

                    if excess_rows:
                        # remove extra covariance matrix rows corresponding to 2nd fission width:
                        tmp = numpy.delete(matrix, excess_rows, axis=0)
                        matrix = numpy.delete(tmp, excess_rows, axis=1)

                assert len(matrix) == len(parameter_values)

                if _type is covarianceEnumsModule.Type.absolute:
                    parameter_values = numpy.array(parameter_values)
                    if not numpy.any(matrix[parameter_values == 0]):
                        # convert to relative matrix
                        parameter_values[parameter_values == 0] = 1
                        matrix = matrix / parameter_values / parameter_values[:,None]
                        _type = covarianceEnumsModule.Type.relative

                # switch to diagonal matrix if possible (much more compact):
                if numpy.all(matrix == (numpy.identity(len(matrix)) * matrix.diagonal())):
                    GNDSmatrix = arrayModule.Diagonal(shape=matrix.shape, data=matrix.diagonal())
                elif matrixClass is arrayModule.Flattened:
                    GNDSmatrix = arrayModule.Flattened.fromNumpyArray(matrix, symmetry=arrayModule.Symmetry.lower)
                else:
                    GNDSmatrix = arrayModule.Full(shape=matrix.shape, data=matrix[numpy.tril_indices(len(matrix))],
                                                  symmetry=arrayModule.Symmetry.lower)

                covmatrix = covarianceModelParametersModule.ParameterCovarianceMatrix(
                    info.style, GNDSmatrix, parameters, type=_type)
                if ENDFconversionFlags:
                    info.ENDFconversionFlags.add(covmatrix, ','.join(ENDFconversionFlags))

            elif LRF == 7:
                dum, dum, IFG, LCOMP, NJS, ISR = funkyFI(dat.next(), logFile=info.logs)
                assert IFG in (0, 1), f"IFG={IFG} not supported"
                if ISR > 0:
                    raise NotImplementedError("scattering radius uncertainty in MF32 LRF7")
                if LCOMP == 1:
                    AWRI_lineNumber = dat.index
                    AWRI, dum, dum, dum, NSRS, NLRS = funkyFI(dat.next(), logFile=info.logs)
                    info.ZA_massLineInfo.add(ZA, AWRI, mt, mf, AWRI_lineNumber, column=0)
                    dum, dum, NJSX, dum, dum, dum = funkyFI(dat.next(), logFile=info.logs)

                    for jdx in range(NJSX):
                        spinGroup = resonances.resolved.evaluated[jdx]
                        dum, dum, NCH, NRB, sixNX, NX = funkyFI(dat.next(), logFile=info.logs)
                        assert sixNX == 6*NX
                        resonanceParams = []
                        nlines = int(math.ceil((NCH + 1) / 6.0))  # Extra "1" is for the Eres column
                        for i in range(max(1,NRB)):
                            vals = []
                            for j in range(nlines):
                                vals += funkyF(dat.next(), logFile=info.logs)
                            if NRB>0: resonanceParams.append(vals[:NCH + 1])

                        # test for consistency with MF2
                        if spinGroup.resonanceParameters.table.data != resonanceParams:
                            raise BadCovariance(
                                "MF32 resonance parameters don't match MF2 parameters for spin group %d!" % jdx)

                    # rest of matrix:
                    dum, dum, dum, dum, N, NPARB = funkyFI(dat.next(), logFile=info.logs)
                    assert N == (NPARB*(NPARB+1))/2
                    matrix = read_LCOMP1(NPARB, N, dat)
                    _type = covarianceEnumsModule.Type.absolute
                    ENDFconversionFlags.append("LCOMP=1")

                elif LCOMP == 2:
                    dum, dum, NPP, NJSX, twelveNPP, twoNPP = funkyFI(dat.next(), logFile=info.logs)
                    assert (twoNPP == 2*NPP) and (twelveNPP == 12*NPP)
                    for idx in range(NPP):
                        # FIXME should check these against MF2 values:
                        MA, MB, ZA, ZB, IA, IB = funkyF(dat.next(), logFile=info.logs)
                        Q, PNT, SHF, MT, PA, PB = funkyF(dat.next(), logFile=info.logs)

                    if NJSX not in (NJS, 0):
                        warningList.append( "WARNING in MF=32: NJSX not consistent with NJS: %d vs %d!" % (NJSX, NJS))

                    allUncerts = []
                    for jdx in range(NJS):
                        spinGroup = resonances.resolved.evaluated[jdx]
                        AJ, PJ, dum, dum, sixNCH, NCH = funkyFI(dat.next(), logFile=info.logs)
                        for cidx in range(NCH):
                            # FIXME should also check these against MF2:
                            PPI, L, SCH, BND, APE, APT = funkyFI(dat.next(), logFile=info.logs)

                        dum, dum, dum, NRSA, twelveNX, NX = funkyFI(dat.next(), logFile=info.logs)
                        if twelveNX != 12*NX:
                            warningList.append("WARNING: incorrect LRF7 header, line %d" % dat.index)
                        if NRSA == 0: dat.next()   # skip empty line
                        resonanceParams = []
                        resonanceUncerts = []
                        nlines = int(math.ceil((NCH+1)/6.0))  # Extra "1" is for the Eres column
                        for i in range(NRSA):
                            vals = []
                            uncerts = []
                            for j in range(nlines):
                                vals += funkyF(dat.next(), logFile=info.logs)
                            for j in range(nlines):
                                uncerts += funkyF(dat.next(), logFile=info.logs)
                            resonanceParams.append(vals[:NCH+1])
                            resonanceUncerts.append(uncerts[:NCH+1])
                            allUncerts += uncerts[:NCH+1]

                        # now test for consistency with MF2
                        if spinGroup.resonanceParameters.table.data != resonanceParams:
                            raise BadCovariance(
                                "MF32 resonance parameters don't match MF2 parameters for spin group %d!" % jdx)
                        J, pi = translateENDFJpi(AJ, PJ)
                        if not J == spinGroup.spin and pi == spinGroup.parity:
                            raise BadCovariance("Inconsistent J/pi for MF2 / MF32 spin group %d" % jdx)

                    # correlations:
                    dum, dum, NDIGIT, NNN, NM, dum = funkyFI(dat.next(), logFile=info.logs)
                    matrix = read_LCOMP2_correlation(NNN, NM, NDIGIT, dat)

                    # now we can either add uncertainty columns to the parameter tables (and store correlation matrix),
                    # or convert to covariance matrix. For now do the latter
                    rsd = numpy.array(allUncerts)
                    matrix = matrix * numpy.outer(rsd, rsd)
                    _type = covarianceEnumsModule.Type.absolute
                    ENDFconversionFlags.append("LCOMP=2")
                    ENDFconversionFlags.append("NDIGIT=%d" % NDIGIT)
                else:
                    raise NotImplementedError("MF32 LRF=7 LCOMP=%d" % LCOMP)

                # store into GNDS (need links to each spinGroup)
                parameters = covarianceModelParametersModule.Parameters()
                parameter_values = []
                startIndex = 0
                for spinGroup in resonances.resolved.evaluated:
                    nParams = spinGroup.resonanceParameters.table.nColumns * spinGroup.resonanceParameters.table.nRows
                    for row in spinGroup.resonanceParameters.table.data:
                        parameter_values += row
                    if nParams == 0: continue
                    parameters.add(covarianceModelParametersModule.ParameterLink(
                        label=spinGroup.label, link=spinGroup.resonanceParameters.table, root="$reactions",
                        matrixStartIndex=startIndex, nParameters=nParams
                    ))
                    startIndex += nParams

                if _type is covarianceEnumsModule.Type.absolute:
                    parameter_values = numpy.array(parameter_values)
                    if not numpy.any(matrix[parameter_values == 0]):
                        # convert to relative matrix
                        parameter_values[parameter_values == 0] = 1
                        matrix /= numpy.outer(parameter_values, parameter_values)
                        _type = covarianceEnumsModule.Type.relative

                # switch to diagonal matrix if possible (much more compact):
                if numpy.all(matrix == (numpy.identity(len(matrix)) * matrix.diagonal())):
                    GNDSmatrix = arrayModule.Diagonal(shape=matrix.shape, data=matrix.diagonal())
                else:
                    GNDSmatrix = arrayModule.Flattened.fromNumpyArray(matrix, symmetry=arrayModule.Symmetry.lower)

                covmatrix = covarianceModelParametersModule.ParameterCovarianceMatrix(
                    info.style, GNDSmatrix, parameters, type=_type)
                if ENDFconversionFlags:
                    info.ENDFconversionFlags.add(covmatrix, ','.join(ENDFconversionFlags))

            else:
                raise KeyError("Unknown LRF %d encountered in MF32" % LRF)

            rowData = covarianceSectionModule.RowData(
                info.reactionSuite.resonances.resolved.evaluated, root='$reactions')
            parameterSection = covarianceModelParametersModule.ParameterCovariance("resolved resonances", rowData)
            parameterSection.add(covmatrix)
            sections.append(parameterSection)

        else:
            # unresolved resonance parameters

            def makeURRcovariance(uncert, energyBounds, conversionFlag=None):
                matrix = arrayModule.Full(shape=(1, 1), data=[uncert])

                axes = axesModule.Axes(3, labelsUnits={
                    0: ('matrix_elements', ''),
                    1: ('column_energy_bounds', 'eV'),
                    2: ('row_energy_bounds', 'eV')})
                axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit,
                                          style=xDataEnumsModule.GridStyle.boundaries,
                                          values=valuesModule.Values(energyBounds[0]))
                if len(energyBounds) == 2:
                    axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit,
                                              style=xDataEnumsModule.GridStyle.boundaries,
                                              values=valuesModule.Values(energyBounds[1]))
                else:
                    axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit,
                                              style=xDataEnumsModule.GridStyle.boundaries,
                                              values=linkModule.Link(link=axes[2].values, relative=True))
                covmatrix = covarianceMatrixModule.CovarianceMatrix(
                    info.style, type=covarianceEnumsModule.Type.relative, matrix=griddedModule.Gridded2d(axes, matrix))
                if conversionFlag:
                    info.ENDFconversionFlags.add(covmatrix, conversionFlag)
                return covmatrix

            URR = resonances.unresolved.evaluated
            LJs = []

            SPI, AP, dum, dum, NLS, dum = funkyFI(dat.next(), logFile=info.logs)
            for lval in range(NLS):
                AWRI_lineNumber = dat.index
                AWRI, dum, L, dum, tmp, NJS = funkyFI(dat.next(), logFile=info.logs)
                info.ZA_massLineInfo.add(ZA, AWRI, mt, mf, AWRI_lineNumber, column=0)
                if tmp != 6*NJS: raise BadCovariance("Incorrect header in MF32 unresolved section!")
                for jval in range(NJS):
                    D, AJ, GNO, GG, GF, GX = funkyF(dat.next(), logFile=info.logs)
                    if AJ.is_integer(): AJ = int(AJ)
                    LJs.append((L, AJ, {'D': D, 'GNO': GNO, 'GG': GG, 'GF': GF, 'GX': GX}))

            # matrix:
            dum, dum, MPAR, dum, tmp, NPAR = funkyFI(dat.next(), logFile=info.logs)
            if tmp != (NPAR*(NPAR+1))/2: raise BadCovariance("Incorrect header in MF32 unresolved section!")
            nlines = int(math.ceil(tmp/6.0))
            data = []
            for line in range(nlines): data += funkyF(dat.next(), logFile=info.logs)
            matrix = numpy.zeros((NPAR, NPAR))
            start, length = 0, NPAR
            for i1 in range(NPAR):
                matrix[i1, i1:] = matrix[i1:, i1] = data[start:start+length]
                start = start+length
                length = length-1
            if numpy.all(matrix == 0):
                warningList.append("ignoring empty unresolved covariance matrix!")
                continue

            # find URR section corresponding to each row in the matrix:
            matrixSections = []
            urr_warnings = []
            for L, J, averageParams in LJs:
                conversionFlag = []
                for key in sorted(averageParams):
                    value = averageParams[key]
                    if value: conversionFlag.append('%s=%s' % (key,value))
                conversionFlag = ','.join(conversionFlag)

                lsections = [lsec for lsec in URR.Ls if lsec.L == L]
                if len(lsections) == 0:
                    urr_warnings.append("No match in MF2 for MF32 URR section with L=%d" % L)
                    continue
                lsection, = lsections
                jsections = [jsec for jsec in lsection.Js if jsec.J == J]
                if len(jsections) == 0:
                    urr_warnings.append("No match in MF2 for MF32 URR section with L=%d, J=%d" % (L, J))
                    continue
                jsection, = jsections
                matrixSections.append([lsection.L, jsection.J, jsection.levelSpacing, conversionFlag])
                if MPAR > len(jsection.widths) + 1:
                    urr_warnings.append(f"MF32 MPAR={MPAR} is too large")
                    continue
                for idx in range(MPAR-1):
                    matrixSections.append([lsection.L, jsection.J, jsection.widths[idx], conversionFlag])

            if urr_warnings:
                print("Errors encountered in MF=32 URR covariance section:")
                for warning in urr_warnings:
                    print("    " + warning)
                warningList += urr_warnings
                info.doRaise.append("Encountered malformed MF=32 URR data")
                continue

            assert len(matrixSections) == len(matrix)

            crossTermCounter = 0
            for sidx, (L, J, width, conversionFlag) in enumerate(matrixSections):

                uncert = matrix[sidx, sidx]
                rowData = covarianceSectionModule.RowData(width, root='$reactions')
                if isinstance(width, unresolvedResonanceModule.LevelSpacing):
                    label = "URR levelSpacing: L=%s J=%s" % (L, float(J))
                else:
                    label = "%s URR: L=%s J=%s" % (width.resonanceReaction, L, float(J))
                covarianceSection = covarianceModelParametersModule.AverageParameterCovariance(
                    label, rowData=rowData )
                covarianceSection.add(makeURRcovariance(uncert, [width.data.domain()], conversionFlag))
                width.data.uncertainty = uncertaintiesModule.Uncertainty(
                            functional=uncertaintiesModule.Covariance(
                                link=covarianceSection[info.style], root="$covariances"))
                sections.append(covarianceSection)

                if numpy.any(matrix[sidx, sidx+1:]): # cross terms present

                    for ctidx, crossTerm in enumerate(matrix[sidx, sidx+1:]):
                        if crossTerm != 0:
                            otherWidth = matrixSections[sidx+ctidx+1][2]
                            columnData = covarianceSectionModule.ColumnData(otherWidth, root='$reactions')
                            label = "URR cross term %d" % crossTermCounter
                            covarianceSection = covarianceModelParametersModule.AverageParameterCovariance(
                                label, rowData=rowData, columnData=columnData)
                            covarianceSection.add(makeURRcovariance(crossTerm, [width.data.domain()]))
                            sections.append(covarianceSection)

                            crossTermCounter += 1
    return sections, []


def readMF34(info, dat, mf, mt, cov_info, warningList):
    """ angular distribution covariances: """

    # dat contains one MFMT section
    dat = MyIter(dat)
    sectionList, linkData = [], []

    AWR_lineNumber = dat.index
    ZA, AWR, dum, LTT, dum, NMT = funkyFI(dat.next(), logFile=info.logs)
    ZA = int(ZA)
    info.ZA_massLineInfo.add(ZA, AWR, mt, mf, AWR_lineNumber)
    info.addMassAWR(ZA, AWR)

    for subsection in range(NMT):
        dum, dum, MAT1, MT1, NL, NL1 = funkyFI(dat.next(), logFile=info.logs)
        if MT1 == mt:
            NSS = NL*(NL+1)//2
        else:
            NSS = NL*NL1
        if MAT1 != 0:
            raise BadCovariance(
                "Cross material angular distribution covariance is not allowed in ENDF format.  Found MAT1 = %s"
                % MAT1)
        if MT1 != mt:
            raise NotImplementedError(
                "Cross reaction covariances in angular distribution covariance data not supported")

        for iNSS in range(NSS):
            dum, dum, L, L1, LCT, NI = funkyFI(dat.next(), logFile=info.logs)
            frame = [xDataEnumsModule.Frame.lab, xDataEnumsModule.Frame.lab, xDataEnumsModule.Frame.centerOfMass][LCT]

            covarsThisL = []

            for NIdx in range(NI):
                dum, dum, LS, LB, NT, NE = funkyFI(dat.next(), logFile=info.logs)
                if (LS not in [1, 0]) and (LB < 7):
                    raise BadCovariance("Unexpected LS%d LB%d in MF34" % (LS, LB))
                # pretend it's MF33 file:
                matrix = readMatrix(info, mf, mt, LS, LB, NT, NE, dat, warningList)
                if LB == 0: raise Exception("LB=0 in Legendre covariances")
                covarsThisL.append(covarianceMatrixModule.CovarianceMatrix(
                    label=info.style, type=covarianceEnumsModule.Type.relative, matrix=matrix, productFrame=frame))

            if len(covarsThisL) > 1:
                sectionForm = covarianceMixedModule.MixedForm(label=info.style, components=covarsThisL)
                for idx in range(len(sectionForm)): sectionForm[idx].label = str(idx)
            elif len(covarsThisL) == 1:
                sectionForm = covarsThisL[0]

            # add unique id to the section:
            idNow, pointers = genID(cov_info, mt, mf)
            rowdat, coldat = pointers
            rowdat.dimension = 2    # axis corresponding to incident energy
            rowdat.slices.add(covarianceSectionModule.Slice(1, domainValue=L))  # 1 = axis corresponding to Legendre order
            if L1 != L:
                coldat = covarianceSectionModule.ColumnData(rowdat.link, root=rowdat.root, ENDF_MFMT="34,2")
                coldat.slices.add(covarianceSectionModule.Slice(1, domainValue=L1))
            idNow += " L%d vs. L%d" % (L, L1)
            section = covarianceSectionModule.CovarianceSection(label=idNow, rowData=rowdat, columnData=coldat)
            section.add(sectionForm)

            sectionList.append(section)
            linkData.append((mt, mf, mt, mf, idNow))

            if LCT == 0:
                info.ENDFconversionFlags.add(sectionForm, "LCT=0")
                cov_info.setdefault('MF34_missingFrames', {}).setdefault(mt, []).append(section)

        if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%d MT%d" % (mf,mt))

    return sectionList, linkData


def readMF35(info, dat, mf, mt, cov_info, warningList):
    """ spectra covariances are fairly simple: """

    # dat contains one MFMT section
    dat = MyIter(dat)
    sectionList, linkData = [], []

    AWR_lineNumber = dat.index
    ZA, AWR, dum, dum, NK, dum = funkyFI(dat.next(), logFile=info.logs)
    ZA = int(ZA)
    info.ZA_massLineInfo.add(ZA, AWR, mt, mf, AWR_lineNumber)
    info.addMassAWR(ZA, AWR)

    for subsection in range(NK):
        # each subsection contains matrix for a different incident energy interval.
        E1, E2, LS, LB, NT, NE = funkyFI(dat.next(), logFile=info.logs)
        if not (LS == 1 and LB == 7):
            raise BadCovariance("Unexpected LS%d LB%d in MF35" % (LS, LB))
        # pretend it's MF33 file:
        LS, LB = 1, 5

        matrix = readMatrix(info, mf, mt, LS, LB, NT, NE, dat, warningList)
        form = covarianceMatrixModule.CovarianceMatrix(
            label=info.style, type=covarianceEnumsModule.Type.absolute, matrix=matrix)

        # add unique id to the section:
        idNow, pointers = genID(cov_info, mt, mf)
        idNow += " energy range %d" % subsection
        rowdat, coldat = pointers
        rowdat.dimension = 1    # axis corresponding to outgoing energy
        rowdat.slices.add(covarianceSectionModule.Slice(2, domainUnit='eV', domainMin=E1, domainMax=E2))
        section = covarianceSectionModule.CovarianceSection(label=idNow, rowData=rowdat, columnData=coldat)
        section.add(form)

        sectionList.append(section)
        linkData.append((mt, mf, mt, mf, idNow))

        # end loop over NK subsections

    if dat.index != dat.length:
        raise BadCovariance("Not all covariance data converted, MF%d MT%d" % (mf, mt))
    return sectionList, linkData


def readMF40(info, dat, mf, mt, cov_info, warningList):
    """ production of radioactive isotopes. Also very similar to MF33 """

    dat = MyIter(dat)
    sectionList, linkData = [], []
    AWR_lineNumber = dat.index
    ZA, AWR, LIS, dum, NS, dum = funkyFI(dat.next(), logFile=info.logs)
    ZA = int(ZA)
    info.ZA_massLineInfo.add(ZA, AWR, mt, mf, AWR_lineNumber, LIS=LIS)
    info.addMassAWR(ZA, AWR)

    # each subsection represents different excited state of residual
    for subsection in range(NS):
        try:
            QM, QI, IZAP, LFS, dum, NL = funkyFI(dat.next(), logFile=info.logs)
        except StopIteration:
            warningList.append('MF40 MT%d lists incorrect number of subsections!' % mt)
            info.doRaise.append(warningList[-1])
            break

        ELFS = QM - QI
        if LFS == 0 and ELFS != 0 and mt != 18:
            warningList.append(
                "MF40 claims non-zero QM-QI = %s for the ground state, MT = %s. Setting QI=QM." % (ELFS, mt))
            info.doRaise.append(warningList[-1])
            QI = QM
        if mt == 18 and IZAP != -1:
            warningList.append("For MT=18 MF=40, expected IZAP=-1 but got %d" % IZAP)
            if not info.acceptBadMF10FissionZAP:
                info.doRaise.append(warningList[-1])
            return [], []

        for subsubsection in range(NL):
            # each subsubsection is a single matrix
            XMF1, XLFS1, MAT1, MT1, NC, NI = funkyFI(dat.next(), logFile=info.logs)
            XMF1, XLFS1 = int(XMF1), int(XLFS1)  # XLFS1: level index

            covarsThisSection = []
            if XMF1 not in (0, 10):
                raise BadCovariance("non-zero XMF1/XLFS1 in covariances not currently handled!")
            if MAT1 != 0:
                warningList.append("cross-material covariance with MAT=%d" % MAT1)

            for NCdx in range(NC):
                dum, dum, dum, LTY, dum, dum = funkyFI(dat.next(), logFile=info.logs)
                if LTY == 0:
                    E1, E2, dum, dum, NCI2, NCI = funkyFI(dat.next(), logFile=info.logs)
                    subsec = []
                    nlines = int(math.ceil(NCI2/6.0))
                    for line in range(nlines): subsec += funkyF(dat.next(), logFile=info.logs)
                    #coefs = subsec[:NCI2][::2]
                    summands = [
                            linkModule.Link('summand', genID( cov_info, int(mtnum), mf )[0],
                                            attributes={'ENDF_MFMT':"%d,%d"%(mf,mtnum), 'coefficient':coef})
                            for coef,mtnum in zip(subsec[:NCI2][::2],subsec[:NCI2][1::2])
                            ]
                    covarsThisSection.append(covarianceSummedModule.SummedCovariance(
                        label=info.style, domainMin=E1, domainMax=E2, domainUnit='eV', summands=summands))
                else:
                    warningList.append('non-zero LTY in MF40')

            for NIdx in range(NI):
                dum, dum, LS, LB, NT, NP = funkyFI(dat.next(), logFile=info.logs)
                matrix = readMatrix(info, mf,mt,LS,LB,NT,NP, dat, warningList)
                if LB not in (0, 1, 5, 6):
                    warningList.append('skipping LB%d section for MF%d MT%d' % (LB, mf, mt))
                    continue
                Type = covarianceEnumsModule.Type.relative
                if LB == 0:
                    Type = covarianceEnumsModule.Type.absolute
                    matrix.axes[0].unit = 'b**2'

                covarsThisSection.append(
                    covarianceMatrixModule.CovarianceMatrix(label=info.style, type=Type, matrix=matrix))

            if MAT1 != 0:
               continue

            # create unique id for each section:
            idNow, pointers = genID(cov_info, mt, mf, MT2=MT1, MF2=(XMF1 or mf), MAT2=MAT1, QI=QI)
            rowdat, coldat = pointers
            section = covarianceSectionModule.CovarianceSection(label=idNow, rowData=rowdat, columnData=coldat)

            if len(covarsThisSection) > 1:
                form = covarianceMixedModule.MixedForm(label=info.style, components=covarsThisSection)
                for idx in range(len(form)):
                    form[idx].label = str(idx)
            elif len(covarsThisSection) == 1:
                form = covarsThisSection[0]
            else:
                #raise Exception("Encountered empty covariance section!!!")
                info.logs.write("Missing covariance data from section!")
                continue
            section.add(form)

            sectionList.append(section)
            linkData.append((mt, mf, MT1, XMF1, idNow))
            # end loop over NL subsections

    if dat.index != dat.length:
        warningList.append("Not all covariance data converted for MF%d MT%d" % (mf, mt))
        info.doRaise.append(warningList[-1])
    return sectionList, linkData


def parseMF6FissionData(info, MT, MF6Data, fissionNeutronsAndGammasDataFromMF6, warningList):

    print("    WARNING: parseMF6FissionData function not complete.")
    dataLine = 0
    ZA, AWR, JP, LCT, NK, dummy = funkyFI(MF6Data[dataLine], logFile=info.logs)
    info.ZA_massLineInfo.add(ZA, AWR, MT, 6, 0)
    dataLine += 1

    JPP = JP // 10
    JPN = JP - 10 * JPP
    for particleIndex in range(NK):
        dataLine = parseMF6FissionParticle(
            info, dataLine, MT, MF6Data, JPN, JPP, fissionNeutronsAndGammasDataFromMF6, warningList)


def parseMF6FissionParticle(info, dataLine, MT, MF6Data, JPN, JPP, fissionNeutronsAndGammasDataFromMF6, warningList):
    # The work on this is incomplete.

    dataLine, productData, multiplicityRegions = endfFileToGNDSMiscModule.getTAB1Regions(
        dataLine, MF6Data, logFile=info.logs, axes=multiplicityAxes)
    ZAP, AWP = int(productData['C1']), productData['C2']
    LIP, LAW, NP = productData['L1'], productData['L2'], productData['NR']
    info.ZA_massLineInfo.add(ZAP, AWP, MT, 6, 0)

    energyDistribution = None
    JP = JPP
    if ZAP == 1: JP = JPN
    if JP != 1: raise ValueError('JP = %s not suppported' % JP)
    fissionNeutronsAndGammasDataFromMF6[ZAP].append([multiplicityRegions, energyDistribution])

    return dataLine


def fillRemainingProductsResidualForBreakup(info, decayChannel, lightIsotopeNames, breakupProducts, residualZA, crossSection):

    residualZA2 = residualZA
    for productName in lightIsotopeNames:
        if productName in breakupProducts:
            multiplicity = breakupProducts[productName]
            product = toGNDSMiscModule.newGNDSParticle(
                info, toGNDSMiscModule.getTypeNameGamma(info, productNameToZA[productName]),
                crossSection, multiplicity=multiplicity)
            decayChannel.products.add(decayChannel.products.uniqueLabel(product))
            if residualZA % 1000 > 0:
                residualZA2 -= multiplicity * productNameToZA[productName]
            else:
                residualZA2 -= multiplicity * (1000 * (productNameToZA[productName] // 1000))
    if residualZA2 != 0:
        decayChannel.products.add(decayChannel.products.uniqueLabel(
            toGNDSMiscModule.newGNDSParticle(info, toGNDSMiscModule.getTypeNameGamma(info, residualZA2), crossSection)))


def parseReaction(info, target, projectileZA, targetZA, MT, MTData, warningList,
                  parseCrossSectionOnly=False, channelProcess=None):
    """
    Translate all available data for the reaction with this MT.
    :return: tuple(crossSection, outputChannel, MFKeys)  where MFKeys contains MF numbers that remain untranslated (should be empty)
    """

    productList = []
    MFKeys = list( MTData.keys( ) )
    info.logs.write( '    %3d %s' % ( MT, sorted( MFKeys ) ) )
    LRProductZAs = None # If not None, special case to treat excitation level of residual for breakup. Cannot calculate 'level' energy here as not all masses are know.

    for MF in [8, 9, 10, 31, 32, 33, 34, 35, 40, 45]:
        if MF in MFKeys: MFKeys.remove(MF)

    if (MT == 3) and (3 not in MFKeys):     # Kludge, for ENDF files that for MT 3 have MF 12 and 13 but not MF 3 data.
        QM, QI, crossSection, LR, breakupProducts = 0, 0, None, 0, None
    else:
        QM, QI, crossSection, LR, breakupProducts = readMF3(info, MT, MTData[3], warningList)
        MFKeys.remove(3)
    if parseCrossSectionOnly:
        channel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process=channelProcess)
        channel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, QM, crossSection))
        return (crossSection, channel, MFKeys, LRProductZAs)

    neutronMFs = []
    fissionNeutronsAndGammasDataFromMF6 = { 0 : [], 1 : [] }        # Keys are ZAP (i.e., 0 for gammas and 1 for neutrons).
    for MF in [4, 5]:
        if MF in MFKeys: neutronMFs.append(MF)

    if neutronMFs != [] and 6 in MFKeys:
        if MT == 18:
            parseMF6FissionData(info, MT, MTData[6], fissionNeutronsAndGammasDataFromMF6, warningList)
            MFKeys.remove(6)
        else:
            raise Exception('MF 6 present and MF 4 and/or 5 present for MT %s: not allowed' % MT) # This should never happen!

    endfMTProductList = endf_endlModule.endfMTtoC_ProductLists[MT]
    compoundZA = calculateZA(targetZA, projectileZA, minus=False)
    lightIsotopeZAs = sorted([productNameToZA[product] for product in lightIsotopeNames])
    lightIsotopeZAsMultiplicity = {}
    for product in lightIsotopeNames:
        lightIsotopeZAsMultiplicity[productNameToZA[product]] = endfMTProductList.productCounts[product]

    if 4 in neutronMFs or (MT == 18 and neutronMFs == [5]): # MT == 18 and neutronMFs == [5] is a special case for bad data (g + Am241).
        ZAP = 1
        if MT not in [2, 5, 18, 19, 20, 21, 38]:                # Not elastic, fission or sumOfRemainingReactions.
            for product in lightIsotopeNames:
                if endfMTProductList.productCounts[product] > 0: break
            ZAP = productNameToZA[product]

    isUndefinedTwoBody = (102 < MT <= 107) or (MT in [91, 649, 699, 749, 799, 849, 999])
    isTwoBody = MT == 2 or 50 <= MT < 91 or ((600 <= MT < 849 or 900 <= MT < 999) and not isUndefinedTwoBody)

    levelIndex, decayChannel, twoBodyResidualZA = None, None, None
    if 50 <= MT < 91:
        levelIndex = MT - 50
    elif 600 <= MT < 649:
        levelIndex = MT - 600
    elif 650 <= MT < 699:
        levelIndex = MT - 650
    elif 700 <= MT < 749:
        levelIndex = MT - 700
    elif 750 <= MT < 799:
        levelIndex = MT - 750
    elif 800 <= MT < 849:
        levelIndex = MT - 800
    elif 875 <= MT < 891:
        levelIndex = MT - 875
    elif 900 <= MT < 999:
        levelIndex = MT - 900
    level = QM - QI                                                 # If level > 0., residual is in an excited state.
    if breakupProducts is not None:
        if isTwoBody:
            if MT not in [50, 600, 650, 700, 750, 800, 900]:
                level = -QI
        elif MT == 91:
            pass
        else :
            print('\nQM, QI', QM, QI)
            print(breakupProducts)
            raise NotImplementedError('breakup for MT %s is not supported' % MT)
    elif MT in [50, 650, 700, 750, 800, 850, 900]:
        level = 0.0

    if isTwoBody or isUndefinedTwoBody:
        if MT == 2:
            ZAP = projectileZA
        else:
            for productName in endfMTProductList.productCounts :
                if endfMTProductList.productCounts[productName] != 0: break
            if productName == IDsPoPsModule.photon:
                ZAP = 0
            else:
                ZAP = productNameToZA[productName]
        twoBodyResidualZA = calculateZA(compoundZA, ZAP)
    undefinedLevelInfo = { 'ZA' : twoBodyResidualZA, 'level' : level, 'levelIndex' : levelIndex, 'count' : 0 }
    if neutronMFs == [4]:                     # This is a two-body reaction with only angular data.
        if not isTwoBody:
            raise ValueError('With only MF = 4 data present, reaction is assumed to be two-body and it is not for MT = %s' % MT)
        product = toGNDSMiscModule.newGNDSParticle(info, toGNDSMiscModule.getTypeNameGamma(info, ZAP), crossSection)
        form = readMF4(info, product, MT, MTData[4], angularModule.TwoBody, warningList)
        MFKeys.remove(4)
        productList.append(product)
    elif neutronMFs == [4, 5] or (neutronMFs == [5] and ZAP == 1):
        # Don't check ZAP if MT=5. Currently this combination, MT=5, MF=4/5 appears only for incident gammas
        if MT != 5 and ZAP != 1: raise ValueError('ZAP = %d != 1 for MFs = [ 4, 5 ] for MT = %d' % (ZAP, MT))
        multiplicity = 1
        if MT not in [2, 5, 18, 19, 20, 21, 38]:                # Not elastic or fission.
            for product in lightIsotopeNames :
                if endfMTProductList.productCounts[product] > 0: break
            ZAP = productNameToZA[product]
            multiplicity = endfMTProductList.productCounts[product]
        else:
            if MT not in [2, 5]: multiplicity = -1            # MT 5 is special case where (g,n') put into MT 5 instead of one of 50-91.
        product = toGNDSMiscModule.newGNDSParticle(info, toGNDSMiscModule.getTypeNameENDF(info, ZAP, undefinedLevelInfo),
                crossSection, multiplicity = multiplicity)

        if neutronMFs == [5]:
            warningList.append("MF=5 found with no accompanying MF=4, assuming angular distribution for MT=%i is isotropic"%MT)
            angularSubform = angularModule.Isotropic2d()             # MF = 5 data is always in lab frame.
        else:
            angularSubform = readMF4(info, product, MT, MTData[4], None, warningList)
            MFKeys.remove(4)

        energySubform , weights = readMF5(info, MT, MTData[5], warningList, product=product)
        MFKeys.remove(5)

        form = uncorrelated(info.style, frames[1], angularSubform, energySubform)  # BRB: is frame right.
        product.distribution.add(form)
        productList.append(product)
    elif 6 in MFKeys:
        compoundZA2 = compoundZA
        isTwoBody = readMF6(MT, info, MTData[6], productList, warningList, undefinedLevelInfo, isTwoBody, crossSection, LR, compoundZA)
        MFKeys.remove(6)
        if isTwoBody and MT == 102:
            ZAP = 0
            undefinedLevelInfo['ZA'] = calculateZA(compoundZA, ZAP)
    elif neutronMFs == []:
        if isTwoBody and False:                 # ????????? Why False
            raise Exception( 'How did we get here.' )
            product = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeNameGamma( info, ZAP ), crossSection )
            residualZA = calculateZA( compoundZA, ZAP )
            if( levelIndex is not None ) :
                if( ( levelIndex <= info.targetLevel ) and ( info.targetZA == residualZA ) ) : levelIndex -= 1
            if( QI != QM ) :     # Residual is in an excited state.
                decayChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                decayChannel.products.add( decayChannel.products.uniqueLabel(
                        toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeNameGamma( info, residualZA ), crossSection ) ) )
            residual = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeNameGamma( info, residualZA, level = level ),
                    crossSection, outputChannel = decayChannel )
            productList.append( product )
            productList.append( residual )
    else:
        pass

    _dummyCrossSection = []
    readMF12_13(info, MT, MTData, productList, warningList, crossSection, _dummyCrossSection)
    _crossSection = crossSection
    if _crossSection is None: _crossSection = _dummyCrossSection[0] # Should only happen for MT=3 with no MF=3.

    for MF in [12, 13, 14, 15]:
        if MF in MFKeys: MFKeys.remove(MF)

    if MT == 5:
        if QM != QI: info.logs.write('    --QM %s != QI = %s\n' % (QM, QI))
        outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.sumOfRemainingOutputChannels)
        outputChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, QM, _crossSection))
    elif MT == 102 and not isTwoBody:
        residualIndex, gammaMissing = -1, False
        for index, product in enumerate(productList):
            if product.pid != IDsPoPsModule.photon: residualIndex = index
            gammaMissing = product.pid == IDsPoPsModule.photon or gammaMissing
        if residualIndex == -1:
            productList.insert(0, toGNDSMiscModule.newGNDSParticle(info,
                    toGNDSMiscModule.getTypeNameENDF(info, calculateZA(compoundZA, 0), undefinedLevelInfo), crossSection))
        if residualIndex > 0: productList.insert(0, productList.pop(residualIndex))
        if not gammaMissing: productList.append(toGNDSMiscModule.newGNDSParticle(info,
                toGNDSMiscModule.getTypeNameENDF(info, 0, undefinedLevelInfo), crossSection))
        outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process=channelProcess)
        outputChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, QM, _crossSection))  # Q????? What about QI?
    elif isTwoBody:
        if QI == 0 and QM != 0: raise Exception("QI = 0, QM = %f for MT=%d" % (QM, MT))
        outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.twoBody, process=channelProcess)
        outputChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, QI, _crossSection))
        if len( [p for p in productList if p.pid != IDsPoPsModule.photon] ) == 0:
            gammaProducts = productList
            productList = []
            for ZA in lightIsotopeZAs :
                if lightIsotopeZAsMultiplicity[ZA] != 0:
                    productList.append(toGNDSMiscModule.newGNDSParticle(
                        info, toGNDSMiscModule.getTypeNameENDF(info, ZA, undefinedLevelInfo), crossSection))
                    break
            if len(productList) == 0:
                if MT != 2: raise Exception("product data for reaction MT = %s needs to be implemented" % MT)
                productList.append(toGNDSMiscModule.newGNDSParticle(
                        info, toGNDSMiscModule.getTypeNameENDF(info, projectileZA, undefinedLevelInfo), crossSection))
            productList += gammaProducts    # later logic assumes gammas are at the end of the list

        if LR == 1: LRProductZAs = [ particleZA( info, product.pid ) for product in productList ]

        decayProductList = productList[1:]
        productList = productList[:1]                               # Assume first product is "b" in "a + A -> b + B" where B is the larger product.
        ZA = particleZA(info, productList[0].pid)
        residualZA = calculateZA(compoundZA, ZA)

        if LR == 1:
            if MT not in [50, 600, 650, 700, 750, 800, 850, 900] and residualZA not in LRProductZAs:
                LRProductZAs = [ZAP, residualZA]
            else:
                LRProductZAs = None

        levelIndex = undefinedLevelInfo['levelIndex']
        if levelIndex is not None:
            if levelIndex <= info.targetLevel and targetZA == residualZA:
                levelIndex -= 1
        undefinedLevelInfo['levelIndex'] = levelIndex

        for index, product in enumerate(decayProductList):
            ZA = particleZA(info, product.pid)
            if residualZA == ZA:
                productList.append(decayProductList.pop(index))
                break
        if len(productList) < 2:
            if MT == 2:
                productList.append(toGNDSMiscModule.newGNDSParticle(info, target, crossSection))
            else:
                if ZA == undefinedLevelInfo['ZA']:
                    undefinedLevelInfo['ZA'] = None
                productList.append(toGNDSMiscModule.newGNDSParticle(info,
                        toGNDSMiscModule.getTypeNameENDF(info, residualZA, undefinedLevelInfo), crossSection))
            if 4 not in MTData.keys():
                info.ENDFconversionFlags.add(productList[-1], 'implicitProduct')
            if info.style in productList[0].distribution:
                recoilForm = angularModule.TwoBody(info.style, xDataEnumsModule.Frame.centerOfMass,
                        angularSubform = angularModule.Recoil( productList[0].distribution[info.style], relative=True))
                productList[-1].distribution.add(recoilForm)

        decayZAs, decayGammaList, decayNonGammaList = 0, [], []
        for decayProduct in decayProductList:
            if decayProduct.pid == IDsPoPsModule.photon:
                decayGammaList.append(decayProduct)
                mult = 1
            else:
                decayNonGammaList.append(decayProduct)
                mult = decayProduct.multiplicity.getConstant()
            decayZAs += particleZA(info, decayProduct.pid) * mult
        if LR == 1:
            if decayZAs != residualZA: raise Exception("decayZAs = %d != residualZA = %d" % (decayZAs, residualZA))
        elif decayZAs == 0:
            if len(decayGammaList) != 0:
                if len(decayNonGammaList) == 0:
                    decayNonGammaList.append(toGNDSMiscModule.newGNDSParticle(
                        info, toGNDSMiscModule.getTypeNameENDF(info, residualZA, None), crossSection))
            elif len(decayNonGammaList) != 0:
                if len(decayGammaList) == 0: decayGammaList.append(toGNDSMiscModule.newGNDSParticle(
                    info, toGNDSMiscModule.getTypeNameENDF(info, 0, None), crossSection))
            decayProductList = decayNonGammaList + decayGammaList
        else:
            raise Exception("decayZAs = %d != 0" % decayZAs)

        if breakupProducts is not None:
            if decayChannel is not None: raise Exception('breakupProducts and decayChannel both not None')
            decayChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
            decayChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, QM - QI, _crossSection))
            fillRemainingProductsResidualForBreakup(info, decayChannel, lightIsotopeNames, breakupProducts,
                particleZA(info, productList[1].pid), crossSection)
            productList[1].addOutputChannel(decayChannel)
        elif len(decayProductList) > 0:                         # At this point, both two bodies are in productList and second one is redisual.
            if QI > QM: raise Exception("Negative decay Q-value for MT%d, QI = %s, QM = %s" % (MT, QI, QM))
            decayChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
            decayChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, QM - QI, _crossSection))  # Q????? Not right?
            for decayProduct in decayProductList: decayChannel.products.add(decayChannel.products.uniqueLabel(decayProduct))
            productList[1].addOutputChannel(decayChannel)

    elif endfMTProductList.isFission:
        useThisQM = QM
        outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process=channelProcess)
        if hasattr(info, 'fissionEnergyReleaseData') and MT == 18:
                FER = getFissionEnergies(info, crossSection.domainMin, crossSection.domainMax, warningList)
                outputChannel.fissionFragmentData.fissionEnergyReleases.add(FER)

                # Check for consistency between polynomial expansion and approximate constant Q.
                ENDF_Q = FER.nonNeutrinoEnergy.data.evaluate(1e-5)
                if abs(QM - ENDF_Q) > 1e-7 * abs(QM): warningList.append("Fission QM inconsistent with energy release data for MT = " + str(MT))

                # Compute 'prompt' Q-value (energy to neutrons + gammas + prompt products) for GNDS:
                useThisQM = FER.promptProductKE.data.evaluate(1e-5) + FER.promptNeutronKE.data.evaluate(1e-5) + FER.promptGammaEnergy.data.evaluate(1e-5)
        outputChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, useThisQM, _crossSection))
        if MT == 18:
            if len(productList) > 0: outputChannel.products.add(outputChannel.products.uniqueLabel(productList.pop(0)))
            if len(outputChannel) == 0:
                multiplicity = multiplicityModule.Unspecified(info.style)
                product = toGNDSMiscModule.newGNDSParticle(info, toGNDSMiscModule.getTypeNameGamma(info, 1), crossSection,
                        multiplicity=multiplicity)
                outputChannel.products.add(outputChannel.products.uniqueLabel(product))
            else:
                for product in outputChannel:
                    if product.pid == IDsPoPsModule.neutron: break
                info.firstFissionNeutron = product
                if promptToken in info.totalOrPromptFissionNeutrons:
                    product.multiplicity.remove(info.style)
                    product.multiplicity.add(info.totalOrPromptFissionNeutrons[promptToken])
                    product.multiplicity.remove('constant')
                elif totalToken in info.totalOrPromptFissionNeutrons:
                    product.multiplicity.remove(info.style)
                    product.multiplicity.add(info.totalOrPromptFissionNeutrons[totalToken])
                    product.multiplicity.remove('constant')
                if hasattr(info, 'delayedFissionDecayChannel'):
                    for i1, (decayRate, _delayedNeutron) in enumerate(info.delayedFissionDecayChannel):
                        product = delayedNeutronModule.Product(IDsPoPsModule.neutron, label=IDsPoPsModule.neutron)
                        product.multiplicity.add(_delayedNeutron.multiplicity[info.style])
                        if len(_delayedNeutron.distribution) > 0:
                            product.distribution.add(_delayedNeutron.distribution[info.style])
                        else:
                            print('FIXME: need delayed neutron distribution')

                        delayedNeutron = delayedNeutronModule.DelayedNeutron(str(i1), product)
                        delayedNeutron.rate.add(rateModule.Double(info.style, decayRate, '1/s'))
                        outputChannel.fissionFragmentData.delayedNeutrons.add(delayedNeutron)
        else:
            if neutronMFs == []:
                pass    # we used to add a reference to total fission nubar and PFNS, but that's not physically correct
                """
                if hasattr(info, 'firstFissionNeutron'):
                    multiplicity = multiplicityModule.Reference(link=info.firstFissionNeutron.multiplicity, label=info.style)
                else:                                              # When singleMTOnly is fission MT != 18.
                    multiplicity = multiplicityModule.Unspecified(label=info.style)
                product = toGNDSMiscModule.newGNDSParticle(info, toGNDSMiscModule.getTypeNameGamma(info, 1),
                        crossSection, multiplicity=multiplicity)
                if hasattr(info, 'firstFissionNeutron'):
                    form = referenceModule.Form(link=info.firstFissionNeutron.distribution, label=info.style)
                    product.distribution.add(form)
                outputChannel.products.add(outputChannel.products.uniqueLabel(product))
                """

        # Some files have distributions for 1stChanceFission etc, but should still link to total nubar:
        for product in productList:
            multiplicity = product.multiplicity[info.style]
            if isinstance(multiplicity, multiplicityModule.Constant1d) and product.multiplicity[info.style].evaluate(0) == -1:
                if hasattr(info, 'firstFissionNeutron'):
                    product.multiplicity.remove(info.style)
                    multiplicity = multiplicityModule.Reference(info.firstFissionNeutron.multiplicity, label=info.style)
                    product.multiplicity.add(multiplicity)

        while len(productList) > 0: outputChannel.products.add(outputChannel.products.uniqueLabel(productList.pop(0)))
    else:
        Q = QI
        if isUndefinedTwoBody: Q = QM
        outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process=channelProcess)
        outputChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, Q, _crossSection))
        if isUndefinedTwoBody: info.ENDFconversionFlags.add(outputChannel.Q, "QI=%s" % QI)

        if MT not in [1, 18, 19, 20, 21, 38] + list(range(201,208)):
            residualZA, ZAsMultiplicities, productAsResidual, biggestProduct = compoundZA, {}, None, 0
            for index, product in enumerate(productList):
                if product.pid == IDsPoPsModule.photon: continue
                ZA = particleZA(info, product.pid)
                multiplicity = product.multiplicity[info.style]
                if isinstance(multiplicity, multiplicityModule.Constant1d):
                    mult = int(multiplicity.value)
                else:
                    info.logs.write('\n\nIncorrect multiplicity in ENDF file! MT = %s\n' % MT)
                    info.logs.write('Multiplicity should be constant but is "%s".\n' % multiplicity.moniker)
                    raise ValueError('Multiplicity should be a constant and it is not.')
                if ZA in lightIsotopeZAs:
                    residualZA = calculateZA(residualZA, mult * ZA, minus=True)
                        # If we have different distributions for both neutrons in (n,2n), n shows up twice in the productList.
                    if ZA in ZAsMultiplicities:
                        ZAsMultiplicities[ZA] += mult
                    else:
                        ZAsMultiplicities[ZA] = mult
                else:
                    if productAsResidual is not None:
                        raise Exception('multiple residuals for MT = %d, %s %s' % (MT, productAsResidual.pid, product.pid))
                    productAsResidual = product

            if residualZA != 0:
                for ZA in lightIsotopeZAs:
                    if ZA not in ZAsMultiplicities: ZAsMultiplicities[ZA] = 0
                    if ZAsMultiplicities[ZA] == lightIsotopeZAsMultiplicity[ZA]: continue       # All this ZA accounted for.
                    if ZAsMultiplicities[ZA] > lightIsotopeZAsMultiplicity[ZA]:
                        raise Exception('negative multiplicity for ZA = %s for MT = %s' % (ZA, MT))
                    multiplicity = lightIsotopeZAsMultiplicity[ZA] - ZAsMultiplicities[ZA]
                    productList.append(toGNDSMiscModule.newGNDSParticle(
                        info, toGNDSMiscModule.getTypeNameENDF(info, ZA, None), crossSection, multiplicity=multiplicity))
                    residualZA = calculateZA(residualZA, multiplicity * ZA, minus=True)
                if productAsResidual is None:
                    if residualZA > 0:
                        productList.append( toGNDSMiscModule.newGNDSParticle(
                            info, toGNDSMiscModule.getTypeNameENDF(info, residualZA, undefinedLevelInfo), _crossSection))

            if breakupProducts is not None:
                if MT == 91:
                    if decayChannel is not None: raise Exception('breakupProducts and decayChannel both not None')
                    decayChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                    decayChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, QM - QI, _crossSection))
                    fillRemainingProductsResidualForBreakup(info, decayChannel, lightIsotopeNames, breakupProducts,
                        particleZA(info, productList[1].pid), crossSection)
                    productList[1].addOutputChannel(decayChannel)
                else:
                    raise Exception('breakup not supported for MT %d' % MT)

    for product in productList:
        if len( product.distribution ) == 0:
            frame = xDataEnumsModule.Frame.lab
            if isTwoBody:
                frame = xDataEnumsModule.Frame.centerOfMass
            form = unspecifiedModule.Form(info.style, productFrame=frame)
            product.distribution.add(form)
            info.ENDFconversionFlags.add(product, 'implicitProduct')
        outputChannel.products.add(outputChannel.products.uniqueLabel(product))

    if outputChannel.genre == enumsModule.Genre.twoBody:
        productID = outputChannel[0].pid
        if productID != info.projectile:
            if productID not in info.missingTwoBodyMasses: info.missingTwoBodyMasses[productID] = []
            info.missingTwoBodyMasses[productID].append([outputChannel[1].pid, QM])

    LRProducts = None
    if LRProductZAs is not None: LRProducts = QI

    return crossSection, outputChannel, MFKeys, LRProducts


def parseCovariances(info, MTDatas, MTdict, singleMTOnly=None, resonances=None,
                     formatVersion=GNDS_formatVersionModule.default, verbose=1):

    evaluation = "%s-%d.%d" % (info.evaluation, info.NVER, info.LREL)
    covarianceSuite = covarianceSuiteModule.CovarianceSuite(
        info.projectile, info.target, evaluation, formatVersion=formatVersion,
        interaction=enumsModule.Interaction.nuclear)
    linkData = []    # which mf/mts need to be linked for each covariance?
    if singleMTOnly: return covarianceSuite, linkData

    # make list of available covariance information:
    warningList = []
    cov_info = {'MTL': {}, 'MTL_2': {}, 'lumpedChannels': {}, 'externalFiles': [], 'mfmts': [], 'MTdict': MTdict,
                'resonances': resonances, 'NC_data': [], 'style': info.style, 'projectile': info.projectile,
                'multiplicitySums': info.reactionSuite.sums.multiplicitySums, 'doRaise': info.doRaise}
    for mt in MTDatas.keys():
        if (singleMTOnly is not None) and (mt!=singleMTOnly): continue
        for mf in MTDatas[mt].keys():
            if mf >= 30:    # all covariance-type data
                cov_info['mfmts'].append((mf,mt))
    cov_info['mfmts'].sort()  # sorting first by MF, then by MT

    for mf,mt in cov_info['mfmts']:
        try:
            if mf in (31,33):
                covars, tmp = readMF31_33( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 32:
                covars, tmp = readMF32( info, MTDatas[151][32], mf, mt, cov_info, warningList )
            elif mf == 34:
                covars, tmp = readMF34( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 35:
                covars, tmp = readMF35( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 40:
                covars, tmp = readMF40( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            else:
                warningList.append( 'MF%d not yet supported' % mf)
                continue
            for cov in covars:
                if mf == 32: covarianceSuite.parameterCovariances.add(cov)
                else:
                    covarianceSuite.covarianceSections.add(cov)
                    if cov.columnData is None:
                        covarianceLink = uncertaintiesModule.Covariance(link=cov[info.style], root="$covariances")
                        centralValue = cov.rowData.link
                        if mt in range(851,872):
                            continue    # lumped channel covariances are handled below
                        if( mf == 35 ) : centralValue = centralValue.data
                        if not hasattr(centralValue, 'uncertainty'):
                            continue   # kludge for ENDF-VIII Ni58 and Ni60 (they have MF33 MT3 but no MF3 MT3)
                        if centralValue.uncertainty.data is not None:
                            if isinstance( centralValue.uncertainty.data, uncertaintiesModule.Covariance ):
                                listOfCovars = uncertaintiesModule.ListOfCovariances()
                                listOfCovars.add( centralValue.uncertainty.data )
                                centralValue.uncertainty = uncertaintiesModule.Uncertainty( functional=listOfCovars )
                            centralValue.uncertainty.data.add( covarianceLink )
                        else:
                            centralValue.uncertainty = uncertaintiesModule.Uncertainty( functional=covarianceLink )
            linkData += tmp
        except BadCovariance as e:
            warningList.append('MF%d MT%d covariance conversion failed with message "%s"' % (mf,mt,e) )
            info.doRaise.append( warningList[-1] )

    # fix links for summed matrices:
    for summedMatrix in cov_info['NC_data']:
        for pointer in summedMatrix:
            pointed_to = [sec[info.style] for sec in covarianceSuite.covarianceSections if sec.columnData is None
                    and sec.rowData.ENDF_MFMT == pointer.ENDF_MFMT]
            if len(pointed_to) != 1:
                thisMFMT = summedMatrix.findAttributeInAncestry("rowData").ENDF_MFMT
                warningList.append( "Covariance for MF,MT=%s attempts to sum over non-existant covariance %s"
                        % (thisMFMT, pointer.ENDF_MFMT) )
                info.doRaise.append( warningList[-1] )
                continue
            pointer.link = pointed_to[0]
            pointer.root = None # internal link

    # fix lumped channel covariances (MT851-871) and summed channels (MT1,3,4,103-107)
    for lumpedKey in cov_info['lumpedChannels'].keys():
        if lumpedKey not in cov_info['MTL'] and lumpedKey not in cov_info['MTL_2']:
            warningList.append( "MT%d lumped sum covariance is present, but no reactions are included in the sum!" % lumpedKey[0] )
            info.doRaise.append( warningList[-1] )

    summedReactions = cov_info['MTL'].copy()
    summedReactions.update( cov_info['MTL_2'] )
    for (mt,mf) in sorted(summedReactions):
        try:
            lumpedChannels = cov_info['lumpedChannels'][(mt, mf)]
        except KeyError as e:
            warningList.append("Cannot find lumped channel %s" % str(e))
            continue

        Qs = []
        for (mt2, mf2) in summedReactions[(mt, mf)]:
            if mt not in range(851, 872) and mt2 not in cov_info['MTdict']: continue
            reacs = [r1 for r1 in info.reactionSuite.reactions if r1.ENDF_MT == mt2]
            if len(reacs) != 1: continue
            reac = reacs[0]
            xsc = reac.crossSection
            lumpedChannels.summands.append(sumsModule.Add(link=xsc))
            Qs.append(reac.thresholdQAs('eV'))
        if len(lumpedChannels.summands) == 0:
            warningList.append("MT%d: trying to sum over empty list (may be caused by earlier errors)" % mt)
            info.doRaise.append(warningList[-1])
        else:
            # mean value is not specified but can be computed (may require resonance reconstruction first)
            if lumpedChannels.ENDF_MT in (1, 3) and info.reactionSuite.supportsResonanceReconstruction():
                # need to generate resonancesWithBackground as 'evaluated' form.
                resWithBkg, other = [], []
                for summand in lumpedChannels.summands:
                    evaluated = summand.link.evaluated
                    if isinstance(evaluated, crossSectionModule.ResonancesWithBackground):
                        resWithBkg.append(evaluated)
                    else:
                        other.append(evaluated)

                def accumulateSum(term1, term2):
                    term2 = term2.toPointwise_withLinearXYs(lowerEps=1e-8)
                    term1, term2 = term1.mutualify(0, 1e-8, 0, term2, 0, 1e-8, 0)
                    return term1 + term2

                backgrounds = [None, None, None]
                for rwb in resWithBkg:
                    for idx, region in enumerate(('resolvedRegion', 'unresolvedRegion', 'fastRegion')):
                        regionData = getattr(rwb.background, region)
                        if regionData is None: continue
                        if backgrounds[idx] is None:
                            backgrounds[idx] = regionData.data.toPointwise_withLinearXYs(lowerEps=1e-8)
                        else:
                            backgrounds[idx] = accumulateSum(backgrounds[idx], regionData.data)

                sum_ = other[0].toPointwise_withLinearXYs(lowerEps=1e-8)
                for idx in range(1, len(other)):
                    sum_ = accumulateSum(sum_, other[idx])
                sum_ = sum_.trim()  # remove extra zeros

                if sum_.domainMin >= 0.999 * backgrounds[2].domainMin:  # allow for round-off error
                    backgrounds[2] = accumulateSum( backgrounds[2], sum_.domainSlice(backgrounds[2].domainMin, None) )
                else:
                    # some background needs to be added to RRR / URR background:
                    for idx, term in enumerate(backgrounds):
                        if term is None: continue
                        slice_ = sum_.domainSlice(term.domainMin, term.domainMax)
                        backgrounds[idx] = accumulateSum( slice_, term )

                resonanceLink = crossSectionModule.ResonanceLink( link = resonances )
                wrappedBackgrounds = []
                for term, wrapper in zip(backgrounds, (crossSectionModule.ResolvedRegion,
                                        crossSectionModule.UnresolvedRegion, crossSectionModule.FastRegion)):
                    if term is None:
                        wrappedBackgrounds.append(None)
                    else:
                        wrappedBackgrounds.append( wrapper(term) )

                background_ = crossSectionModule.Background( *wrappedBackgrounds )
                newXsc = crossSectionModule.ResonancesWithBackground(
                    info.style, resonanceLink, background_) # , backgroundForm.uncertainty )

            else:
                newXsc = lumpedChannels.sumSummands()
                newXsc.label = info.style

            covSec, = [sec for sec in covarianceSuite.covarianceSections if sec.columnData is None and
                       sec.rowData.ENDF_MFMT == '33,%d' % mt]
            newXsc.uncertainty = uncertaintiesModule.Uncertainty( functional=uncertaintiesModule.Covariance(
                link = covSec[info.style], root='$covariance' ) )
            lumpedChannels.crossSection.add( newXsc )
            lumpedChannels.Q.add(QModule.Constant1d(max(Qs), newXsc.domainMin, newXsc.domainMax, axes=QModule.defaultAxes('eV'), label=info.style))
        info.reactionSuite.sums.crossSectionSums.add( lumpedChannels )

    for exReac in sorted(set(cov_info['externalFiles'])):
        covarianceSuite.externalFiles.add( externalFileModule.ExternalFile(exReac,'FIXME/need/path/to/file') )

    if cov_info.get('MF34_missingFrames'):
        for (mt, LegendreLVals) in cov_info['MF34_missingFrames'].items():
            matchingreactions=[tmp for tmp in info.reactionSuite.reactions if tmp.ENDF_MT==mt]
            if matchingreactions:
                reaction = matchingreactions[0]
                # MF=34 is only for neutrons, so only need to look at neutron product:
                frame = reaction.outputChannel.getProductWithName(IDsPoPsModule.neutron).distribution[ info.style ].productFrame
                for lval in LegendreLVals:
                    evaluated = lval[info.style]
                    if isinstance(evaluated, covarianceMatrixModule.CovarianceMatrix):
                        evaluated.productFrame = frame
                    else:
                        for subsec in evaluated:
                            subsec.productFrame = frame

    sys.stdout.flush( )
    if( verbose > 0 ) :
        for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )
    return covarianceSuite, linkData
