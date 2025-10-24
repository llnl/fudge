# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
reactionSuite.py contains the 'reactionSuite' class that in turn holds all reactions for a given target/projectile.
reactionSuite is the top-level class for the GNDS structure.
"""

import re
import sys
import os
import copy

import pathlib
import numpy

from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule
from LUPY.hdf5 import HDF5_present, h5py
from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from LUPY import checksums as checksumsModule

from PoPs import IDs as IDsPoPsModule
from PoPs import database as PoPsModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.chemicalElements import chemicalElement as chemicalElementPoPsModule
from PoPs.decays import misc as decaysMiscPoPsModule
from PoPs.families import unorthodox as unorthodoxPoPsModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import vector as vectorModule
from xData import matrix as matrixModule
from xData import productArray as productArrayModule
from xData import XYs1d as XYs1dModule

from . import enums as enumsModule
from fudge.core.utilities import guessInteraction as guessInteractionModule
from fudge.processing import group as groupModule
from fudge.processing import nuclearPlusCoulombInterference as nuclearPlusCoulombInterferenceModule
from fudge.processing.deterministic import multiGroupUpScatter as multiGroupUpScatterModule
from fudge.processing.deterministic import particles as particlesModule
from fudge.processing import transporting as transportingModule

from . import suites as suitesModule
from . import sums as sumsModule
from . import styles as stylesModule
from . import warning as warningModule
from . import product as productModule
from . import institution as institutionModule
from .resonances import resonances as resonancesModule

from .reactionData import crossSection as crossSectionModule
from .reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule
from .outputChannelData import Q as QModule
from .productData import averageProductEnergy as averageProductEnergyModule
from .productData.distributions import angular as angularModule

from brownies.legacy.toENDF6 import ENDFconversionFlags as ENDFconversionFlagsModule


class NullDevice:
    """For internal use. Used in methods to set logFile when one is not entered."""

    def write(self, **kwargs): pass


class ReactionSuite(ancestryModule.AncestryIO):
    """
    This is the main class for a gnds projectile/target object. It contains
        * optional external files (e.g. for covariances or binary data stored in another file)
        * a 'styles' section describing all evaluated and processed data styles, documentation, etc.
        * a 'PoPs' Database with a list of all particles encountered in the file
        * optional resonance data
        * a list of reactions
        * additional lists: summed reactions (e.g. total), production reactions, 'orphan products', etc.
    """

    moniker = 'reactionSuite'
    ancestryMembers = ('externalFiles', 'styles', 'PoPs', 'resonances', 'reactions', 'orphanProducts', 'productions',
                       'incompleteReactions', 'fissionComponents', 'sums', 'applicationData')
    childNodeOrder = {
            GNDS_formatVersionModule.version_1_10: (suitesModule.ExternalFiles.moniker,                stylesModule.Styles.moniker, 
                                                 'documentations',                                  PoPsModule.Database.moniker,
                                                 resonancesModule.Resonances.moniker,               suitesModule.Reactions.moniker, 
                                                 suitesModule.OrphanProducts.moniker,               sumsModule.Sums.moniker, 
                                                 suitesModule.FissionComponents.moniker,            suitesModule.Productions.moniker,
                                                 suitesModule.IncompleteReactions.moniker,          suitesModule.ApplicationData.moniker ),
            GNDS_formatVersionModule.version_2_0_LLNL_4: (suitesModule.ExternalFiles.moniker,                 stylesModule.Styles.moniker,
                                                        PoPsModule.Database.moniker,
                                                        resonancesModule.Resonances.moniker,        suitesModule.Reactions.moniker, 
                                                        suitesModule.OrphanProducts.moniker,                sumsModule.Sums.moniker, 
                                                        suitesModule.FissionComponents.moniker,             suitesModule.Productions.moniker,
                                                        suitesModule.IncompleteReactions.moniker,           suitesModule.ApplicationData.moniker),
            GNDS_formatVersionModule.version_2_0:      (suitesModule.ExternalFiles.moniker,                 stylesModule.Styles.moniker,
                                                        PoPsModule.Database.moniker,
                                                        resonancesModule.Resonances.moniker,        suitesModule.Reactions.moniker, 
                                                        suitesModule.OrphanProducts.moniker,                sumsModule.Sums.moniker, 
                                                        suitesModule.FissionComponents.moniker,             suitesModule.Productions.moniker,
                                                        suitesModule.IncompleteReactions.moniker,           suitesModule.ApplicationData.moniker)}

    def __init__(self, projectile, target, evaluation, interaction=None, formatVersion=GNDS_formatVersionModule.default,
                 style=None, projectileFrame=xDataEnumsModule.Frame.lab, MAT=None, PoPs=None, sourcePath=None):
        """
        Creates a new reactionSuite object, containing all reaction data for a projectile/target combination.

        :param projectile:      Projectile ID (string)
        :param target:          Target ID (string)
        :param evaluation:      The evaluation for this protare (string)
        :param interaction:     The GNDS interaction for the protare (string)
        :param formatVersion:   GNDS format version (string)
        :param style:           Indicates style of data in this reactionSuite (evaluated, processed, etc.).  Instance of styles.style
        :param projectileFrame: Frame of the projectile. Allowed values are 'lab' or 'centerOfMass'
        :param MAT:             Integer ENDF MAT number (only needed for writing back to ENDF-6)
        :param PoPs:            A PoPs Database that is copied to self's PoPs Database
        :return:
        """

        ancestryModule.AncestryIO.__init__(self)

        if not isinstance(projectile, str): raise TypeError('projectile ID not a string')
        self.__projectile = projectile

        if not isinstance(target, str): raise TypeError('target ID not a string')
        self.__target = target

        if not isinstance(evaluation, str): raise TypeError('evaluation not a string')
        self.__evaluation = evaluation

        if interaction == enumsModule.Interaction.legacyTNSL:
            interaction = enumsModule.Interaction.TNSL
        self.interaction = interaction
        self.__guessedInteraction = False

        if formatVersion not in GNDS_formatVersionModule.allowedPlus:
            raise Exception("Unsupported GNDS PoPs structure '%s'!" % str(formatVersion))
        self.formatVersion = formatVersion

        self.__projectileFrame = xDataEnumsModule.Frame.checkEnumOrString(projectileFrame)

        self.MAT = MAT

        self.__sourcePath = sourcePath

        self.__externalFiles = suitesModule.ExternalFiles()
        self.__externalFiles.setAncestor(self)

        self.__styles = stylesModule.Styles()
        self.__styles.setAncestor(self)
        if style is not None: self.styles.add(style)

        if PoPs is not None:
            self.__PoPs = PoPs.copy()
        else :
            self.__PoPs = PoPsModule.Database('protare_internal', '1.0')
        self.__PoPs.setAncestor(self)

        self.__resonances = None

        self.__reactions = suitesModule.Reactions()
        self.__reactions.setAncestor(self)

        self.__orphanProducts = suitesModule.OrphanProducts()
        self.__orphanProducts.setAncestor(self)

        self.__sums = sumsModule.Sums()
        self.__sums.setAncestor(self)

        self.__productions = suitesModule.Productions()
        self.__productions.setAncestor(self)

        self.__incompleteReactions = suitesModule.IncompleteReactions()
        self.__incompleteReactions.setAncestor(self)

        self.__fissionComponents = suitesModule.FissionComponents()
        self.__fissionComponents.setAncestor(self)

        self.__applicationData = suitesModule.ApplicationData()
        self.__applicationData.setAncestor(self)

        self._externalLinks = []        # Keep track of links to external files.
        self._loadedCovariances = None   # If loadCovariances called, store results here.

    def __iter__(self):

        for reaction in self.reactions: yield reaction
        for orphanProduct in self.orphanProducts: yield orphanProduct
        for _sum in self.sums.crossSectionSums: yield _sum
        for reaction in self.fissionComponents: yield reaction
        for reaction in self.productions: yield reaction
        for reaction in self.incompleteReactions: yield reaction

    def __str__(self):
        """See method toString."""

        string = self.inputParticlesToReactionString() + (' (%s)' % self.evaluation)
        
        return string

    @property
    def compound_nucleus(self):
        """
        Compute the compound nucleus id corresponding to this reaction.

        Note, we don't check to see if the compound nucleus is a valid concept, so
        this will fail on electro-atomic, photo-atomic and natural element target data
        """

        ZAproj = chemicalElementMiscPoPsModule.ZA(self.PoPs[self.projectile])
        ZAtarg = chemicalElementMiscPoPsModule.ZA(self.PoPs[self.target])
        Zcompound, Acompound = divmod(ZAproj + ZAtarg, 1000)
        return chemicalElementMiscPoPsModule.idFromZAndA(Zcompound, Acompound)

    @property
    def format(self):

        return self.formatVersion

    @property
    def projectile(self):
        """Returns self's projectile."""

        return self.__projectile

    @property
    def target(self):
        """Returns self's target."""

        return self.__target

    @property
    def projectileFrame(self):

        return self.__projectileFrame

    @property
    def evaluation(self):
        """Returns self's evaluation."""

        return self.__evaluation

    @evaluation.setter
    def evaluation(self, value):

        if not isinstance(value, str):
            raise TypeError("Value must be a string")

        self.__evaluation = value

    @property
    def interaction(self):
        """Returns self's interaction."""

        return self.__interaction

    @interaction.setter
    def interaction(self, value):

        if value is None:
            value = enumsModule.Interaction.nuclear
            print(f'Need to specify interaction when calling reactionSuite.__init__. Setting it to "{value}".'
                  'Please update your code as in the future this will execute a raise.' )

        self.__interaction = enumsModule.Interaction.checkEnumOrString(value)

    @property
    def guessedInteraction(self):
        """
        Returns the guessedInteraction value. This is for pre GNDS 2.0 use only and should not be used with GNDS 2.0 or higher
        (i.e., it should always be False for GNDS 2.0 or higher.
        """

        return self.__guessedInteraction

    @guessedInteraction.setter
    def guessedInteraction(self, value):
        """
        Set self's guessedInteraction member. This is for pre GNDS 2.0 use only and should not be used with GNDS 2.0 or higher.
        """

        self.__guessedInteraction = value

    @property
    def sourcePath(self):
        """Returns the sourcePath member which is the path to the reactionSuite file for self if self is from a file."""

        return self.__sourcePath

    @property
    def externalFiles(self):
        """Returns the externalFiles instance."""

        return self.__externalFiles

    @property
    def styles(self):
        """Returns the PoPs instance."""

        return self.__styles
    
    @property
    def PoPs(self):
        """Returns the PoPs instance."""

        return self.__PoPs

    @property
    def resonances(self):
        return self.__resonances

    @resonances.setter
    def resonances(self, resonances):
        """Set resonances to the specified resonances instance."""

        if not isinstance(resonances, resonancesModule.Resonances):
            raise TypeError("Must be Resonances instance, not '%s'" % type(resonances))
        self.__resonances = resonances
        self.__resonances.setAncestor(self)

    @property
    def reactions(self):
        """Returns the reactions instance."""

        return self.__reactions

    @property
    def orphanProducts(self):
        """Returns the orphanProducts instance."""

        return self.__orphanProducts

    @property
    def sums(self):
        """Returns the sums instance."""

        return self.__sums

    @property
    def productions(self):
        """Returns the productions instance."""

        return self.__productions

    @property
    def incompleteReactions(self):
        """Returns the incompleteReactions instance."""

        return self.__incompleteReactions

    @property
    def fissionComponents(self):
        """Returns the fissionComponents instance."""

        return self.__fissionComponents

    @property
    def applicationData(self):
        """Returns the applicationData instance."""

        return self.__applicationData

    @property
    def domainMin(self):

        return self.reactions[0].crossSection.domainMin

    @property
    def domainMax(self):

        return self.reactions[0].crossSection.domainMax

    @property
    def domainUnit(self):

        return self.reactions[0].crossSection.domainUnit

    def isThermalNeutronScatteringLaw(self):

        for reaction in self.reactions:
            if reaction.isThermalNeutronScatteringLaw(): return True
        return False

    def photonBranchingData(self):
        """
        For each isotope in *self.PoPs*, calls decaysMiscPoPsModule.photonBranchingData.

        :return:        Dictionary with keys being each isotope in *self.PoPs* and their values being that returned by
                            decaysMiscPoPsModule.photonBranchingData.
        """

        branchingData = {}
        for chemicalElement in self.PoPs.chemicalElements:
            for isotope in chemicalElement:
                if isotope.symbol in self.PoPs:
                    branchingData[isotope.symbol] = decaysMiscPoPsModule.photonBranchingData(self.PoPs, isotope.symbol)
        return branchingData

    def convertUnits(self, unitMap):
        """
        Recursively searches for units within all sections, converting units that appear in the unitMap.
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        self.styles.convertUnits(unitMap)
        if self.resonances is not None: self.resonances.convertUnits(unitMap)
        self.reactions.convertUnits(unitMap)
        self.orphanProducts.convertUnits(unitMap)
        self.sums.convertUnits(unitMap)
        self.productions.convertUnits(unitMap)
        self.incompleteReactions.convertUnits(unitMap)
        self.PoPs.convertUnits(unitMap)
        self.fissionComponents.convertUnits(unitMap)
        self.applicationData.convertUnits(unitMap)

    def fixDomains(self, energyMax):
        """
        For each reaction type and summed data, makes all projectile energy domains be consistent. This only modifies *evaluated* data and
        not any *processed* data.
        """

        numberOfFixes, labels = self.styles.fixDomains(energyMax)
        numberOfFixes += self.reactions.fixDomains(labels, 0.0, energyMax)               # The energyMin value is ignored by all "reaction" exclusive suites.
        numberOfFixes += self.orphanProducts.fixDomains(labels, 0.0, energyMax)
        numberOfFixes += self.productions.fixDomains(labels, 0.0, energyMax)
        numberOfFixes += self.fissionComponents.fixDomains(labels, 0.0, energyMax)
        numberOfFixes += self.incompleteReactions.fixDomains(labels, 0.0, energyMax)
        numberOfFixes += self.sums.fixDomains(labels, energyMax)

        return numberOfFixes

    def check(self, **kwargs):
        """
        Check all data in the reactionSuite, returning a gnds.warning.context object with list of warnings.

        Currently supported options::

            'branchingRatioSumTolerance' 1e-6        # branching ratios must sum to 1 (within this tolerance)
            'dQ'                         '1e-3 MeV'  # tolerance for comparing tabulated / calculated Q values
            'dThreshold'                 '1e-3 MeV'
            'crossSectionEnergyMax'      '20 MeV'    # warn if cross section ends before this energy
            'crossSectionOnly':          False
            'crossSectionMaxDiff':       1e-3        # for comparing crossSectionSum to summands
            'multiplicityMaxDiff':       1e-3        # for comparing multiplicitySum to summands
            'transportables'             ('n',)      # distribution required for these products
            'normTolerance':             1e-5        # for checking distribution normalization
            'checkEnergyBalance'         True
            'reconstructResonances'      True
            'dEnergyBalanceRelative'     1e-3
            'dEnergyBalanceAbsolute'     '1 eV'
            'fissionEnergyBalanceLimit'  0.15        # at least 85% of available energy should go to fission products
            'failOnException'            False       # if True, crash instead of converting exceptions to warnings
            'cleanUpTempFiles'           True        # remove derived data that was produced during checking
            'verbose'                    False       # be verbose while running checks

        Currently unused options::

            'checkForEnDepData'         False
            'allowZeroE'                False
            'xCloseEps'                 None
            'maxAbsFloatValue'          None
            'checkEMin'                 True
            'checkMissing'              True
            'maxGammaMultiplicity'      100.
          """

        from fudge import warning

        options = {
                'branchingRatioSumTolerance': 1e-6,
                'dQ': '1e-3 MeV',
                'dThreshold': '1e-3 MeV',
                'crossSectionEnergyMin': f'{self.domainMin} {self.domainUnit}',
                'crossSectionEnergyMax': f'{self.domainMax} {self.domainUnit}',
                'crossSectionOnly': False,
                'crossSectionMaxDiff': 1e-3,
                'multiplicityMaxDiff': 1e-3,
                'transportables': ('n', 'photon', 'H1', 'H2', 'H3', 'He3', 'He4'),
                'normTolerance': 1e-5,
                'checkEnergyBalance': True,
                'reconstructResonances': True,
                'dEnergyBalanceRelative': 1e-3,
                'dEnergyBalanceAbsolute': '1 eV',
                'fissionEnergyBalanceLimit': 0.15,
                'failOnException': False,
                'cleanUpTempFiles': True,
                'verbose': False,
                # currently unused options:
                'checkForEnDepData': False,
                'allowZeroE': False,
                'xCloseEps': None,
                'maxAbsFloatValue': None,
                'checkEMin': True,
                'checkMissing': True,
                'maxGammaMultiplicity': 100.,
                }
        for key in kwargs:
            if key in options:
                options[key] = kwargs[key]
            else:
                raise KeyError("check() received unknown keyword argument '%s'" % key)

        warnings = []

        if isinstance(self.target, PoPsModule.unorthodoxModule.Particle):
            warnings.append(warning.UnorthodoxParticleNotImplemented(self.target))
            return warning.Context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)

        evaluatedStyle = self.styles.getEvaluatedStyle()
        if evaluatedStyle is None:
            warnings.append(warning.NotImplemented("Checking currently only supported for 'evaluated' style"))
            return warning.Context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)

        # assemble some useful info, to be handed down to children's check() functions:
        incidentEnergyUnit = self.domainUnit
        # FIXME warn if incidentEnergyUnit not one of 'eV', 'MeV'?
        massUnit = f'{incidentEnergyUnit}/c**2'
        projectile = self.PoPs[self.projectile]
        projectileZ, projectileA, projectileZA, projectileLevelIndex = chemicalElementMiscPoPsModule.ZAInfo(projectile)
        target = self.PoPs.final(self.target)
        targetZ, targetA, targetZA, targetLevelIndex = chemicalElementMiscPoPsModule.ZAInfo(target)
        if targetA == 0:
            compoundZA = 1000 * (projectileZ + targetZ)     # elemental target
        else:
            compoundZA = targetZA + projectileZA
        CoulombChannel = (targetZ != 0) and (projectileZ != 0)
        elementalTarget = targetA == 0
        if not elementalTarget:
            projectileMass_amu = projectile.getMass('amu')
            targetMass_amu = target.getMass('amu')
            kinematicFactor = (targetMass_amu + projectileMass_amu) / targetMass_amu    # (M+m)/M

            projectileMass = projectile.getMass(massUnit)
            targetMass = target.getMass(massUnit)
            availableEnergy = projectileMass + targetMass
        else:
            # For elemental targets, calculating these factors doesn't make sense since there is no defined target mass
            if self.interaction is enumsModule.Interaction.nuclear:
                warnings.append(warning.ElementalTarget(self.target))
            kinematicFactor = 1.0
            availableEnergy = None

        info = {'reactionSuite': self, 'kinematicFactor': kinematicFactor, 'compoundZA': compoundZA,
                'availableEnergy': availableEnergy, 'CoulombChannel': CoulombChannel, 'style': evaluatedStyle,
                'reconstructedStyleName': None, 'elementalTarget': elementalTarget, 'interaction': self.interaction}
        info.update(options)
        if elementalTarget:
            # For elemental targets, calculating energy balance isn't possible
            info['checkEnergyBalance'] = False

        externalFileWarnings = []
        for externalFile in self.externalFiles:
            fullPath = os.path.join(
                os.path.split(self.sourcePath)[0], externalFile.path
            )
            if not os.path.exists(fullPath):
                warnings.append(warning.MissingExternalFile(externalFile.path))
                continue

            if externalFile.checksum is not None:
                algorithm = externalFile.algorithm
                checksum = checksumsModule.checkers[algorithm].from_file(fullPath)
                if checksum != externalFile.checksum:
                    warnings.append(
                        warning.WrongExternalFileChecksum(externalFile.path,
                                                          externalFile.checksum,
                                                          checksum))

        if externalFileWarnings:
            warnings.append(warning.Context('externalFiles', externalFileWarnings))

        # does evaluation start at an appropriate incident energy?
        # FIXME these limits may differ by library. Move to config file?
        expectedDomainMin = {
            IDsPoPsModule.neutron: '1e-5 eV',
            IDsPoPsModule.photon: '1 MeV',
            'H1': '1 MeV',
            'H2': '1 MeV',
            'H3': '1 MeV',
            'He3': '1 MeV',
            'He4': '1 MeV',
        }.get(self.projectile, '1 MeV')
        if self.interaction is enumsModule.Interaction.atomic:
            expectedDomainMin = {
                IDsPoPsModule.photon: '1 eV',
                IDsPoPsModule.electron: '10 eV'
            }.get(self.projectile, '10 eV')
        if PQUModule.PQU(self.domainMin, self.domainUnit) > PQUModule.PQU(expectedDomainMin):
            warnings.append(warning.EvaluationDomainMinTooHigh(expectedDomainMin, self))

        if self.resonances is not None:
            resonanceWarnings = self.resonances.check(info)
            if resonanceWarnings:
                warnings.append(warning.Context('resonances', resonanceWarnings))
            if options['reconstructResonances'] and self.resonances.reconstructCrossSection:
                info['reconstructedStyleName'] = self.styles.getTempStyleNameOfClass(stylesModule.CrossSectionReconstructed)
                # convert resonance parameters to pointwise data, interpolable to .1 percent:
                reconstructedStyle = stylesModule.CrossSectionReconstructed(
                    info['reconstructedStyleName'], derivedFrom=evaluatedStyle.label)
                try:
                    self.reconstructResonances(reconstructedStyle, 0.001, thin=False, verbose=False)
                except Exception as e:
                    warnings.append(warning.ExceptionRaised("when reconstructing resonances: %s" % e))
                    if info['failOnException']: raise

        if info['checkEnergyBalance']:
            # setup options for calculating average product energy and momentum
            averageProductDataStyle = stylesModule.AverageProductData(
                    label=self.styles.getTempStyleNameOfClass(stylesModule.AverageProductData),
                    derivedFrom=evaluatedStyle.label)

            self.styles.add(averageProductDataStyle)
            info['averageProductDataStyle'] = averageProductDataStyle
            info['averageProductDataArgs'] = {
                'verbosity': 1 if options['verbose'] else 0,   # additional required arguments
                'incrementalIndent': '  ', 'energyAccuracy': 1e-6, 'momentumAccuracy': 1e-6,
                'incidentEnergyUnit': incidentEnergyUnit, 'massUnit': massUnit,
                'momentumUnit': incidentEnergyUnit + '/c', 'projectileMass': projectileMass, 'targetMass': targetMass,
                'reactionSuite': self, 'photonBranchingData': self.photonBranchingData(),
                'isInfiniteTargetMass': isinstance(target, chemicalElementPoPsModule.ChemicalElement)}

        if self.projectile == 'n' and self.target != 'n':
            # test Wick's limit: 0-degree elastic xsc >= ( total xsc * k/4pi )^2
            elastic = self.getReaction('elastic')
            total = self.getReaction('total')
            if total is None:
                warnings.append(warning.TestSkipped("Wick's limit", "can't find reaction 'total'"))
            else:
                try:
                    if info['reconstructedStyleName'] in elastic.crossSection:
                        elastic_xsc = elastic.crossSection[info['reconstructedStyleName']]
                        total_xsc = total.crossSection[info['reconstructedStyleName']]
                    else:
                        elastic_xsc = elastic.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
                        total_xsc = total.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
                    elastic_distribution = elastic.outputChannel.getProductWithName('n').distribution[evaluatedStyle.label].angularSubform
                    if isinstance(elastic_distribution, angularModule.XYs2d):
                        linearized = elastic_distribution.toPointwise_withLinearXYs()
                        forward_scattering = crossSectionModule.XYs1d(data=[], axes=crossSectionModule.XYs1d.defaultAxes())
                        for energy_in in linearized.getEnergyArray():
                            forward_scattering.setValue(energy_in, linearized.interpolateAtValue(energy_in).evaluate(1.0))
                    elif isinstance(elastic_distribution, angularModule.Regions2d):
                        forward_scattering = crossSectionModule.Regions1d(axes=crossSectionModule.XYs1d.defaultAxes())
                        for region in elastic_distribution:
                            ptw = crossSectionModule.XYs1d(data=[])
                            linearized = region.toPointwise_withLinearXYs()
                            for energy_in in linearized.getEnergyArray():
                                ptw.setValue(energy_in, linearized.interpolateAtValue(energy_in).evaluate(1.0))
                            forward_scattering.append(ptw)
                        forward_scattering = forward_scattering.toPointwise_withLinearXYs(
                            accuracy=1e-5, lowerEps=1e-8, upperEps=1e-8)
                    else:
                        raise NotImplementedError(f"Checking Wick's limit for {type(elastic_distribution)} distribution")

                    mutualDomain = list(zip(
                        *[dat.domain() for dat in (elastic_xsc, total_xsc, forward_scattering)]))
                    mutualDomain = (max(mutualDomain[0]), min(mutualDomain[1]))
                    egrid = [e for e in forward_scattering.domainGrid if mutualDomain[0] <= e <= mutualDomain[1]]

                    # get probability at mu=1.0:
                    forward_scattering = [forward_scattering.evaluate(e) for e in egrid]
                    elastic_xsc = [elastic_xsc.evaluate(e) for e in egrid]
                    total_xsc = [total_xsc.evaluate(e) for e in egrid]

                    wlcons = 3.05607e-8 * kinematicFactor**2   # ( sqrt(2 * neutronMass) / (4 * pi * hbar) * (M+m)/M )^2 in 1/(eV*b)
                    for i1 in range(len(egrid)):
                        if forward_scattering[i1] * elastic_xsc[i1] < wlcons * egrid[i1] * total_xsc[i1]**2:
                            ratio = (forward_scattering[i1] * elastic_xsc[i1]) / (wlcons * egrid[i1] * total_xsc[i1]**2)
                            warnings.append(warning.WicksLimitError( 1-ratio, egrid[i1]))
                except Exception as e:
                    warnings.append(warning.ExceptionRaised("when checking Wick's limit: %s" % e))
                    if info['failOnException']: raise

        particleWarnings = self.PoPs.check(**info)
        if particleWarnings: warnings.append(warning.Context('PoPs', particleWarnings))

        info['isInfiniteTargetMass'] = isinstance(target, chemicalElementPoPsModule.ChemicalElement)
        for reaction in self:
            reactionWarnings = reaction.check(info)
            if reactionWarnings:
                warnings.append(
                    warning.Context('%s label %s (MT%d)' % (reaction.moniker, reaction.label, reaction.ENDF_MT), reactionWarnings))

        result = warning.Context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)
        result.info = info

        if options['cleanUpTempFiles']:
            if info['reconstructedStyleName'] is not None:
                self.removeStyle(info['reconstructedStyleName'])

        return result

    def diff(self, other):

        diffResults1 = warningModule.DiffResults(self, other)

        if not isinstance(other, ReactionSuite): raise ValueError('Other not a reactionSuite instance.')
        if self.domainUnit != other.domainUnit: raise ValueError('Energy units not the same.')

        for style in self.styles:
            if not isinstance(style, (stylesModule.Evaluated, stylesModule.CrossSectionReconstructed)):
                raise TypeError('Style "%s" not allowed when diffing.' % style.moniker)

        for style in other.styles:
            if not isinstance(style, (stylesModule.Evaluated, stylesModule.CrossSectionReconstructed)):
                raise TypeError('Style "%s" not allowed when diffing.' % style.moniker)

        self.reactions.diff(other.reactions, diffResults1)

        return diffResults1

    def elasticReactionLabel(self):
        """
        Returns the label for the elastic channel. When constructing the label, if target is an alias, use the alias' pid; otherwise, use the target id.
        """

        return '%s + %s' % (self.projectile, self.targetID(protareID=False))

    def targetID(self, protareID=False):
        """
        Returns the ID for the target. If the target is a nuclear meta-stable, then the target ID as listed by the 'target' attribute for the
        reactionSuite (i.e., protare) will be the meta-stable's ID while everywhere else it will be the nuclide excitation ID. Otherwise,
        the target ID as listed by the 'target' attribute for the reactionSuite is the same as it is everywhere else. For example, for the
        protare 'n + Am242', the target will be listed as 'Am242' everywhere, while for the protare 'n + Am242_m1', the target will be listed
        as 'Am242_m1' by the 'target' attribute for the reactionSuite and 'Am242_e2' everywhere else. Note, for the protare 'n + Am242_m1'
        the PoPs.aliases instance must have the alias 'Am242_m1' for 'Am242_e2'. If protareID is True, the target ID as listed by the 'target'
        attribute for the reactionSuite will be returned, otherwise, the target ID list everywhere else will be returned.
        """

        if protareID: return self.target
        return self.PoPs.final(self.target).id

    def findEntity(self, entityName, attribute=None, value=None):
        """
        Overrides ancestry.findEntity. reactionSuite contains several different lists,
        so may need to descend into those to find desired entity.
        """

        if entityName in ('reaction', 'summedReaction', 'fissionComponent', 'production', 'reactionSum'):
            for entity in getattr(self, entityName+'s'):
                if getattr(entity, attribute, None) == value:
                    return entity
        return ancestryModule.Ancestry.findEntity(self, entityName, attribute, value)

    def findLinks(self):

        links = []
        for ancestryMember in self.ancestryMembers:
            if hasattr(self, ancestryMember):
                item = getattr(self, ancestryMember)
                if hasattr(item, 'findLinks'): getattr(item, 'findLinks')(links)
        return links

    def fixLinks(self):

        ENDFconversions = []
        if 'LLNL' in self.applicationData:
            ENDFconversions = [LLNL_data for LLNL_data in self.applicationData['LLNL']
                               if isinstance(LLNL_data, ENDFconversionFlagsModule.ENDFconversionFlags)]
            if len(ENDFconversions) > 1:
                raise Exception('Too many ENDFconversionFlags in applicationData["LLNL"] section.')
            if len(ENDFconversions) == 1:
                ENDFconversions = ENDFconversions[0]

        links = self.findLinks()
        for instance, link, xlink in links:
            link2 = link
            if xlink[0] == '/':          # Cannot look for relative xlink as it is relative to link which is not known.
                try:
                    link2 = self.followXPath(instance.path)
                except ancestryModule.XPathNotFound:
                    link2 = None
                except:
                    raise

            if link2 is None:
                if isinstance(instance, ENDFconversionFlagsModule.Conversion):
                    ENDFconversions.flags.pop(ENDFconversions.flags.index(instance))
                elif isinstance(instance, sumsModule.Add):
                    summands = instance.ancestor
                    summands.pop(summands.index(instance))
            elif link is not link2:
                instance.link = link2

    def hasParticle(self, name):

        return name in self.PoPs

    def heatCrossSections(self, style, lowerlimit=None, upperlimit=None, interpolationAccuracy=0.001,
                          heatAllPoints=False, doNotThin=True, heatBelowThreshold=True, heatAllEDomain=True,
                          setThresholdToZero=False, verbose=0):

        EMin = self.domainMin
        kwargs = {
            'lowerlimit': lowerlimit, 'upperlimit': upperlimit, 'interpolationAccuracy': interpolationAccuracy,
            'heatAllPoints':  heatAllPoints, 'doNotThin': doNotThin, 'heatBelowThreshold': heatBelowThreshold,
            'heatAllEDomain': heatAllEDomain, 'setThresholdToZero': setThresholdToZero, 'verbose': verbose}

        for reaction in self.reactions:
            reaction.heatCrossSection(style, EMin, **kwargs)

        for orphanProduct in self.orphanProducts:
            orphanProduct.heatCrossSection(style, EMin, **kwargs)

        for production in self.productions:
            production.heatCrossSection(style, EMin, **kwargs)

        for _sum in self.sums.crossSectionSums:
            _sum.heatCrossSection(style, EMin, **kwargs)

        for fissionComponent in self.fissionComponents:
            fissionComponent.heatCrossSection(style, EMin, **kwargs)

    def getParticle(self, name):

        return self.PoPs[name]

    def getMassRatio(self):

        if not hasattr(self, '__massRatio'):
            M = self.PoPs[self.target].getMass('amu')
            m = self.PoPs[self.projectile].getMass('amu')
            self.__massRatio = (M / (M+m))
        return self.__massRatio

# BRB6 This does not belong here.
    def getIsotopeName(self, *args):
        """
        Return name of compound nucleus formed by nuc1+nuc2+...
        if a nucleus in 'args' starts with '-', subtract instead.
        """

# CALEB: This logic need to use PoPs stuff.
        def parse(name):
            # FIXME remove '_natural'? Now using 'Os0' instead to indicate elemental target
            if IDsPoPsModule.photon in name: return 1, 0, 0, 1
            sign, mult, symbol, A = re.search(r"([-]?)([0-9]+)?([a-zA-Z]+)(_natural|[0-9]*)", name).groups()
            if not mult: mult = 1
            if (symbol == IDsPoPsModule.neutron) and not A: A = 1
            if A != '_natural': A = int(A)
            if symbol == IDsPoPsModule.neutron:
                Z = 0
            else:
                Z = chemicalElementMiscPoPsModule.ZFromSymbol[symbol]
            return int(sign + '1'), Z, A, int(mult)

        retZ, retA = 0, None
        try:
            for nucleus in args:
                if 'photon' in nucleus: continue
                sign, Z, A, mult = parse(nucleus)
                retZ += sign*mult*Z
                if retA is None:
                    retA = A
                elif A == 0 or retA == 0:
                    retA = 0
                else:
                    retA += sign*mult*A
            return '%s%s' % (chemicalElementMiscPoPsModule.symbolFromZ[retZ], retA)
        except:
            print("      WARNING: couldn't extract isotope name from product list!")

# BRB6 many hardwired in this method.
    def getReaction(self, channel):
        """
        Search list of reactions for a specified channel.
        The 'channel' argument should be either a reaction type ('elastic','capture','fission', etc) or
        the list of outgoing particles ('n + Pu239' for example).

        If 'channel' is an int, then we assume it's an ENDF MT

        Returns None if no matching reaction is found.
        """

        if channel == 'capture':
            channel = 102

        # If channel is an int, it might be an ENDF MT, so check those first
        if isinstance(channel, (int, numpy.integer)):   # must be an ENDF MT
            for reactionList in self.reactions, self.sums.crossSectionSums, self.fissionComponents, self.productions:
                for reaction in reactionList:
                    if reaction.ENDF_MT == channel: return reaction
            return None

        # translate special channel names:
        if channel == 'elastic':
            channel = channel_tr = self.elasticReactionLabel()

            # FIXME we're not dealing with metastables or other aliases very well right now, need to rework this.
            # special treatment to check for elastic reactions with labels like 'n + Am242_m1' (instead of Am242_e2):
            if self.target in self.PoPs.aliases:
                metastableElastic = f"{self.projectile} + {self.target}"
                if metastableElastic in self.reactions:
                    channel = channel_tr = metastableElastic

        elif channel == 'capture':
            channel_tr = 'z,gamma'
        else:
            channel_tr = channel

        # check if 'channel' == one of the fudge reactions:
        chStrings, reacs = [], []
        for reaction in self:
            if str(reaction) == channel: return reaction
            chStrings.append(str(reaction))
            reacs.append(reaction)

        # make list containing a set() of products for each reaction.
        chStrings = [tmp.replace('(', '').replace(')', '').replace('->', '+') for tmp in chStrings]
        chSets = [set(tmp.split(' + ')) for tmp in chStrings]
        chSetsNoPhoton = [tmp.copy() for tmp in chSets]
        for tmp in chSetsNoPhoton: tmp.discard('photon')
        if set(channel.split(' + ')) in chSets:
            return reacs[chSets.index(set(channel.split(' + ')))]
        if set(channel.split(' + ')) in chSetsNoPhoton:
            return reacs[chSetsNoPhoton.index(set(channel.split(' + ')))]

        if 'fission' in channel.lower():
            channel_fiss = channel.lower()

            def matchlist(*args):
                return any([a in channel_fiss for a in args])

            if channel_fiss == 'fission' or 'total' in channel_fiss:
                fissionGenre = enumsModule.FissionGenre.total
            elif matchlist('first', '1st'):
                fissionGenre = enumsModule.FissionGenre.firstChance
            elif matchlist('second', '2nd'):
                fissionGenre = enumsModule.FissionGenre.secondChance
            elif matchlist('third', '3rd'):
                fissionGenre = enumsModule.FissionGenre.thirdChance
            elif matchlist('fourth', '4th'):
                fissionGenre = enumsModule.FissionGenre.fourthChance
            else:
                print('Cannot determine fission genre from "%s".' % channel_fiss)
                fissionGenre = enumsModule.FissionGenre.none
            retVal  = [r1 for r1 in self.reactions         if r1.fissionGenre == fissionGenre]
            retVal += [r1 for r1 in self.fissionComponents if r1.fissionGenre == fissionGenre]
            if len(retVal) == 1:
                return retVal[0]
        else:
            # check if channel is in form '(z,2na)', 'n,4n' or similar:
            patt = re.match('^[(]?[znpdthag],([a-zA-Z0-9]+)[)]?$', channel_tr)
            if patt:
                thisChannelSet = set()
                match = re.findall('([1-9]?)(gamma|He3|[npagdt]?)[+]?', patt.groups()[0] )
                for mul, prod in match:
                    if not prod: continue
                    prod = {'g': IDsPoPsModule.photon, 'gamma': IDsPoPsModule.photon, 'n': 'n', 'p': 'H1',
                            'd': 'H2', 't': 'H3', 'h': 'He3', 'a': 'He4'}[prod]
                    if mul: prod = mul + prod
                    if prod in thisChannelSet:
                        raise KeyError("Please specify multiplicity explicitly ('z,2n' instead of 'z,nn')")
                    thisChannelSet.add(prod)
                if not thisChannelSet:
                    return None
                # also add residual to the set:
                proj, target = self.projectile, self.target
                thisChannelSet.add( self.getIsotopeName(*([proj, target] + ['-'+a for a in thisChannelSet])))
                if thisChannelSet in chSets:
                    return reacs[chSets.index(thisChannelSet)]
                if thisChannelSet in chSetsNoPhoton:
                    return reacs[chSetsNoPhoton.index(thisChannelSet)]

        return None

    def calculateAverageProductData( self, style, indent = '', additionalReactions = [], **kwargs ) :
        """
        Calculate average energy and momentum data for all products of all reactions.
        Resulting data are stored within each product. Example usage is:

        from fudge import reactionSuite as reactionSuiteModule
        from fudge import styles as stylesSuiteModule
        reactionSuite = reactionSuiteModule.ReactionSuite.readXML_file( "16.xml" )

        style = stylesModule.AverageProductData( label = 'productData', 'eval' )
        reactionSuite.calculateAverageProductData( style )

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        if( not( isinstance( style, stylesModule.AverageProductData ) ) ) : raise TypeError( 'Invalid style' )

        verbosity = kwargs.get( 'verbosity', 0 )
        kwargs['verbosity'] = verbosity

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        kwargs['incrementalIndent'] = incrementalIndent
        indent2 = indent + incrementalIndent

        logFile = kwargs.get( 'logFile', NullDevice( ) )
        kwargs['logFile'] = logFile

        energyAccuracy = kwargs.get( 'energyAccuracy', 1e-5 )
        kwargs['energyAccuracy'] = energyAccuracy
        momentumAccuracy = kwargs.get( 'momentumAccuracy', 1e-3 )
        kwargs['momentumAccuracy'] = momentumAccuracy

        kwargs['incidentEnergyUnit'] = self.domainUnit
        kwargs['momentumUnit'] = kwargs['incidentEnergyUnit'] + '/c'
        kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'

        kwargs['projectileMass'] = self.PoPs[self.projectile].getMass( kwargs['massUnit'] )
        target = self.PoPs[self.target]
        kwargs['isInfiniteTargetMass'] = isinstance( target, chemicalElementPoPsModule.ChemicalElement )
        if( kwargs['isInfiniteTargetMass'] ) :
            kwargs['targetMass'] = None
        else :
            if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.id]
            kwargs['targetMass'] = target.getMass( kwargs['massUnit'] )

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))

        kwargs['reactionSuite'] = self

        for reaction in self.reactions:
            try:
                reaction.calculateAverageProductData(style, indent=indent2, **kwargs)
            except CoulombPlusNuclearElasticModule.CoulombDepositionNotSupported:
                if not kwargs.get('skipCoulombPlusNuclearElasticError', False):
                    raise
                print('INFO: Skipping Coulomb scattering reaction.')
            except:
                raise
        for reaction in self.orphanProducts : reaction.calculateAverageProductData( style, indent = indent2, **kwargs )
        for reaction in additionalReactions : reaction.calculateAverageProductData( style, indent = indent2, **kwargs )

        if institutionModule.photoAtomicIncoherentDoppler in self.applicationData:
            MT504 = self.getReaction(504)
            totalCrossSection = MT504.crossSection[0]

            photon = MT504.outputChannel.products[IDsPoPsModule.photon]
            totalAverageProductEnergy = photon.averageProductEnergy[style.label]
            totalAverageProductMomentum = photon.averageProductMomentum[style.label]

            reactions = self.applicationData[institutionModule.photoAtomicIncoherentDoppler][0]
            for reaction in reactions:
                crossSection = reaction.crossSection[0]
                MT504_data = {}

                averageProductEnergy = crossSection * totalAverageProductEnergy / totalCrossSection
                averageProductEnergy = averageProductEnergy.thin(2e-3)
                MT504_data['averageProductEnergy'] = averageProductEnergy

                averageProductMomentum = crossSection * totalAverageProductMomentum  / totalCrossSection
                averageProductMomentum = averageProductMomentum.thin(2e-3)
                MT504_data['averageProductMomentum'] = averageProductMomentum

                kwargs['MT504_data'] = MT504_data
                reaction.calculateAverageProductData(style, indent=indent2, **kwargs)
                del kwargs['MT504_data']

    def buildCrossSectionSum( self, sumLabel, ENDF_MT, styleLabel, reactionKeys, Q ) :
        """
        For the list of indices and/or labels in reactionKeys, calcualate the summed crosstions and create a crossSectionSum
        for it. If insert is True, the crossSectionSum is added to self's sums/crossSectionSums node. If N keys match the
        same reaction in this, that reaction will be added N times.

        :param sumLabel: Label for the constructed crossSectionSum node
            :type sumLabel: str
        :param ENDF_MT: The ENDF MT for the crossSectionSum node
            :type ENDF_MT: int
        :param styleLabel: The label for the function1d added to the Q and crossSection nodes
            :type styleLabel: str
        :param reactionKeys: All reactions with an index or label matching an item of this list are added to cross section sum.
            :type reactionKeys: iterator, list or tuple
        :param Q: The Q value for the cross section sum
            :type Q: float

        :return: A crossSectionSum instance representing the requested reaction cross section sums
        :rtype: crossSectionSum instance
        """

        crossSectionSum = sumsModule.CrossSectionSum( sumLabel, ENDF_MT )

        sum = XYs1dModule.XYs1d( axes = crossSectionModule.defaultAxes( self.domainUnit ) )
        for key in reactionKeys :
            reaction = self.reactions[key]
            sum += reaction.crossSection.toPointwise_withLinearXYs( lowerEps = 1e-6, upperEps = 1e-6 )
            crossSectionSum.summands.append( sumsModule.Add( link = reaction.crossSection ) )

        crossSectionSum.crossSection.add( crossSectionModule.XYs1d( data = sum, label = styleLabel, axes = crossSectionModule.defaultAxes( self.domainUnit ) ) )

        crossSectionSum.Q.add( QModule.Constant1d( Q, domainMin = sum.domainMin, domainMax = sum.domainMax, axes = QModule.defaultAxes( self.domainUnit ), label = styleLabel ) )

        return( crossSectionSum )

    def inputParticlesToReactionString( self, prefix = "", suffix = "" ) :

        return( "%s%s + %s%s" % ( prefix, self.projectile, self.target, suffix ) )

    def cullProcessedData(self, include_reconstructed_data=False):
        """
        This method removes processed data from *self*.
        If include_reconstructed_data is True, then reconstructed data is also removed.
        """

        from fudge import institution as institutionModule
        from fudge.resonances import probabilityTables as probabilityTablesModule
        from fudge.processing.deterministic import tokens as deterministicTokensModule

        preProcessedStyles = self.styles.preProcessingStyles(not include_reconstructed_data)

        stylesToRemove = []
        for style in self.styles:
            if not isinstance(style, preProcessedStyles):
                stylesToRemove.append(style.label)
        self.removeStyles(stylesToRemove)

        applicationDataLabels = (deterministicTokensModule.multiGroupReactions,
            deterministicTokensModule.multiGroupDelayedNeutrons,
            deterministicTokensModule.multiGroupIncompleteProducts,
            probabilityTablesModule.LLNLProbabilityTablesToken)

        for label in applicationDataLabels: # applicationData requires special handling.
            if label in self.applicationData:
                self.applicationData.pop(label)

        try:
            photoAtomicIncoherentDoppler = self.applicationData[institutionModule.photoAtomicIncoherentDoppler]
            photoAtomicIncoherentDoppler.data[0].removeStyles(stylesToRemove)
        except:
            pass

    def cullStyles( self, styleName, removeDocumentation, removeApplicationData ) :
        """
        Removes all forms from each component except for the requested style or, if it is not present
        its nearest derived from style.
        """

        style = self.__styles[styleName]
        if( isinstance( style, stylesModule.HeatedMultiGroup ) ) :
            styleList = [ style ]
        elif( isinstance( style, stylesModule.GriddedCrossSection ) ) :
            styleList = []
            while( style is not None ) :
                styleList.append( style )
                style = style.derivedFromStyle
        else :
            raise TypeError( 'style name "%s" is for type "%s" which is not supported for culling.' % ( styleName, style.moniker ) )

        self.__reactions.cullStyles( styleList )
        self.__orphanProducts.cullStyles( styleList )
        self.__sums.cullStyles( styleList )
        self.__productions.cullStyles( styleList )
        self.__incompleteReactions.cullStyles( styleList )
        self.__fissionComponents.cullStyles( styleList )
# FIXME What about self.resonances

    def CoulombPlusNuclearMuCutoffs( self ) :

        muCutoffs = []
        for style in self.styles :
            if( isinstance( style, stylesModule.CoulombPlusNuclearElasticMuCutoff ) ) : muCutoffs.append( style.muCutoff )

        return( list( set( muCutoffs ) ) )

    def hasFission( self ) :

        for reaction in self.__reactions:
            if reaction.isFission(): return True

        return False

    def listOfProducts(self, finalOnly=False, includeQualifier=True):
        '''
        Returns, as a python set, the list of PoPs ids for all products (i.e., outgoing particles) for *self*.

        :param finalOnly:           If `True`, only final product ids are returned, otherwise, all are returned..
        :param includeQualifier:    If `True`, particle qualifiers are include in ids, otherwise, they are stripped from ids.

        :return:                    A python set instance.
        '''

        products = set()
        for reaction in self.reactions:
            products.update(reaction.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))
        for reaction in self.orphanProducts:
            products.update(reaction.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))
        for reaction in self.fissionComponents:
            products.update(reaction.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))
        for reaction in self.productions:
            products.update(reaction.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))
        for reaction in self.incompleteReactions:
            products.update(reaction.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))

        return products

    def addCovariance(self, covarianceSuite):
        '''
        Addes a CovarianceSuite instance to *self*'s internal *_loadedCovariances* member.

        :param covarianceSuite:     CovarianceSuite instance to add to member *_loadedCovariances*.
        '''

        from . import externalFile as externalFileModule
        from fudge.covariances import covarianceSuite as covarianceSuiteModule

        if not isinstance(covarianceSuite, covarianceSuiteModule.CovarianceSuite):
            raise TypeError('Instance not a CovarianceSuite instance.')

        if self._loadedCovariances is None:
            self._loadedCovariances = []

        defaultPath = 'dummy%s.xml'
        covariancePath = covarianceSuite.sourcePath 

        if covariancePath is None:
            covariancePath = defaultPath % '-covar'
            covarianceSuite.sourcePath = covariancePath

        selfPath = self.sourcePath
        if selfPath is None:
            selfPath = defaultPath % ''
            self.__sourcePath = selfPath

        covarianceSuite.externalFiles.add(externalFileModule.ExternalFile('reactions', selfPath))
        self.externalFiles.add(externalFileModule.ExternalFile('covariances', covariancePath))

        self._loadedCovariances.append(covarianceSuite)

    def covarianceExternalFiles(self):
        '''
        Returns a list of all ExternalFile instanaces in *self* that point to covariance files.

        :return:    The list of ExternalFile instanaces that point to covariance files.
        '''

        from fudge import GNDS_file as GNDS_fileModule
        from fudge.covariances import covarianceSuite as covarianceSuiteModule

        externalFiles = []
        for externalFile in self.externalFiles:
            # TODO submit proposal to include file type attribute in the externalFile GNDS node so we don't need to run GNDS_fileModule.type
            if not os.path.exists(externalFile.realpath()):
                print("  WARNING: skipping externalFile '%s' which points to non-existing file %s." %
                      (externalFile.label, externalFile.realpath()))
                continue
            gndsType, metadata = GNDS_fileModule.type(externalFile.realpath())
            if gndsType == covarianceSuiteModule.CovarianceSuite.moniker:
                externalFiles.append(externalFile)

        return externalFiles

    def loadCovariances(self, skipConvertingUnits=False):
        """
        Load all external files of type 'covarianceSuite', and resolve links between self and covarianceSuites.

        :param skipConvertingUnits:  If False (default), convert all covarianceSuites to have the same domain unit as self.
        :return: list of loaded covarianceSuites
        """

        from fudge.covariances import covarianceSuite as covarianceSuiteModule

        if self._loadedCovariances is None:
            kwargs = {'reactionSuite': self}
            self._loadedCovariances = []
            for covariance in self.covarianceExternalFiles():
                covariance.instance = covarianceSuiteModule.read(covariance.realpath(), **kwargs)
                if not skipConvertingUnits:
                    covariance.instance.convertUnits({covariance.instance.domainUnit: self.domainUnit})
                self._loadedCovariances.append(covariance.instance)

        return self._loadedCovariances

    def partialProductionIntegral( self, productID, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.lab,
                LegendreOrder = 0, **kwargs ) :

        partialProductionIntegralSum = 0.0

        for reaction in self.reactions :
            partialProductionIntegralSum += reaction.partialProductionIntegral( self, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                    phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )
        for orphanProduct in self.orphanProducts :
            partialProductionIntegralSum += orphanProduct.partialProductionIntegral( self, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                    phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )

        return( partialProductionIntegralSum )

    def processCoulombPlusNuclearMuCutoff( self, style, excludeRutherfordScattering = False ) :

        nuclearPlusCoulombInterference = None

        for reaction in self.reactions :
            RutherfordScatteringExcluded = reaction.processCoulombPlusNuclearMuCutoff( style, excludeRutherfordScattering = excludeRutherfordScattering )

            if( RutherfordScatteringExcluded is not None ) :
                crossSection, function2d = RutherfordScatteringExcluded
                firstProduct = reaction.outputChannel[0]
                secondProduct = reaction.outputChannel[1]

                nuclearPlusCoulombInterference = nuclearPlusCoulombInterferenceModule.NuclearPlusCoulombInterference(style.label)

                try :
                    LLNL = self.applicationData['LLNL']
                except :
                    LLNL = institutionModule.Institution('LLNL')
                    self.applicationData.add( LLNL )
                LLNL.append( nuclearPlusCoulombInterference )

                nuclearPlusCoulombInterference.reaction.crossSection.add( crossSection )

                Q = reaction.outputChannel.Q[0].copy( )
                nuclearPlusCoulombInterference.reaction.outputChannel.Q.add( Q )

                product1 = productModule.Product(firstProduct.pid, firstProduct.label)
                product1.multiplicity.add( firstProduct.multiplicity[0].copy( ) )
                angularTwoBody1 = angularModule.TwoBody( crossSection.label, xDataEnumsModule.Frame.centerOfMass, function2d )
                product1.distribution.add( angularTwoBody1 )
                nuclearPlusCoulombInterference.reaction.outputChannel.products.add( product1 )

                product2 = productModule.Product(secondProduct.pid, secondProduct.label)
                product2.multiplicity.add( secondProduct.multiplicity[0].copy( ) )
                recoil = angularModule.Recoil( link = product1.distribution[style.label], relative = True )
                angularTwoBody2 = angularModule.TwoBody( crossSection.label, xDataEnumsModule.Frame.centerOfMass, recoil )
                product2.distribution.add( angularTwoBody2 )
                nuclearPlusCoulombInterference.reaction.outputChannel.products.add( product2 )

                nuclearPlusCoulombInterference.reaction.outputChannel.process = 'no Rutherford'
                nuclearPlusCoulombInterference.reaction.updateLabel( )

        return( nuclearPlusCoulombInterference )

    def processThermalNeutronScatteringLaw( self, style, indent = '', verbosity = 0 ) :

        kwargs = { 'reactionSuite' : self }
        kwargs['product'] = self.PoPs[IDsPoPsModule.neutron]
        kwargs['incidentEnergyUnit'] = self.domainUnit
        kwargs['temperatureUnit'] = kwargs['incidentEnergyUnit'] + '/k'
        kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'
        kwargs['momentumUnit'] = kwargs['incidentEnergyUnit'] + '/c'
        kwargs['neutronMass'] = self.PoPs[self.projectile].getMass( kwargs['massUnit'] )
        kwargs['productMass'] = kwargs['neutronMass']
        kwargs['accuracy'] = 0.001
        kwargs['epsilon'] = 1e-6
        kwargs['verbosity'] = verbosity

        energyMax = self.styles.getEvaluatedStyle().projectileEnergyDomain.max
        if( energyMax > 0.0 ) : kwargs['energyMax'] = energyMax

        for reaction in self.reactions : reaction.processThermalNeutronScatteringLaw( style, kwargs )

    def processMC_cdf( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', additionalReactions = [] ) :

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))

        tempInfo = { 'reactionSuite' : self }
        tempInfo['verbosity'] = verbosity
        tempInfo['incrementalIndent'] = incrementalIndent
        tempInfo['brokenLinks'] = []

        for reaction in self.reactions : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )
        for reaction in self.orphanProducts : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )
        for reaction in additionalReactions : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )

        if institutionModule.photoAtomicIncoherentDoppler in self.applicationData:
            reactions = self.applicationData[institutionModule.photoAtomicIncoherentDoppler][0]
            for reaction in reactions:
                reaction.processMC_cdf(style, tempInfo, indent + incrementalIndent)

        for original, newReference in tempInfo['brokenLinks']:
            newReference.link = original.referenceInstance.ancestor[ newReference.label ]

    def processGriddedCrossSections( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', additionalReactions = [] ) :
        """
        Generate a union grid for all cross sections, for faster continuous energy Monte Carlo sampling.
        The grid spans the same domain as the elastic scattering cross section.
        """

        def mutualDomain( style, reaction, thresholds ) :

            crossSection = style.findFormMatchingDerivedStyle( reaction.crossSection )
            if( isinstance( crossSection, crossSectionModule.Regions1d ) ) :
                crossSection = crossSection.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-6 )
            crossSection = crossSection.domainSlice( domainMin, domainMax )                 # truncate to match elastic domain
            if( ( crossSection.domainMin > domainMin ) and ( crossSection[0][1] != 0 ) ) :
                print('    WARNING, domains not mutual, setting first y-value to 0.')
                crossSection[0] = [ crossSection[0][0], 0 ]
            if( ( crossSection.domainMax < domainMax ) and ( crossSection[-1][1] != 0 ) ) :
                print('    WARNING, domains not mutual, setting last y-value to 0.')
                crossSection[-1] = [ crossSection[-1][0], 0 ]
            thresholds.add( crossSection[0][0] )
            return( crossSection )

        def photoAtomic( style, reaction, thresholds ) :

            crossSection = style.findFormMatchingDerivedStyle( reaction.crossSection )
            if( isinstance( crossSection, crossSectionModule.Regions1d ) ) :
                values = []
                x1, y1 = crossSection[0][0]
                for iRegion, region in enumerate( crossSection ) :
                    domainGrid = region.domainGrid
                    values += domainGrid
            else :
                x1, y1 = crossSection[0]
                values = crossSection.domainGrid
            if( ( x1 > domainMin ) and ( y1 != 0 ) ) :
                xEpsilon = 1e-6 * x1 + x1
                if( xEpsilon < values[1] ) : values.insert( 1, xEpsilon )

            return( set( values ) )

        evaluationStyle = style.findDerivedFromStyle( stylesModule.Evaluated )
        projectileEnergyDomain = evaluationStyle.projectileEnergyDomain
        domainMin = projectileEnergyDomain.min
        domainMax = projectileEnergyDomain.max

        axes = style.findFormMatchingDerivedStyle( self.reactions[0].crossSection ).axes.copy( )

        isPhotoAtomic = False
        if( not( self.isThermalNeutronScatteringLaw( ) ) ) :
            target = self.PoPs[self.target]
            isChemicalElement = isinstance( target, chemicalElementPoPsModule.ChemicalElement )
            isUnorthodox = isinstance( target, unorthodoxPoPsModule.Particle )
            isPhotoAtomic = ( self.projectile == IDsPoPsModule.photon ) and ( isChemicalElement or isUnorthodox )

        thresholds = set( )
        total = crossSectionModule.XYs1d( data = [], axes = axes )
        if( isPhotoAtomic ) :
            values = set( )
            for reaction in self.reactions : values.update( set( photoAtomic( style, reaction, thresholds ) ) )
            values = sorted( values )
        else :
            for reaction in self.reactions : total += mutualDomain( style, reaction, thresholds )
            for reaction in self.orphanProducts : total += mutualDomain( style, reaction, thresholds )
            values = total.domainGrid

        tooCloseXsIndicies = []
        x1 = -1
        for i1, x2 in enumerate( values ) :
            if( abs( x2 - x1 ) < 1e-12 * x2 ) : tooCloseXsIndicies.append( i1 )
            x1 = x2
        tooCloseXsIndicies.reverse( )
        for index in tooCloseXsIndicies :                   # Remove close points but any threshold points.
            if( values[index] in thresholds ) : index -= 1
            del values[index]

        values = valuesModule.Values( values )
        grid = axesModule.Grid( total.axes[1].label, 1, total.axes[1].unit, xDataEnumsModule.GridStyle.points, values )
        style.grid = grid

        for reaction in self.reactions :
            reaction.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent, isPhotoAtomic = isPhotoAtomic )

        if( not( isPhotoAtomic ) ) :
            for reaction in self.orphanProducts :
                reaction.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

            for production in self.productions :
                production.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

            for sum in self.sums.crossSectionSums :
                sum.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

            for fissionComponent in self.fissionComponents :
                fissionComponent.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

        for reaction in additionalReactions :
            reaction.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

        if institutionModule.photoAtomicIncoherentDoppler in self.applicationData:
            reactions = self.applicationData[institutionModule.photoAtomicIncoherentDoppler][0]
            for reaction in reactions:
                reaction.processGriddedCrossSections(style, verbosity=verbosity, indent=indent, incrementalIndent=incrementalIndent)

    def processMultiGroup( self, style, legendreMax, verbosity = 0, indent = '', incrementalIndent = '  ', logFile = None, workDir = None,
                restart = False, additionalReactions = [] ) :
        """
        Generate multi-group processed data, including cross sections, weighted multiplicities and Q-values and transfer matrices.
        Processed results will be stored as new forms inside the reactionSuite.

        :param style: fudge.styles.HeatedMultiGroup instance
        :param legendreMax: number of Legendre orders to produce in transfer matrices.
        :param verbosity: verbosity level
        :param indent: indentation level for verbose output
        :param incrementalIndent: amount to increase verbose output indentation when calling processMultiGroup on child classes
        :param logFile: open file where log will be written.  FIXME logFile is a required argument, but default value = None
        :param workDir: directory to save working files generated during processing. Note that this directory can become very large.
        :param restart: load previous Merced outputs (if found) rather than recomputing, useful if a processProtare run times out or crashes.
        :return: None, except ValueError is raised if any errors occurred during processing
        """

        from LUPY import times as timesModule

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))
        if( not( isinstance( style, stylesModule.HeatedMultiGroup ) ) ) : raise( 'Instance is not a HeatedMultiGroup style.' )

        t0 = timesModule.Times( )

        kwargs = { 'reactionSuite' : self }
        kwargs['failures'] = 0
        if( len( self.reactions ) > 0 ) :
            kwargs['verbosity'] = verbosity
            kwargs['incrementalIndent'] = incrementalIndent
            kwargs['logFile'] = logFile
            kwargs['incidentEnergyUnit'] = self.domainUnit
            kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'
            kwargs['restart'] = restart
            kwargs['legendreMax'] = legendreMax

            projectile = self.PoPs[self.projectile]
            kwargs['projectile'] = projectile
            kwargs['projectileZA'] = chemicalElementMiscPoPsModule.ZA( projectile )
            kwargs['projectileMass'] = projectile.getMass( kwargs['massUnit'] )

            kwargs['targetMass'] = None
            if( not( self.isThermalNeutronScatteringLaw( ) ) ) :
                target = self.PoPs[self.target]
                kwargs['isInfiniteTargetMass'] = isinstance( target, chemicalElementPoPsModule.ChemicalElement )
                if( not( isinstance( target, chemicalElementPoPsModule.ChemicalElement ) ) ) :
                    if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.pid]
                    kwargs['targetMass'] = target.getMass( kwargs['massUnit'] )
                kwargs['target'] = target
                kwargs['targetZA'] = chemicalElementMiscPoPsModule.ZA( target )

            kwargs['masses'] = { 'Projectile'   : kwargs['projectileMass'],
                                 'Product'      : None,
                                 'Residual'     : None,
                                 'Target'       : kwargs['targetMass'] }

            kwargs['zeroPerTNSL'] = self.isThermalNeutronScatteringLaw( )
            projectileEnergyDomain = style.projectileEnergyDomain
            transportable = style.transportables[self.projectile]
            for index, boundary in enumerate( transportable.group.boundaries.values ) :
                if( boundary > projectileEnergyDomain.max ) : break
                kwargs['maximumProjectileGroupIndex'] = index

            if( workDir is None ) : workDir = 'Merced.work'
            kwargs['workDir'] = workDir
            kwargs['workFile'] = []

            style.processMultiGroup( style, kwargs, indent + incrementalIndent )
            kwargs['groupedFlux'] = [ x for x in style.multiGroupFlux.array.constructArray( )[:,0] ]

# BRB FIXME, must have test to determine if reconstructResonances is needed.
#        self.reconstructResonances( styleName = 'reconstructed', 1e-3, verbose = False )
            for i1, reaction in enumerate( self.reactions ) :
                kwargs['reactionIndex'] = "%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.orphanProducts ) :
                kwargs['reactionIndex'] = "o%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.productions ) :
                kwargs['reactionIndex'] = "p%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.sums.crossSectionSums ) :
                kwargs['reactionIndex'] = "s%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.fissionComponents ) :
                kwargs['reactionIndex'] = "f%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( additionalReactions ) :
                kwargs['reactionIndex'] = "a%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )

            if institutionModule.photoAtomicIncoherentDoppler in self.applicationData:
                reactions = self.applicationData[institutionModule.photoAtomicIncoherentDoppler][0]
                for reaction in reactions:
                    reaction.processMultiGroup(style, kwargs, indent + incrementalIndent)

        if logFile is not None: logFile.write( '    ' + str( t0 ) + '\n' )
        if( kwargs['failures'] > 0 ) : raise ValueError( "kwargs['failures'] = %d  > 0" % kwargs['failures'] )

    def multiGroupVector(self, childSuite, multiGroupVectorMethod, multiGroupSettings, temperatureInfo, **kwargs):
        r"""
        General method to return multiGroupVector output from a specified multiGroupVectorMethod for specified reactionSuite child nodes.

        :param childSuite: Suite which to iterate over.
        :param multiGroupVectorMethod: Method to be executed for each entry in the suite of reactionSuite child nodes.
        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param \*\*kwargs: Additional arguments for the selected method in the entries in the selected suite of reactionSuite child nodes.
        """

        vectorResult = vectorModule.Vector()
        for suiteEntry in childSuite:
            vectorMethod = getattr(suiteEntry, multiGroupVectorMethod)
            vectorResult += vectorMethod(multiGroupSettings, temperatureInfo, **kwargs)
        
        return vectorResult

    def multiGroupMatrix(self, childSuite, multiGroupMatrixMethod, multiGroupSettings, temperatureInfo, **kwargs):
        """
        General method to return multiGroupMatrix output from a specified multiGroupVectorMethod for specified reactionSuite child nodes.

        :param childSuite: Suite which to iterate over.
        :param multiGroupMatrixMethod: Method to be executed for each entry in the suite of reactionSuite child nodes.
        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param kwargs: Additional arguments for the selected method in the entries in the selected suite of reactionSuite child nodes.
        """

        matrixResult = matrixModule.Matrix()
        for suiteEntry in childSuite:
            matrixMethod = getattr(suiteEntry, multiGroupMatrixMethod)
            matrixResult += matrixMethod(multiGroupSettings, temperatureInfo, **kwargs)
        
        return matrixResult

    def multiGroupCrossSection(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, total cross section for the requested label. This is summed over all reactions.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        return self.multiGroupVector(self.reactions, 'multiGroupCrossSection', multiGroupSettings, temperatureInfo)

    def multiGroupInverseSpeed(self, temperatureInfo):
        """
        Returns the inverse speeds for the requested label. The label must be for a heated multi-group style.

        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        styleLabel = temperatureInfo.heatedMultiGroup
        return self.styles[styleLabel].inverseSpeed.data.constructVector()

    def multiGroupQ(self, multiGroupSettings, temperatureInfo, includeFinalProducts):
        """
        Returns the multi-group, total Q for the requested label. This is a cross section weighted multiplicity summed over all reactions.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param includeFinalProducts: If true, the Q is calculated for all output channels, including those for products.
        """

        return self.multiGroupVector(self.reactions, 'multiGroupQ', multiGroupSettings, temperatureInfo, includeFinalProducts=includeFinalProducts)

    def multiGroupMultiplicity(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the multi-group, total multiplicity for the requested label for the requested product. This is a cross section weighted multiplicity.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        multiplicity = self.multiGroupVector(self.reactions, 'multiGroupMultiplicity', multiGroupSettings, temperatureInfo, productID=productID)
        multiplicity += self.multiGroupVector(self.orphanProducts, 'multiGroupMultiplicity', multiGroupSettings, temperatureInfo, productID=productID)

        return multiplicity

    def multiGroupFissionNeutronMultiplicity(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, total fission neutron multiplicity for the requested label. This is a cross section weighted multiplicity.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        multiplicity = vectorModule.Vector()
        for reaction in self.reactions:
            if reaction.isFission():
                multiplicity += reaction.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, IDsPoPsModule.neutron)

        return multiplicity

    def multiGroupFissionGammaMultiplicity(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, total fission photon multiplicity for the requested label. This is a cross section weighted multiplicity.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        multiplicity = vectorModule.Vector()
        for reaction in self.reactions:
            if reaction.isFission():
                multiplicity += reaction.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, IDsPoPsModule.photon)

        return multiplicity

    def multiGroupAvailableEnergy(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, total available energy for the requested label. This is a cross section weighted available energy summed over all reactions.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        return self.multiGroupVector(self.reactions, 'multiGroupAvailableEnergy', multiGroupSettings, temperatureInfo)

    def multiGroupAverageEnergy(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the multi-group, total average energy for the requested label for the requested product. This is a cross section weighted average energy summed over all reactions.
        
        This is a cross section weighted average energy summed over all products for this reaction.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        averageEnergy = self.multiGroupVector(self.reactions, 'multiGroupAverageEnergy', multiGroupSettings, temperatureInfo, productID=productID)
        averageEnergy += self.multiGroupVector(self.orphanProducts, 'multiGroupAverageEnergy', multiGroupSettings, temperatureInfo, productID=productID)

        return averageEnergy

    def multiGroupDepositionEnergy(self, multiGroupSettings, temperatureInfo, particleIDs):
        """
        Returns the multi-group, total deposition energy for the requested label.
        
        This is a cross section weighted deposition energy summed over all reactions. The deposition energy is calculated by subtracting 
        the average energy from each transportable particle from the available energy. The list of transportable particles is specified 
        via the list of particle specified in the *particleIDs* argument.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        """

        depositionEnergy = self.multiGroupAvailableEnergy(multiGroupSettings, temperatureInfo)
        for particleID in particleIDs:
            depositionEnergy -= self.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, particleID)

        return depositionEnergy 

    def multiGroupAvailableMomentum(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, total available momentum for the requested label.
        
        This is a cross section weighted available momentum summed over all reactions.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        return self.multiGroupVector(self.reactions, 'multiGroupAvailableMomentum', multiGroupSettings, temperatureInfo)

    def multiGroupAverageMomentum(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the multi-group, total average momentum for the requested label for the requested product.
        
        This is a cross section weighted average momentum summed over all reactions.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        averageMomentum = self.multiGroupVector(self.reactions, 'multiGroupAverageMomentum', multiGroupSettings, temperatureInfo, productID=productID)
        averageMomentum += self.multiGroupVector(self.orphanProducts, 'multiGroupAverageMomentum', multiGroupSettings, temperatureInfo, productID=productID)

        return averageMomentum

    def multiGroupDepositionMomentum(self, multiGroupSettings, temperatureInfo, particleIDs):
        """
        Returns the multi-group, total deposition momentum for the requested label.
        
        This is a cross section weighted deposition momentum summed over all reactions. The deposition momentum is calculated by 
        subtracting the average momentum from each transportable particle from the available momentum. The list of transportable 
        particles is specified via the list of particle specified in the *particleIDs* argument.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        """

        depositionMomentum = self.multiGroupAvailableMomentum(multiGroupSettings, temperatureInfo)
        for particleID in particleIDs:
            depositionMomentum -= self.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, particleID)

        return depositionMomentum

    def multiGroupGain(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the multi-group, gain for the requested particle and label.
        
        This is a cross section weighted gain summed over all reactions.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        particleGain = self.multiGroupVector(self.reactions, 'multiGroupGain', multiGroupSettings, temperatureInfo, \
            productID=productID, projectileID=self.projectile)

        return particleGain

    def multiGroupTransportCorrection(self, multiGroupSettings, temperatureInfo, particleIDs, legendreOrder, \
        temperature, transportCorrectionType=transportingModule.TransportCorrectionType.none):
        """
        Returns the multi-group transport correction for the requested label.
        
        The transport correction is calculated from the transfer matrix for the projectile id for the Legendre order of legendreOrder + 1.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param legendreOrder: Requested Legendre order.
        :param temperature: The temperature of the flux to use when collapsing.
        :param transportCorrectionType: Requested transport correction type.
        """

        if transportCorrectionType == transportingModule.TransportCorrectionType.Pendlebury:
            multiGroupProductMatrix = self.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, \
                self.projectile, legendreOrder + 1)

            matrixCollapsed = particlesModule.CollapseMatrix(multiGroupProductMatrix, multiGroupSettings, particleIDs, temperature, self.projectile)
            transportCorrection = vectorModule.Vector(values=matrixCollapsed.matrix.diagonal())

        elif transportCorrectionType is transportingModule.TransportCorrectionType.none:
            transportCorrection = vectorModule.Vector()

        else:
            raise TypeError('Unsupported transport correction: only None and Pendlebury (i.e., Pendlebury/Underhill) are currently supported.')

        return transportCorrection

    def multiGroupProductMatrix(self, multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder):
        """
        Returns the multi-group, total product matrix for the requested label for the requested product id for the requested Legendre order.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        :param legendreOrder: Requested Legendre order.
        """

        productMatrix = self.multiGroupMatrix(self.reactions, 'multiGroupProductMatrix', multiGroupSettings, temperatureInfo, \
            particleIDs=particleIDs, productID=productID, legendreOrder=legendreOrder)

        productMatrix += self.multiGroupMatrix(self.orphanProducts, 'multiGroupProductMatrix', multiGroupSettings, temperatureInfo, \
            particleIDs=particleIDs, productID=productID, legendreOrder=legendreOrder)

        return productMatrix

    def multiGroupFissionMatrix(self, multiGroupSettings, temperatureInfo, particleIDs, legendreOrder):
        fissionMatrix = self.multiGroupMatrix(self.reactions, 'multiGroupFissionMatrix', multiGroupSettings, temperatureInfo, \
            particleIDs=particleIDs, legendreOrder=legendreOrder)

        return fissionMatrix

    def multiGroupProductArray(self, multiGroupSettings, temperatureInfo, particleIDs, productID):
        """
        Returns the full multi-group, total product array for the requested label for the requested product id.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        """

        productArray = productArrayModule.ProductArray()

        for reaction in self.reactions:
            productArray += reaction.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

        for reaction in self.orphanProducts:
            productArray += reaction.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

        return productArray

    def addMultiGroupSums(self, replace=False):

        from .processing.deterministic import addMultiGroupSums as addMultiGroupSumsModule

        addMultiGroupSumsModule.addMultiGroupSums(self, replace=replace)

    def processSnElasticUpScatter(self, style, legendreMax, verbosity = 0, indent = '',
                                  incrementalIndent = '  ', logFile = None, workDir = None,
                                  restart = False):

        if self.target == IDsPoPsModule.neutron: return

        if not isinstance(style, stylesModule.SnElasticUpScatter):
            raise TypeError("style must be an instance of stylesModule.SnElasticUpScatter, not %s"
                            % type(style))

        tempInfo = {
            'reactionSuite' : self,
            'verbosity': verbosity,
            'incrementalIndent': incrementalIndent,
            'logFile': logFile,
            'temperature': style.temperature,
            'legendreOrder': legendreMax,
            'zeroPerTNSL': False,
            'workDir': workDir,
            'restart': restart,
        }

        incidentEnergyUnit = self.domainUnit
        tempInfo['minEval'] = PQUModule.PQU(1.e-11, 'MeV').getValueAs(incidentEnergyUnit)
        tempInfo['maxEval'] = PQUModule.PQU(20., 'MeV').getValueAs(incidentEnergyUnit)

        massUnit = incidentEnergyUnit + '/c**2'
        target = self.PoPs.final(self.target)
        tempInfo['targetMassRatio'] = target.getMass(massUnit) / self.PoPs[self.projectile].getMass(massUnit)
        tempInfo['groupBoundaries'] = style.transportables[self.projectile].group

        elasticReaction = self.getReaction('elastic')
        product = elasticReaction.outputChannel.getProductWithName(IDsPoPsModule.neutron)

        TM_1, TM_E, averageEnergy, maxIncidentGroup = multiGroupUpScatterModule.SnElasticUpScatter(
            style, tempInfo)
        style.upperCalculatedGroup = maxIncidentGroup

                # multiply constant cross section into the distribution (upscatter code assumes sigma(E) = 1 )
        groupedCrossSec = product.findAttributeInAncestry('crossSection')[style.derivedFromStyle.label].array.values
        TM_1 = multiGroupUpScatterModule.rescaleCrossSection(groupedCrossSec, TM_1)

        distribution = product.distribution                             # Put matrix into GNDS structures.
        derivedForm = distribution[style.derivedFrom]
        multiGroupSubform = derivedForm.multiGroupSubform
        array = multiGroupSubform.array.constructArray()
        for incidentGroup in range(style.upperCalculatedGroup + 1, len(array)):
            if incidentGroup not in TM_1: TM_1[incidentGroup] = {}
            for outgoingGroup in range(len(array[incidentGroup])):
                if outgoingGroup not in TM_1[incidentGroup]: TM_1[incidentGroup] = {}
                TM_1[incidentGroup][outgoingGroup] = list(array[incidentGroup][outgoingGroup])

        tempInfo['productName'] = product.pid
        tempInfo['groupedFlux'] = [x for x in style.derivedFromStyle.multiGroupFlux.array.constructArray()[:,0]]
        tempInfo['incidentEnergyUnit'] = incidentEnergyUnit
        tempInfo['crossSection']  = style.derivedFromStyle.findFormMatchingDerivedStyle(elasticReaction.crossSection)
        multiGroup = groupModule.TMs2Form(style, tempInfo, TM_1, TM_E)
        distribution.add(multiGroup)

        averageProductEnergy = product.averageProductEnergy
        axes = averageProductEnergyModule.defaultAxes(energyUnit = incidentEnergyUnit)
        averageEnergy = averageProductEnergyModule.XYs1d(data = averageEnergy, axes = axes, label = style.label)
        averageEnergyMultiGroup = averageEnergy.processMultiGroup(style, tempInfo, indent)
        grouped = averageEnergyMultiGroup.array.constructArray()
        array = averageProductEnergy[style.derivedFrom].array.constructArray()
        for incidentGroup in range(style.upperCalculatedGroup + 1, len(array)): grouped[incidentGroup] = array[incidentGroup]
        averageProductEnergy.add(groupModule.toMultiGroup1d(averageProductEnergyModule.Gridded1d, style, 
                tempInfo, averageEnergyMultiGroup.axes, grouped, zeroPerTNSL=tempInfo['zeroPerTNSL']))

        return averageEnergy

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__styles.removeStyles( styleLabels )
        self.__reactions.removeStyles( styleLabels )
        self.__orphanProducts.removeStyles( styleLabels )
        self.__productions.removeStyles( styleLabels )
        self.__incompleteReactions.removeStyles( styleLabels )
        self.__fissionComponents.removeStyles( styleLabels )
        self.__sums.removeStyles( styleLabels )

    def supportsResonanceReconstruction( self ):
        """
        Checks whether self.resonances contains any resonance parameters that can be reconstructed into
        cross sections and/or product distributions.
        :return: boolean
        """
        if self.resonances is None: return False
        return self.resonances.reconstructCrossSection

    def reconstructResonances( self, style, accuracy, thin = True, significantDigits = None, verbose = False ):
        """
        Turn resonance parameters into pointwise cross sections, then merge the results with
        tabulated pointwise cross sections. Resulting pointwise cross sections are stored
        alongside the original 'resonancesWithBackground' data in the reactionSuite.

        :param style: style - an instance of style CrossSectionReconstructed.
        :param accuracy: float - target accuracy during reconstruction. For example, 0.001.
        :param thin: boolean - enable/disable thinning after reconstruction (disabling makes summed cross sections consistency checks easier).
        :param significantDigits: int - energy grid will only include points that can be exactly represented using specified number of significant digits.
        :param verbose: boolean - turn on/off verbosity
        """

        if( self.resonances is None ) : return
        if not self.resonances.reconstructCrossSection:
            return # nothing to do
        from fudge.processing.resonances import reconstructResonances

        if not isinstance( style, stylesModule.CrossSectionReconstructed ):
            raise TypeError("style must be an instance of CrossSectionReconstructed, not %s" % type(style))

        xsecs = reconstructResonances.reconstructResonances(self, accuracy, verbose = verbose, significantDigits = significantDigits )
        epsilon = 1e-8  # for joining multiple regions together

        # for each reaction, add tabulated pointwise data (ENDF MF=3) to reconstructed resonances:
        possibleChannels = { 'elastic' : False, 'capture' : True, 'fission' : True, 'total' : False, 'nonelastic' : False }
        derivedFromLabel = ''
        for reaction in self :
            if isinstance( reaction, sumsModule.MultiplicitySum ): continue
            evaluatedCrossSection = reaction.crossSection.evaluated
            if not isinstance( evaluatedCrossSection, crossSectionModule.ResonancesWithBackground ):
                continue
            # which reconstructed cross section corresponds to this reaction?
            if( derivedFromLabel == '' ) : derivedFromLabel = evaluatedCrossSection.label
            if( derivedFromLabel != evaluatedCrossSection.label ) :
                print ('WARNING derivedFromLabel = "%s" != "%s"' % (derivedFromLabel, evaluatedCrossSection.label))
            RRxsec = None
            if str( reaction ) in xsecs:
                RRxsec = xsecs[ str( reaction ) ]
            else :
                for possibleChannel in possibleChannels :
                    if( possibleChannels[possibleChannel] ) :
                        if( possibleChannel in str( reaction ) ) : RRxsec = xsecs[possibleChannel]
                    if( RRxsec is None ) :
                        if( reaction is self.getReaction( possibleChannel ) ) : RRxsec = xsecs[possibleChannel]
                    if( RRxsec is not None ) : break
            if RRxsec is None:
                found = False
                if isinstance(reaction, sumsModule.CrossSectionSum):
                    from fudge.reactions import base as baseModule
                    RRxsec = 0
                    for summand in reaction.summands:
                        rlabel = summand.link.findClassInAncestry(baseModule.Base_reaction).label
                        if rlabel in xsecs:
                            RRxsec += xsecs[rlabel]
                            found = True

                if not found:
                    raise ValueError("Couldn't find appropriate reconstructed cross section to add to reaction '%s'!"
                                     % reaction)

            background = evaluatedCrossSection.background
            try:
                background = background.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = epsilon, upperEps = epsilon )
            except Exception as ex:
                print("Encountered an error while linearizing background cross section for reaction '%s'!" % reaction.label)
                raise ex
            RRxsec = RRxsec.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = epsilon, upperEps = epsilon, removeOverAdjustedPoints=True )
            RRxsec.convertUnits( {RRxsec.domainUnit: background.domainUnit,  RRxsec.rangeUnit: background.rangeUnit } )

            background, RRxsec = background.mutualify(0,0,0, RRxsec, -epsilon,epsilon,True)
            RRxsec = background + RRxsec    # result is a crossSection.XYs1d instance
            if( RRxsec.rangeMin < 0 ) :
                # turn any negative xsc to 0
                RRxsec = RRxsec.clip( rangeMin=0 )
                if verbose:
                    print( "Warning: negative cross section encountered for %s; changed to 0 b" % reaction )
            if( RRxsec[0][1] == 0 ) :
                # special handling for threshold reaction:
                #  trim leading zeros + points with cross section < 1e-90 b,
                #  thin aggressively up to xsc = 1e-10 b
                firstNonZero = next((i for i, xy in enumerate(RRxsec) if xy[1]), None)
                RRxsec = RRxsec[firstNonZero-1:]

                points = RRxsec.copyDataToXYs()
                if points[0][0] > evaluatedCrossSection.domainMin:
                    points.insert(0, [evaluatedCrossSection.domainMin, 0])

                indices = []
                priorCrossSection = 0.0
                factor = 1e5
                for idx, xy in enumerate(points):
                    if idx==0: continue
                    if xy[1] > 1e-10: break
                    if xy[1] > 1e-20: factor = 10
                    if xy[1] < 1e-90 or xy[1] < factor * priorCrossSection:
                        indices.append(idx)
                    else:
                        priorCrossSection = xy[1]

                for idx in reversed(indices):
                    del points[idx]

                RRxsec.setData(points)

            if thin:
                RRxsec = RRxsec.thin( accuracy )
            RRxsec.label = style.label
            reaction.crossSection.add( RRxsec )

        self.styles.add( style )

    def reconstructResonancesAngularDistributions( self, styleName, overwrite=False, accuracy = None, thin = False, verbose = False ):
        """
        Turn resonance parameters into Legendre-expanded angular distributions, then merge the results with
        tabulated angular distributions (overwriting original values in the resonance region). Resulting pointwise
        angular distributions are stored alongside the original angular distributions in the reactionSuite.

        :param styleName: string - label for reconstructed distribution style
        :param overwrite: boolean - if style already in use for another distribution, overwrite it with the reconstructed distribution
        :param accuracy: float - giving target accuracy during reconstruction. For example, 0.001
        :param thin: boolean - enable/disable thinning after reconstruction
        :param verbose: boolean - turn on/off verbosity
        """

        if accuracy is not None: raise NotImplementedError("Refining interpolation grid for angular distributions still TBD")
        if thin: raise NotImplementedError("Thinning for angular distributions still TBD")

        if( self.resonances is None ) : return
        from fudge.processing.resonances import reconstructResonances

        newStyle = stylesModule.AngularDistributionReconstructed( styleName, 'eval' )       # BRB6, FIXME - 'eval' should not be hardwired.

        distributions = reconstructResonances.reconstructAngularDistributions( self, tolerance=accuracy, verbose=verbose )

        for key in distributions:
            reaction = self.getReaction( key )
            # Kludgy way to discern the product's name so we can look it up
# BRB6 hardwired.
            if key in 'elastic': productName = 'n'
            else:
# BRB6 hardwired.
                pairs = [ ( p.particle.getMass('amu'), p.pid ) for p in reaction.outputChannel.products ]
                productName = min(pairs)[1] # pick's the lightest product
            product = reaction.outputChannel.getProductWithName(productName)
            original = product.distribution.evaluated.angularSubform
            reconstructed = distributions[key]

            merged = angularModule.Regions2d( axes = original.axes )     # Caleb, FIXME - check that this is correct.
            merged.append( reconstructed )

            if isinstance( original, angularModule.XYs2d ):
                newregion = angularModule.XYs2d( axes = original.axes, interpolation=original.interpolation,
                        interpolationQualifier=original.interpolationQualifier )
                newregion.append( original.getAtEnergy( reconstructed.domainMax ) )
                for val in original:
                    if( val.outerDomainValue <= reconstructed.domainMax ) : continue
                    newregion.append( val )
                merged.append( newregion )
            elif isinstance( original, angularModule.Regions2d ):
                for idx,region in enumerate(original):
                    if( region.domainMax > reconstructed.domainMax ) : break
                if( region.domainMin != reconstructed.domainMax ) :
                    newregion = angularModule.XYs2d( axes = region.axes, interpolation=region.interpolation,
                            interpolationQualifier=region.interpolationQualifier )
                    newregion.append( region.getAtEnergy( reconstructed.domainMax ) )
                    for val in region:
                        if( val.outerDomainValue <= reconstructed.domainMax ) : continue
                        newregion.append( val )
                    merged.append( newregion )
                    idx += 1
                for region in original[idx:]:
                    merged.append( region )

            newForm = angularModule.TwoBody( label = newStyle.label,
                    productFrame = product.distribution.evaluated.productFrame, angularSubform = merged )
            if overwrite and newStyle.label in product.distribution:
                product.distribution.remove( newStyle.label )
            product.distribution.add( newForm )
        if newStyle.label not in self.styles: self.styles.add( newStyle )

    def removeStyle( self, style ) :
        """
        Remove the given style everywhere it appears in this reactionSuite
        :param style: may be a style instance or its name
        :return:
        """

        if isinstance(style, stylesModule.Style): style = style.label

        def removeStyleFromComponent( styleName, component ):
            if styleName in component:
                component.remove( styleName )

        for reaction in self :
            # currently checks in cross section and deposition data, should also check multiplicity and distributions
            if hasattr(reaction, 'crossSection'):
                removeStyleFromComponent( style, reaction.crossSection )

            if not hasattr(reaction, 'outputChannel'): continue
            for product in reaction.outputChannel:
                removeStyleFromComponent( style, product.averageProductEnergy )
                removeStyleFromComponent( style, product.averageProductMomentum )
                if product.outputChannel is None: continue
                for dproduct in product.outputChannel:
                    removeStyleFromComponent( style, dproduct.averageProductEnergy )
                    removeStyleFromComponent( style, dproduct.averageProductMomentum )
        self.styles.remove( style )

    def saveAllToFile(self, fileName, hybrid=False, **kwargs):
        """
        Writes *self* and its covariance data to files. The name of the covariance file will be the same as the *fileName* but
        it will be put into a "Covariance" sub-directory inside the directory where *self* is written and "-covar" will be
        added before the extension. If hybrid is *True*, self will be written out in default hybrid mode with the HDF5 data 
        written into a file with the same name as *self* but different extension and in a sub-directory named "HDF5".

        :param fileName:        The name of the file to write *self* to.
        :param hybrid:          If *True* data are written in hybrid XML/HDF5 files.
        :param kwargs:          Additional key-word arguments that are passed to internal calls.

        :return:                The list of paths for all covariance files written.
        """

        from . import externalFile as externalFileModule
        covarianceDir = pathlib.Path(kwargs.get('covarianceDir', 'Covariances'))
        if covarianceDir.is_absolute():
            raise Exception('Absolute covariance directory path not supported.')

        fileName = pathlib.Path(fileName)
        covariancePathRelative = covarianceDir / (fileName.stem + '-covar' + fileName.suffix)
        covariancePath = fileName.parent / covariancePathRelative                       # Path relative to self's path.
        selfsPathInCovarianceFile = pathlib.Path(os.path.relpath(fileName.parent, covariancePath.parent)) / fileName.name

        covariances = self.loadCovariances()
        originalPaths = {}
        covariancePaths = []
        if len(covariances) > 0:
            if self.sourcePath is None:                 # Self was generated directly, not parsed from a file. Add externalFiles before saving.
                covariances[0].externalFiles.add(externalFileModule.ExternalFile("reactions", path=selfsPathInCovarianceFile))
                covariances[0].saveToFile(covariancePath)
                covariancePaths.append(covariancePath)

                sha1sum = checksumsModule.Sha1sum.from_file(covariancePath)
                self.externalFiles.add(externalFileModule.ExternalFile("covariances", covariancePath, checksum=sha1sum))
            else:                                       # Locate and update the ExternalFiles pointing between the two files.
                if len(covariances) > 1:
                    raise Exception('Currently, only one covariance file is supported.')
                covarianceSuite = covariances[0]
                if covarianceSuite.domainUnit != self.domainUnit:
                    covarianceSuite.convertUnits({covarianceSuite.domainUnit: self.domainUnit})

                selfsExternalFile4Covariance = None                                     # Self's ExternalFile that points to the covariance file.
                for externalFile in self.externalFiles:
                    if externalFile.label == 'covariances':
                        selfsExternalFile4Covariance = externalFile
                if selfsExternalFile4Covariance is None:
                    raise Exception('Could not find externalFile for convarince in self.')

                if not pathlib.Path(selfsExternalFile4Covariance.path).is_absolute():   # Only write if path to covariance file is not absolute.
                                                                                        # Probably should have option to write even if absolute.
                    externalFile4ReactionSuite = None                                   # Covariance's ExternalFile that points to the reactionSuite file.
                    for externalFile in covarianceSuite.externalFiles:
                        if externalFile.label == 'reactions':
                            externalFile4ReactionSuite = externalFile
                    if externalFile4ReactionSuite is None:
                        raise Exception('Could not find externalFile for self in convarince.')

                    originalPaths[externalFile4ReactionSuite] = externalFile4ReactionSuite.path
                    externalFile4ReactionSuite.path = selfsPathInCovarianceFile
                    covarianceSuite.saveToFile(covariancePath, **kwargs)
                    covariancePaths.append(covariancePath)

                    originalPaths[selfsExternalFile4Covariance] = selfsExternalFile4Covariance.path
                    oldPath = selfsExternalFile4Covariance.path
                    selfsExternalFile4Covariance.path = covariancePath
                    selfsExternalFile4Covariance.path = str(covariancePathRelative)
                    if oldPath != covariancePath or selfsExternalFile4Covariance.checksum is None:
                        oldSourcePath = self.__sourcePath                                       # Need to revert back before returning.
                        self.__sourcePath = fileName                                            # Needed by updateChecksum() below.
                        selfsExternalFile4Covariance.updateChecksum()
                        self.__sourcePath = oldSourcePath


        if hybrid:
            self.saveToHybrid(fileName, **kwargs)
        else:
            self.saveToFile(fileName, **kwargs)

        for externalFile in originalPaths:
            externalFile.path = originalPaths[externalFile]

        return covariancePaths

    def saveToHybrid( self, xmlName, hdfName=None, hdfSubDir='HDF5', minLength=3, flatten=True, compress=False, **kwargs ):
        """
        Save the reactionSuite to a hybrid layout, with the data hierarchy in xml
        but most actual data saved in an associated HDF file

        :param xmlName: file name to save.
        :param hdfName: if None, will be the same as xmlName with extension changed to .h5
        :param hdfSubDir: if not None, the sub-directoty to add to path hdfName.
        :param minLength: minimum number of values before moving a dataset from XML to HDF
            Need to determine best default setting
        :param flatten: if True, GNDS datasets are concatenated into flattened HDF5 datasets
        :param compress: enable gzip + shuffle compression for HDF5 datasets
        :param kwargs:
        :return:
        """

        from . import externalFile as externalFileModule

        if not HDF5_present:
            raise ImportError("Missing module 'h5py' must be installed to support generating HDF5 files.")

        xmlName = pathlib.Path(xmlName)

        if hdfName is None:
            hdfName = xmlName
            if hdfName.suffix != '.xml':
                hdfName = hdfName.with_suffix(hdfName.suffix + '.xml')
            hdfName = hdfName.with_suffix('.h5')
        else:
            hdfName = pathlib.Path(hdfName)
        if hdfSubDir is not None:
            hdfName = hdfName.parent / hdfSubDir / pathlib.Path(hdfName).name

        if not hdfName.parent.exists():
            hdfName.parent.mkdir(parents=True)
        with h5py.File(hdfName, "w") as h5:
            HDF_opts = {
                'h5file': h5,
                'index': 0,     # assign unique name to each dataset
                'minLength': minLength,
                'flatten': flatten,
                'iData': [],
                'dData': [],
                'compression': {}
            }
            if compress:
                HDF_opts['compression'] = {'compression':'gzip', 'shuffle':'true', 'chunks':True}

            checksumAlgorithm = checksumsModule.Sha1sum
            relHdfName = hdfName
            try:
                relHdfName = pathlib.Path(hdfName).relative_to(xmlName.parent)
            except:
                if relHdfName == relHdfName.resolve():
                    relHdfName = relHdfName.name
            # add filler external file, actual checksum computed below
            self.externalFiles.add(externalFileModule.ExternalFile("HDF", str(relHdfName), checksum = "deadbeef", algorithm=checksumAlgorithm.algorithm))

            xmlString = self.toXML_strList( HDF_opts = HDF_opts, **kwargs )

            if len(HDF_opts['iData']) > 0:
                iData = numpy.array( HDF_opts['iData'], dtype=numpy.int32 )
                h5.create_dataset('iData', data=iData, **HDF_opts['compression'])
            if len(HDF_opts['dData']) > 0:
                dData = numpy.array( HDF_opts['dData'], dtype=numpy.float64 )
                h5.create_dataset('dData', data=dData, **HDF_opts['compression'])

        sha1sum = checksumAlgorithm.from_file(hdfName)
        for idx, line in enumerate(xmlString):
            if 'externalFile' in line and 'checksum="deadbeef"' in line:
                xmlString[idx] = line.replace("deadbeef", sha1sum)
                break

        if not xmlName.parent.exists():
            xmlName.parent.mkdir(parents=True)
        with open( xmlName, "w" ) as fout :
            fout.write( '<?xml version="1.0" encoding="UTF-8"?>\n')
            fout.write( '\n'.join( xmlString ) )
            fout.write( '\n' )

        self.externalFiles.pop("HDF")

    def thermalNeutronScatteringLawTemperatures( self ) :

        temperatures = {}
        for reaction in self.reactions : reaction.thermalNeutronScatteringLawTemperatures( temperatures )
        return( temperatures )

    def toXML_strList(self, indent="", **kwargs):
        """
        Returns a list of GNDS/XML strings representing self. To see all child suites created, including empty ones, 
        call with argument 'showEmptySuite=True'.
        """

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        if formatVersion in (GNDS_formatVersionModule.version_2_0_LLNL_3, GNDS_formatVersionModule.version_2_0_LLNL_4):
            print('INFO: converting GNDS format from "%s" to "%s".' % (formatVersion, GNDS_formatVersionModule.version_2_0))
            formatVersion = GNDS_formatVersionModule.version_2_0
        kwargs['formatVersion'] = formatVersion
        if formatVersion not in GNDS_formatVersionModule.allowed:
            raise Exception("Unsupported GNDS structure '%s'!" % str(formatVersion))

        interaction = self.interaction
        if interaction == enumsModule.Interaction.TNSL:
            interaction = enumsModule.Interaction.getTNSL_interaction(formatVersion)
        xmlString = ['%s<%s projectile="%s" target="%s" evaluation="%s" format="%s" projectileFrame="%s" interaction="%s">'
            % (indent, self.moniker, self.projectile, self.target, self.evaluation, formatVersion, self.projectileFrame, interaction)]

        xmlString += self.externalFiles.toXML_strList(indent2, **kwargs)
        xmlString += self.styles.toXML_strList(indent2, **kwargs)
        xmlString += self.PoPs.toXML_strList(indent2, **kwargs)

        if self.resonances is not None:
            xmlString += self.resonances.toXML_strList(indent2, **kwargs)

        xmlString += self.reactions.toXML_strList(indent2, **kwargs)
        xmlString += self.orphanProducts.toXML_strList(indent2, **kwargs)
        xmlString += self.sums.toXML_strList(indent2, **kwargs)
        xmlString += self.fissionComponents.toXML_strList(indent2, **kwargs)
        xmlString += self.productions.toXML_strList(indent2, **kwargs)
        xmlString += self.incompleteReactions.toXML_strList(indent2, **kwargs)
        xmlString += self.applicationData.toXML_strList(indent2, **kwargs)

        xmlString.append('%s</%s>' % (indent, self.moniker))
        return xmlString

    def toString( self, indent = '' ) :
        """Returns a string representation of an reactionSuite."""

        s = indent + self.inputParticlesToReactionString( suffix = ' --> ' )
        indent2 = '\n' + len( s ) * ' '
        indent3 = ''
        for reaction in self.reactions :
            s += reaction.toString( indent = indent3 )
            indent3 = indent2
        return( s )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Translates a <reactionSuite> xml node into a reactionSuite instance. Users should use the 'readXML_file' function instead."""

        xPath = ['reactionSuite']   # Keep track of location in the tree, in case errors come up.

        verbosity = kwargs.get('verbosity', 1)
        numberOfBrokenLinksToPrint = kwargs.get('numberOfBrokenLinksToPrint', 5)

        try :
            sourcePath = kwargs.get('sourcePath', './')     # This (i.e., './') may not be the correct default.
            if not isinstance(sourcePath, str):
                sourcePath = sourcePath.name
            path = pathlib.Path(sourcePath).resolve().parent

            formatVersion = node.get( 'format' )
            kwargs['formatVersion'] = formatVersion

            projectile = node.get( 'projectile' )
            target = node.get( 'target' )
            evaluation = node.get( 'evaluation' )
            projectileFrame = node.get( 'projectileFrame' )

            interaction = node.get('interaction')
            if interaction == enumsModule.Interaction.legacyTNSL or interaction == enumsModule.Interaction.TNSL:
                interaction = interaction = enumsModule.Interaction.TNSL
                guessedInteraction = False
            else:
                guessedInteraction, interaction = guessInteractionModule.guessInteraction(interaction, projectile, target)


            rs = cls(projectile, target, evaluation, formatVersion=formatVersion, projectileFrame=projectileFrame, interaction=interaction, sourcePath=sourcePath)

            if formatVersion == GNDS_formatVersionModule.version_1_10:
                rs.guessedInteraction = guessedInteraction
            elif guessedInteraction:
                raise Exception('Function guessInteraction had to guess at the interaction attribute for a non GNDS %s file.' % formatVersion)

            linkData = { 'reactionSuite' : rs, 'unresolvedLinks' : [], 'format' : formatVersion }

            externalFilesNode = node.find(suitesModule.ExternalFiles.moniker)
            if externalFilesNode is not None: rs.externalFiles.parseNode(externalFilesNode, xPath, linkData, **kwargs)

            if "HDF" in rs.externalFiles:   # FIXME hard-coded label
                if not HDF5_present:
                    raise ImportError("Missing module 'h5py' must be installed for HDF5 support.")
                HDF_externalFile = rs.externalFiles.pop('HDF')
                h5FilePath = path / HDF_externalFile.path
                if not path.exists():
                    raise IOError("External HDF file %s not found" % h5FilePath)

                h5File = h5py.File( h5FilePath, 'r', libver='latest' )
                linkData["HDF"] = { "h5File" : h5File }

            resonancesNode = node.find(resonancesModule.Resonances.moniker)
            if resonancesNode is not None:
                rs.resonances = resonancesModule.Resonances.parseNodeUsingClass(resonancesNode, xPath, linkData, **kwargs)

            kwargs['membersToSkip'] = [rs.externalFiles.moniker, resonancesModule.Resonances.moniker]

            childNodesNotParse, membersNotFoundInNode = rs.parseAncestryMembers(node, xPath, linkData, **kwargs)
            childNodesNotParse.pop(rs.externalFiles.moniker, None)
            childNodesNotParse.pop(resonancesModule.Resonances.moniker, None)

            oldDocumentationElement = None
            oldDocumentationElement_3 = None

            if formatVersion == GNDS_formatVersionModule.version_1_10: oldDocumentationElement = childNodesNotParse.pop('documentations', None)
            if formatVersion == GNDS_formatVersionModule.version_2_0_LLNL_3: oldDocumentationElement_3 = childNodesNotParse.pop('documentation', None)

            if oldDocumentationElement is not None:
                if len(oldDocumentationElement) > 1:
                    print( 'Too many children of old documentations node: number of children = %s.' % len( oldDocumentationElement ) )
                else :
                    oldDocumentationElement = oldDocumentationElement[0]
                    text = oldDocumentationElement.text.lstrip( '\n' )
                    if len( text.split( '\n' ) ) < 2:
                        rs.styles.getEvaluatedStyle().documentation.body.body = text
                    else:
                        rs.styles.getEvaluatedStyle().documentation.endfCompatible.body = text
            elif oldDocumentationElement_3:                         # Support some old 2.0.LLNL_3 files of YAHFC to GNDS for now. Hopefully, only temporary.
                rs.styles[0].documentation.parseNode(oldDocumentationElement_3, xPath, linkData, **kwargs)

            if len(childNodesNotParse) > 0: print("Warning: encountered unexpected child nodes '%s' in reactionSuite!" % ', '.join(list(childNodesNotParse.keys())))

        except Exception :
            print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
            raise

        # Fix links.
        unresolvedLinksCounter = 0
        for quant in linkData['unresolvedLinks'] :
            try :
                if quant.path.startswith( '/' + ReactionSuite.moniker ) :
                    quant.link = quant.follow( rs )
                elif quant.root:
                    rs._externalLinks.append( quant )
                else:   # relative link
                    quant.link = quant.follow( quant )
            except ancestryModule.XPathNotFound :
                if verbosity > 0:
                    if unresolvedLinksCounter < numberOfBrokenLinksToPrint: print( '    Cannot resolve link "%s".' % quant )
                    unresolvedLinksCounter += 1
            except :
                raise

        if( unresolvedLinksCounter > 0 ) : print( 'WARNING: %s unresolved links found while parsing file "%s".' % ( unresolvedLinksCounter, rs.sourcePath ) )

        return rs

    @staticmethod
    def read(fileName, **kwargs):
        """
        Reads in the file name *fileName* and returns a **ReactionSuite** instance.
        """

        return ReactionSuite.readXML_file(fileName, **kwargs)

def read(fileName, **kwargs):
    """
    Reads in the file name *fileName* and returns a **ReactionSuite** instance.
    """

    return ReactionSuite.read(fileName, **kwargs)

def niceSortOfMTs(MTs, verbose = 0, logFile=sys.stderr):
    '''
    Sorts the list of ENDF MTs into a sensible, defined order.
    '''

    def removeGetIfPresent(MT, MTs):
        '''
        For internal use only.
        '''

        if MT not in MTs:
            return []

        MTs.remove(MT)

        return [MT]

    MTs = sorted(list(copy.deepcopy(MTs)))

    newMTs = removeGetIfPresent(  2, MTs)                       # Elastic reaction.

    for MT in range( 50, 92):                                   # (z,n) reactions.
        newMTs += removeGetIfPresent(MT, MTs)
    newMTs += removeGetIfPresent(  4, MTs)

    newMTs += removeGetIfPresent( 16, MTs)                      # (z,2n) reactions.
    for MT in range(875, 892):
        newMTs += removeGetIfPresent(MT, MTs)

    newMTs += removeGetIfPresent( 17, MTs)
    newMTs += removeGetIfPresent( 37, MTs)

    newMTs += removeGetIfPresent( 18, MTs)                      # Fission reactions.
    newMTs += removeGetIfPresent( 19, MTs)
    newMTs += removeGetIfPresent( 20, MTs)
    newMTs += removeGetIfPresent( 21, MTs)
    newMTs += removeGetIfPresent( 38, MTs)

    newMTs += removeGetIfPresent( 28, MTs)
    newMTs += removeGetIfPresent( 32, MTs)
    newMTs += removeGetIfPresent( 33, MTs)

    for MT in range(600, 650):                                  # (z,p) reactions.
        newMTs += removeGetIfPresent(MT, MTs)
    newMTs += removeGetIfPresent(103, MTs)

    for MT in range(650, 700):                                  # (z,d) reactions.
        newMTs += removeGetIfPresent(MT, MTs)
    newMTs += removeGetIfPresent(104, MTs)

    for MT in range(700, 750):                                  # (z,t) reactions.
        newMTs += removeGetIfPresent(MT, MTs)
    newMTs += removeGetIfPresent(105, MTs)

    for MT in range( 750, 800):                                 # (z,He3) reactions.
        newMTs += removeGetIfPresent(MT, MTs)
    newMTs += removeGetIfPresent(106, MTs)

    for MT in range( 800, 850):                                 # (z,a) reactions.
        newMTs += removeGetIfPresent(MT, MTs)
    newMTs += removeGetIfPresent(107, MTs)

    for MT in range( 900, 1000):                                # (z,g) reactions.
        newMTs += removeGetIfPresent(MT, MTs)
    newMTs += removeGetIfPresent(102, MTs)

    MT5 = removeGetIfPresent(  5, MTs)                          # (z,everything else).

    MTAtomics = []
    for MT in range(500, 573):
        MTAtomics += removeGetIfPresent(MT, MTs)

    skippingMTs = []
    for MT in [10, 27, 101, 151]:
        skippingMTs += removeGetIfPresent(MT, MTs)
    for MT in range(201, 600):
        skippingMTs += removeGetIfPresent(MT, MTs)
    for MT in range(850, 875):
        skippingMTs += removeGetIfPresent(MT, MTs)

    if verbose > 0 and len(skippingMTs) > 0:
        logFile.write('Skipping MTs = %s\n' % skippingMTs)

    return newMTs + MTs + MT5 + MTAtomics
