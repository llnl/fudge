# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import vector as vectorModule
from xData import matrix as matrixModule
from xData import productArray as productArrayModule
from xData.Documentation import documentation as documentationModule

from .. import enums as enumsModule
from .. import outputChannel as outputChannelModule 
from ..reactionData.doubleDifferentialCrossSection import doubleDifferentialCrossSection as doubleDifferentialCrossSectionModule
from ..reactionData import crossSection as crossSectionModule
from ..reactionData import availableEnergy as availableEnergyModule
from ..reactionData import availableMomentum as availableMomentumModule

class Base_reaction(ancestryModule.AncestryIO):
    """Base class for all types of reaction."""

    ancestryMembers = ('crossSection', 'outputChannel')
    keyName = 'label'

    def __init__(self, label, genre, ENDF_MT):

        ancestryModule.AncestryIO.__init__( self )
        self.__label = label

        self.__documentation = None
        self.ENDF_MT = int( ENDF_MT )

        self.__crossSection = crossSectionModule.Component( )
        self.__crossSection.setAncestor( self )

        self.__outputChannel = outputChannelModule.OutputChannel(genre)
        self.__outputChannel.setAncestor( self )

    def __str__( self ) :

        return( self.label )

    @property
    def documentation(self):
        return self.__documentation

    @documentation.setter
    def documentation(self, _documentation):
        assert isinstance(_documentation, (documentationModule.Documentation, type(None)))
        if _documentation is not None:
            _documentation.setAncestor(self)
        self.__documentation = _documentation

    @property
    def crossSection( self ) :

        return( self.__crossSection )

    @property
    def outputChannel( self ) :

        return( self.__outputChannel )

    @property
    def label( self ) :
        """Returns the reaction's label."""

        if( self.__label is None ) : self.updateLabel( )
        return( self.__label )

    def updateLabel( self ) :
        """Sets the reaction's label from outputChannel products."""

        label = self.__outputChannel.toString( MT = self.ENDF_MT )
        if hasattr(self, 'fissionGenre'):
            if self.fissionGenre is not enumsModule.FissionGenre.none:
                if '[' in label: raise Exception('Label already contains a process: "%s"' % label)
                if len(label) == 0:
                    label = '%s fission' % self.fissionGenre
                else:
                    label += ' [%s fission]' % self.fissionGenre
        self.__label = label

    def findLinks(self, links):

        for ancestryMember in self.ancestryMembers :
            if hasattr(self, 'doubleDifferentialCrossSection'):
                if ancestryMember in ('doubleDifferentialCrossSection',): continue
            getattr(self, ancestryMember).findLinks(links)

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ):

        return( False )

    def isCoulombReaction(self):
        '''Returns **True** if *self* is a Coulomb reaction and **False** otherwise.'''

        if len(self.crossSection) > 0:
            return isinstance(self.crossSection[0], crossSectionModule.CoulombPlusNuclearElastic)

        return False

    def isFission( self ):

        if not hasattr(self, 'fissionGenre'):
            return False

        return self.fissionGenre != enumsModule.FissionGenre.none

    def isThermalNeutronScatteringLaw(self):

        if hasattr(self, 'doubleDifferentialCrossSection'):
            for form in self.doubleDifferentialCrossSection:
                if form.isThermalNeutronScatteringLaw(): return True

        return False

    def isPairProduction(self):
        """
        Return a boolean indicating whether the reaction is pair production. 
        """

        return self.ENDF_MT in [515, 516, 517]

    def check(self, info):
        """
        This method is usually not called directly. Use reactionSuite.check() instead.

        Checks cross section and outputChannel (and differential cross sections if available). Checks include:
          do Z/A balance? do Q-value and thresholds agree?
          Does cross section domain agree with each product distribution/multiplicity domain?
          Does energy balance?

        :param info: dict
        :return list of warnings
        """

        from fudge import warning
        from . import production as productionModule, orphanProduct as orphanProductModule, \
            incompleteReaction as incompleteReactionModule, fissionComponent as fissionComponentModule
        from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic.CoulombPlusNuclearElastic \
            import CoulombDepositionNotSupported
        warnings = []

        reactionSuite = self.rootAncestor
        energyUnit = reactionSuite.domainUnit

        def particleZA(particleID):

            particle = reactionSuite.PoPs.final(particleID)
            return chemicalElementMiscPoPsModule.ZA(particle)


        info['crossSectionDomain'] = self.crossSection.domainMin, self.crossSection.domainMax
        info['isTwoBody'] = self.__outputChannel.genre == enumsModule.Genre.twoBody
        try:
            info['Q'] = self.getQ(energyUnit, final=False)
        except ValueError:
            # Q is energy-dependent. Compute Q at threshold
            info['Q'] = self.outputChannel.Q.evaluated.evaluate(self.crossSection.domainMin)

        cpcount = sum([(particleZA(prod.pid) // 1000) > 0 for prod in self.__outputChannel])
        info['CoulombOutputChannel'] = cpcount > 1
        info['ContinuumOutputChannel'] = self.outputChannel.process == outputChannelModule.Processes.continuum

        if hasattr(self, 'doubleDifferentialCrossSection'):
            differentialCrossSectionWarnings = self.doubleDifferentialCrossSection.check(info)
            if differentialCrossSectionWarnings:
                warnings.append(warning.Context('doubleDifferentialCrossSection:', differentialCrossSectionWarnings))

        crossSectionWarnings = self.crossSection.check(info)
        if crossSectionWarnings:
            warnings.append(warning.Context("Cross section:", crossSectionWarnings))

        del info['CoulombOutputChannel']
        del info['ContinuumOutputChannel']

        if info['crossSectionOnly']:
            return warnings             # otherwise continue to check outputs

        # compare calculated and listed Q-values:
        if not isinstance(self, (productionModule.Production, fissionComponentModule.FissionComponent,
                                 incompleteReactionModule.IncompleteReaction, orphanProductModule.OrphanProduct)):
            try:
                Q = info['Q']
                Qcalc = info['availableEnergy']
                if Qcalc is None: raise ValueError  # caught below. Skips Q-value check for elemental targets
                for prod in self.__outputChannel.products:
                    try:
                        productMass = prod.getMass(f'{energyUnit}/c**2')
                    except Exception:
                        warnings.append(warning.UnknownMass(prod.pid))
                        raise ValueError("Unknown mass")

                    try:
                        productMultiplicity = prod.multiplicity.getConstant()
                    except Exception:   # multiplicity is not constant
                        if prod.pid == IDsPoPsModule.photon: continue
                        raise ValueError("Non-constant multiplicity")

                    Qcalc -= productMass * productMultiplicity

                if abs(Q-Qcalc) > PQUModule.PQU(info['dQ']).getValueAs(energyUnit):
                    if self.__outputChannel.process != outputChannelModule.Processes.continuum:
                        warnings.append(
                            warning.Q_mismatch(PQUModule.PQU(Qcalc, energyUnit), PQUModule.PQU(Q, energyUnit), self))
            except ValueError:
                pass    # this test only works if multiplicity and Q are both constant for all non-gamma products

        if not (self.__outputChannel.genre == enumsModule.Genre.sumOfRemainingOutputChannels or self.isFission() or
                isinstance(self, (productionModule.Production, orphanProductModule.OrphanProduct,
                                  incompleteReactionModule.IncompleteReaction))):
            # check that ZA balances:
            ZAsum = 0
            for product in self.__outputChannel.products:
                if product.pid == IDsPoPsModule.photon: continue
                try:
                    mult = product.multiplicity.getConstant()
                except Exception as ex:
                    warnings.append(warning.NonConstantMultiplicity(self))
                    # likely also causes ZAbalanceWarning
                else:
                    ZAsum += particleZA(product.pid) * mult
            if info['elementalTarget']:
                ZAsum = (ZAsum // 1000) * 1000
            if ZAsum != info['compoundZA']:
                warnings.append( warning.ZAbalanceWarning( self ) )

        if not isinstance(self, productionModule.Production):
            for product in self.__outputChannel.products:
                productWarnings = product.check(info)
                if productWarnings:
                    warnings.append(warning.Context("Product: %s" % product.label, productWarnings))

        fissionFragmentWarnings = self.__outputChannel.fissionFragmentData.check(info)
        if fissionFragmentWarnings:
            warnings.append(warning.Context("Fission fragment info:", fissionFragmentWarnings))

        del info['crossSectionDomain']
        del info['isTwoBody']
        del info['Q']

        if (info['checkEnergyBalance']
                and not isinstance(self, productionModule.Production)
                and self.__outputChannel.genre is not enumsModule.Genre.sumOfRemainingOutputChannels):
            # Calculate energy deposition data for all products, and for decay products.
            # Then recursively check the product list for energy balance at each step of the reaction/decay
            # At each step, if we have energy deposition for *every* product in the list, we can rigorously check
            # energy balance. Otherwise, can only check that we don't exceed available energy.
            try:
                self.calculateAverageProductData(style=info['averageProductDataStyle'], **info['averageProductDataArgs'])
            except CoulombDepositionNotSupported:
                warnings.append(warning.SkippedCoulombElasticEnergyDeposition(self))
            except Exception as e:
                warnings.append(warning.EnergyDepositionExceptionRaised(str(e), self))
                if info['failOnException']: raise
                return warnings

            def checkProductsForEnergyBalance(products, Qs, fission=False, decay=False):
                # sample usage: for the reaction n + F19 -> n + (F19_c -> F19 + gamma), this function
                # should be called twice, to test energy balance at each step of the reaction.
                # First call: products = [n, F19_c] and Qs = [Q_reaction],
                # second call: products = [n, F19, gamma] and Qs = [Q_reaction, Q_decay].
                edepWarnings = []

                averageProductDataLabel = info['averageProductDataStyle'].label
                energyDep = []
                for prod in products:
                    if averageProductDataLabel in prod.averageProductEnergy:
                        energyDep.append([prod.label, prod.averageProductEnergy[averageProductDataLabel]])
                if energyDep:
                    totalEDep = energyDep[0][1]
                    for idx in range(1, len(energyDep)):
                        if (totalEDep.domainMin != energyDep[idx][1].domainMin
                                or totalEDep.domainMax != energyDep[idx][1].domainMax):
                            upperEps = 0
                            if totalEDep.domainMax != energyDep[idx][1].domainMax:
                                upperEps = 1e-8
                            try:
                                totalEDep, energyDep[idx][1] = totalEDep.mutualify(
                                        1e-8, upperEps, 0, energyDep[idx][1], 1e-8, upperEps, 0)
                            except Exception as e:
                                edepWarnings.append(warning.EnergyDepositionExceptionRaised(str(e), self))
                                if info['failOnException']: raise
                                return warnings
                        totalEDep = totalEDep + energyDep[idx][1]
                else: totalEDep = []
                Qsum = sum(Qs)

                def energyDepositedPerProduct(energyDep, ein):
                    """ if energy doesn't balance, determine how much is deposited in each product """
                    result = []
                    availableEnergy = ein + Qsum
                    if availableEnergy == 0: availableEnergy = sum([edep.evaluate(ein) for p, edep in energyDep])
                    for prod, edep in energyDep:
                        edep = edep.evaluate(ein)
                        if edep is None: result.append((prod, 0))
                        else: result.append((prod, 100.0 * edep/availableEnergy))
                    return sorted(result, key=lambda foo: foo[1])[::-1]

                # now we have total energy deposition for all particles, use to check energy balance.
                # a few special cases to consider:
                if fission:
                    # fission products aren't listed (so far, anyway), so about 85% of available energy should be missing:
                    for i, (ein, edep) in enumerate(totalEDep):
                        if edep > abs((ein + Qsum) * info['fissionEnergyBalanceLimit']):
                            edepWarnings.append(warning.FissionEnergyImbalance(PQUModule.PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self))
                elif len(products) == len(energyDep):
                    # have full energy dep. data for all products, so we can rigorously check energy balance:
                    for i, (ein, edep) in enumerate(totalEDep):
                        if (abs(edep - (ein+Qsum)) > abs((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and abs(edep - (ein+Qsum)) > PQUModule.PQU( info['dEnergyBalanceAbsolute'])
                                        .getValueAs(totalEDep.axes[0].unit)):
                            edepWarnings.append(warning.EnergyImbalance(PQUModule.PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self))
                else:
                    # missing some products, so just check that outgoing energy doesn't exceed incoming:
                    for i, (ein, edep) in enumerate(totalEDep):
                        if ((edep - (ein+Qsum)) > ((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and (edep - (ein+Qsum)) > PQUModule.PQU(info['dEnergyBalanceAbsolute'])
                                        .getValueAs(totalEDep.axes[0].unit)):
                            edepWarnings.append(
                                warning.EnergyImbalance(PQUModule.PQU(ein, totalEDep.axes[0].unit),
                                                        i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self))

                if edepWarnings:
                    productList = [prod.pid for prod in products]
                    photons = productList.count("photon")
                    if photons > 2:
                        # simplify context message if many photons are produced
                        photon1 = productList.index('photon')
                        productList[productList.index("photon")] = f"{photons}*photon"
                        productList = [p for p in productList if p != "photon"]
                    context = "Energy balance"
                    if decay: context += " (after decay)"
                    context += " for products: " + ', '.join(productList)
                    warnings.append(warning.Context(context, edepWarnings))

                # now recursively check decay products, if any:
                for pidx, currentProd in enumerate(products):
                    if currentProd.outputChannel is not None:
                        checkProductsForEnergyBalance(
                                products[:pidx] + [p for p in currentProd.outputChannel] + products[pidx+1:],
                                Qs + [currentProd.outputChannel.Q.getConstant()],
                                decay=True, fission=False        # FIXME what about spontaneous fission decay?
                                )
                # end of helper function checkProductsForEnergyBalance

            try:
                Q = self.__outputChannel.Q.getConstant()
                checkProductsForEnergyBalance(products=[p1 for p1 in self.__outputChannel], Qs=[Q],
                                              fission=self.isFission(), decay=False)
            except ValueError:
                pass    # FIXME this test currently disabled when Q non-constant

        return warnings

    @property
    def domainMin( self ) :

        return( self.crossSection.domainMin )

    @property
    def domainMax( self ) :

        return( self.crossSection.domainMax )

    @property
    def domainUnit( self ) :

        return( self.crossSection.domainUnit )

    def domainUnitConversionFactor( self, unitTo ) :

        return( self.crossSection.domainUnitConversionFactor( unitTo ) )

    def convertUnits(self, unitMap):
        """See documentation for reactionSuite.convertUnits."""

        from . import reaction as reactionModule

        if hasattr(self, 'doubleDifferentialCrossSection'):
            self.doubleDifferentialCrossSection.convertUnits(unitMap)
        self.__crossSection.convertUnits(unitMap)
        self.__outputChannel.convertUnits(unitMap)

        if( isinstance( self, reactionModule.Reaction ) ) :
            self.availableEnergy.convertUnits( unitMap )
            self.availableMomentum.convertUnits( unitMap )

    def amendForPatch(self, fromLabel, toLabel):

        from . import reaction as reactionModule

        if hasattr(self, 'doubleDifferentialCrossSection'):
            self.doubleDifferentialCrossSection.amendForPatch(fromLabel, toLabel)
        self.__crossSection.amendForPatch(fromLabel, toLabel)
        self.__outputChannel.amendForPatch(fromLabel, toLabel)

        if isinstance(self, reactionModule.Reaction):
            self.availableEnergy.amendForPatch(fromLabel, toLabel)
            self.availableMomentum.amendForPatch(fromLabel, toLabel)

    def diff( self, other, diffResults ) :

        self.__crossSection.diff( other.crossSection, diffResults )
        self.__outputChannel.diff( other.outputChannel, diffResults )

    def effectiveThreshold(self):
        '''
        Some reactions have a cross section with multiple zero values at their beginning. In such a case, the effective 
        threshold is the point before the first non-zero cross section value. Otherwise, it is the energy of the first point.
        '''

        return self.crossSection.effectiveThreshold()

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** method on the members *crossSection*, *doubleDifferentialCrossSection* and *outputChannel*.
        The *energyMin* argument is ignored. Instead, it is calculated from the minimum energy from the cross sections.
        """

        energyMin = self.__crossSection.domainMin
        energyMax = min(self.__crossSection.domainMax, energyMax)
        numberOfFixes  = self.__crossSection.fixDomains(labels, energyMin, energyMax)
        if hasattr(self, 'doubleDifferentialCrossSection'):
            numberOfFixes += self.doubleDifferentialCrossSection.fixDomains(labels, energyMin, energyMax)
        numberOfFixes += self.__outputChannel.fixDomains(labels, energyMin, energyMax)

        return numberOfFixes

    def heatCrossSection(self, temperature, EMin, lowerlimit=None, upperlimit=None, interpolationAccuracy=0.001, heatAllPoints=False,
            doNotThin=True, heatBelowThreshold=True, heatAllEDomain=True, setThresholdToZero=False, verbose=0):

        if hasattr(self, 'doubleDifferentialCrossSection'):
            if len(self.doubleDifferentialCrossSection) > 0:
                if verbose > 0: print('    Skipping doubleDifferentialCrossSection reaction.')
                return

        crossSection = self.crossSection.heat(temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin,
                heatBelowThreshold, heatAllEDomain, setThresholdToZero=setThresholdToZero, addToSuite=True)

        return crossSection

    def thresholdQAs(self, unit, final=True, pops=None) :

        return self.__outputChannel.thresholdQAs(unit, final=final, pops=pops)

    def getQ( self, unit, final = True ) :
        """Returns the Q-value for this reaction. Converted to float if possible, otherwise a string value is returned."""

        return( self.thresholdQAs( unit, final = final ) )

    def getReactionSuite( self ) :

        from .. import reactionSuite as reactionSuiteModule
        return( self.findClassInAncestry( reactionSuiteModule.ReactionSuite ) )

    def cullStyles( self, styleList ) :
        """ See documentation for reactionSuite.cullStyles. """

        from . import reaction as reactionModule

#        if hasattr(self, 'doubleDifferentialCrossSection'):
#            self.doubleDifferentialCrossSection.cullStyles( styleList )
        self.__crossSection.cullStyles( styleList )
        self.__outputChannel.cullStyles( styleList )
        if( isinstance( self, reactionModule.Reaction ) ) :
            self.availableEnergy.cullStyles( styleList )
            self.availableMomentum.cullStyles( styleList )

    def calculateAverageProductData(self, style, indent='', **kwargs):
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        verbosity = kwargs['verbosity']
        indent2 = indent + kwargs['incrementalIndent']

        if verbosity > 0:
            print('%s%s' % (indent, self.__outputChannel.toString(simpleString=True, MT=self.ENDF_MT)))

        kwargs['reaction'] = self
        kwargs['EMin'] = self.domainMin
        kwargs['EMax'] = self.domainMax
        self.__outputChannel.calculateAverageProductData(style, indent2, **kwargs)
        if hasattr(self, 'availableEnergy'):
            axes = availableEnergyModule.defaultAxes(self.domainUnit)
            QToPointwiseLinear = self.__outputChannel.QToPointwiseLinear(final=True)
            domainMin = QToPointwiseLinear.domainMin
            if self.isCoulombReaction():
                domainMin = min(self.domainMin, domainMin, PQUModule.PQU(1e-4, 'MeV').getValueAs(self.domainUnit))
                QToPointwiseLinear.setData([[domainMin, 0.0], [QToPointwiseLinear.domainMax, 0.0]])

            availableEnergy = availableEnergyModule.XYs1d(data=[[domainMin, QToPointwiseLinear.domainMin], [self.domainMax, self.domainMax]], axes=axes)
            availableEnergy += QToPointwiseLinear
            self.availableEnergy.add(availableEnergyModule.XYs1d(data=availableEnergy, axes=axes, label=style.label))

            massInE = kwargs['projectileMass']
            availableMomentum = availableMomentumModule.calculateMomentumPoints(style, massInE, domainMin, QToPointwiseLinear.domainMax, self.domainUnit)
            self.availableMomentum.add(availableMomentum)

    def listOfProducts(self, finalOnly=False, includeQualifier=True):
        '''
        Returns, as a python set, the list of PoPs ids for all products (i.e., outgoing particles) for *self*.

        :param finalOnly:           If `True`, only final product ids are returned, otherwise, all are returned..
        :param includeQualifier:    If `True`, particle qualifiers are include in ids, otherwise, they are stripped from ids.

        :return:                    A python set instance.
        '''

        products = set()
        if hasattr(self, 'doubleDifferentialCrossSection'):
            products.update(self.doubleDifferentialCrossSection.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))
        products.update(self.__outputChannel.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))

        return products

    def partialProductionIntegral( self, reaction_suite, productID, energyIn, energyOut = None, muOut = None, phiOut = None, 
            frame = xDataEnumsModule.Frame.lab, LegendreOrder = 0, **kwargs ) :

        if( isinstance( self.crossSection[0], crossSectionModule.CoulombPlusNuclearElastic ) ) :
            issueCounters = kwargs.get( 'issueCounters', None )
            if( issueCounters is not None ) : issueCounter = issueCounters.get( 'partialProductionIntegral:CoulombPlusNuclearElastic', 1, True )
            status = issueCounter.increment( )
            if( status ) : print( '    WARNING: partialProductionIntegral: skipping elastic Coulomb reaction %s.' % self )
            return( 0.0 )

        crossSection = self.crossSection.evaluate( energyIn )
        if( crossSection == 0.0 ) : return( 0.0 )
        _partialProductionIntegral = self.__outputChannel.partialProductionIntegral( reaction_suite, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )

        return( crossSection * _partialProductionIntegral )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        if( verbosity > 0 ) : print( '%s%s' % (indent, self.__outputChannel.toString( simpleString = True, MT = self.ENDF_MT ) ) )

        self.__outputChannel.processMC_cdf( style, tempInfo, indent2 )

    def processGriddedCrossSections( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', isPhotoAtomic = False ) :

        self.crossSection.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent, isPhotoAtomic = isPhotoAtomic )

    def processMultiGroup(self, style, tempInfo, indent):
        """ See documentation for reactionSuite.processMultiGroup. """

        from . import reaction as reactionModule

        tempInfo['workFile'].append( 'r%s' % tempInfo['reactionIndex'] )
        tempInfo['transferMatrixComment'] = tempInfo['reactionSuite'].inputParticlesToReactionString( suffix = " --> " ) + self.toString( )

        tempInfo['crossSectionModified'] = False
        reprocessingOpts = tempInfo['partialReprocessingOptions']
        if reprocessingOpts is not None:
            reprocessingOpts['otherStyleLabel'], = [
                s.heatedMultiGroup for s in reprocessingOpts['processedDataSource'].styles.temperatures()
                if s.temperature == style.temperature.value]

            if self in reprocessingOpts['crossSectionChanged']:
                tempInfo['crossSectionModified'] = True

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        if verbosity > 0: print('%s%s' % (indent, self.__outputChannel.toString(simpleString = True, MT = self.ENDF_MT)))

        tempInfo['reaction'] = self
        norm = tempInfo['groupedFlux']
        tempInfo['groupedFlux'] = None

        crossSection = style.findFormMatchingDerivedStyle( self.crossSection )
# BRB FIXME The next line is a kludge, see note on crossSection.resonancesWithBackground.processMultiGroup.
        if isinstance(crossSection, crossSectionModule.Reference):
            crossSection = crossSection.crossSection
        if isinstance(crossSection, crossSectionModule.ResonancesWithBackground):
            crossSection = crossSection.ancestor['recon']
        if not isinstance(crossSection, crossSectionModule.XYs1d):
            crossSection = crossSection.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
        tempInfo['crossSection'] = crossSection
        tempInfo['multiGroupCrossSection'] = self.crossSection.processMultiGroup(style, tempInfo, indent2)
        self.crossSection.remove(style.label)                             # Not normalized by tempInfo['groupedFlux'] so remove.

        tempInfo['groupedFlux'] = norm
        tempInfo['multiGroupCrossSectionNormed'] = self.crossSection.processMultiGroup(style, tempInfo, indent2)     # Normalized by tempInfo['groupedFlux'].

        if isinstance(self, reactionModule.Reaction):
            self.availableEnergy.processMultiGroup(style, tempInfo, indent2)
            self.availableMomentum.processMultiGroup(style, tempInfo, indent2)
        self.__outputChannel.processMultiGroup(style, tempInfo, indent2)

        if reprocessingOpts is not None:
            del reprocessingOpts['otherStyleLabel']
        del tempInfo['workFile'][-1]
        del tempInfo['crossSectionModified']

    def multiGroupCrossSection(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, cross section for the requested label in the reaction.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        crossSection = multiGroupSettings.formAsVector(self.crossSection, temperatureInfo)
        
        return crossSection

    def multiGroupQ(self, multiGroupSettings, temperatureInfo, includeFinalProducts):
        """
        Returns the multi-group, total Q for the requested label. This is a cross section weighted Q summed over all reactions.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param includeFinalProducts: Boolean value indicating whether to include the contriibution from the final fission product data.
        """
        
        return self.outputChannel.multiGroupQ(multiGroupSettings, temperatureInfo, includeFinalProducts)

    def multiGroupMultiplicity(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the multi-group, total multiplicity for the requested label for the requested product.
        
        This is a cross section weighted multiplicity.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        multiplicity = vectorModule.Vector()
        if self.isPairProduction():
            if productID == IDsPoPsModule.photon:
                multiplicity += 2*self.multiGroupCrossSection(multiGroupSettings, temperatureInfo)
            else:
                multiplicity += vectorModule.Vector()
        else:
            multiplicity += self.outputChannel.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, productID)

        return multiplicity

    def multiGroupAvailableEnergy(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, total available energy for the requested label. This is a cross section weighted available energy.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        from . import reaction as reactionModule

        if not isinstance(self, reactionModule.Reaction):
            return vectorModule.Vector()

        return  multiGroupSettings.formAsVector(self.availableEnergy, temperatureInfo)

    def multiGroupAverageEnergy(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the multi-group, total average energy for the requested label for the requested product.
        
        This is a cross section weighted average energy summed over all products for this reaction.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        averageEnergy = vectorModule.Vector()
        if self.isPairProduction():
            if productID == IDsPoPsModule.photon:
                elecronMass = PQUModule.PQU(1, 'me * c**2').getValueAs(self.domainUnit)
                averageEnergy += 2*elecronMass*self.multiGroupCrossSection(multiGroupSettings, temperatureInfo)

        else:
            averageEnergy += self.outputChannel.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, productID)

        return averageEnergy

    def multiGroupDepositionEnergy(self, multiGroupSettings, temperatureInfo, particleIDs):
        '''
        Returns the multi-group, deposition energy for the requested label.

        This is a cross section weighted deposition energy for the reaction. The deposition energy is calculated by subtracting
        the average energy from each transportable particle from the available energy. The list of transportable particles is specified
        via the list of particle specified in the *particleIDs* argument.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        '''

        from . import reaction as reactionModule

        if not isinstance(self, reactionModule.Reaction):
            return vectorModule.Vector()

        depositionEnergy = self.multiGroupAvailableEnergy(multiGroupSettings, temperatureInfo)
        for particleID in particleIDs:
            depositionEnergy -= self.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, particleID)

        return depositionEnergy

    def multiGroupAvailableMomentum(self, multiGroupSettings, temperatureInfo):
        """
        Returns the multi-group, total available momentum for the requested label. This is a cross section weighted available momentum.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        from . import reaction as reactionModule

        if not isinstance(self, reactionModule.Reaction):
            return vectorModule.Vector()

        return  multiGroupSettings.formAsVector(self.availableMomentum, temperatureInfo)

    def multiGroupAverageMomentum(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the multi-group, total average momentum for the requested label for the requested product.
        
        This is a cross section weighted average momentum.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        if self.isPairProduction():
            averageMomentum = vectorModule.Vector()
        else:
            averageMomentum = self.outputChannel.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, productID)

        return averageMomentum

    def multiGroupDepositionMomentum(self, multiGroupSettings, temperatureInfo, particleIDs):
        '''
        Returns the multi-group, total deposition momentum for the requested label.

        This is a cross section weighted deposition momentum for the reaction. The deposition momentum is calculated by
        subtracting the average momentum from each transportable particle from the available momentum. The list of transportable
        particles is specified via the list of particle specified in the *particleIDs* argument.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        '''

        from . import reaction as reactionModule

        if not isinstance(self, reactionModule.Reaction):
            return vectorModule.Vector()

        depositionMomentum = self.multiGroupAvailableMomentum(multiGroupSettings, temperatureInfo)
        for particleID in particleIDs:
            depositionMomentum -= self.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, particleID)

        return depositionMomentum 

    def multiGroupGain(self, multiGroupSettings, temperatureInfo, productID, projectileID):
        '''
        Returns multi-group, gain for the requested particle and label.

        This is a cross section weighted gain summed over all reactions. If productID and projectileID are the same, then the multi-group 
        cross section is subtracted for the returned value to indicate that the projectileID as been absorted.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param producID: The particle PoPs' id for the whose gain is to be calculated.
        :param projectileID: The particle PoPs' id for the projectile.
        '''

        particleGain = self.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, productID)
        if productID == projectileID:
            particleGain -= self.multiGroupCrossSection(multiGroupSettings, temperatureInfo)

        return particleGain

    def multiGroupProductMatrix(self, multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder):
        """
        Returns the multi-group, product matrix for the requested label for the requested productID for the requested Legendre order.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        :param legendreOrder: Requested Legendre order.
        """

        if self.isPairProduction():
            if productID == IDsPoPsModule.photon and legendreOrder == 0:
                productCrossSection = 2 * self.multiGroupCrossSection(multiGroupSettings, temperatureInfo)
                photonParticle = particleIDs[IDsPoPsModule.photon]
                electronMass = PQUModule.PQU(1, 'me * c**2').getValueAs(self.domainUnit)
                multiGroupIndexFromEnergy = photonParticle.multiGroupIndexFromEnergy(electronMass, True)
                productMatrix = matrixModule.Matrix(productCrossSection.size, productCrossSection.size)
                for index, productionSection in enumerate(productCrossSection):
                    productMatrix[index, multiGroupIndexFromEnergy] = productionSection

            else:
                productMatrix = matrixModule.Matrix()
            
        else:
            productMatrix = self.outputChannel.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder)

        return productMatrix

    def multiGroupFissionMatrix(self, multiGroupSettings, temperatureInfo, particleIDs, legendreOrder):
        """
        Returns the multi-group, fission neutron transfer matrix for the requested label for the requested Legendre order.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs:  The list of particles to be transported.
        :param legendreOrder: Requested Legendre order.
        """

        if self.isFission():
            fissionMatrix = self.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, IDsPoPsModule.neutron, legendreOrder)
    
        else:
            fissionMatrix = matrixModule.Matrix()

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

        if self.isPairProduction():
            if productID == IDsPoPsModule.photon:
                matrix = self.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, productID, 0)
                productArray = productArrayModule.ProductArray(matrix.matrix[:,:,None])
        else:
            productArray = self.outputChannel.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

        return productArray

    def removeStyles( self, styleLabels ) :

        from . import reaction as reactionModule

        if hasattr(self, 'doubleDifferentialCrossSection'):
            self.doubleDifferentialCrossSection.removeStyles( styleLabels )
        self.__crossSection.removeStyles( styleLabels )
        self.__outputChannel.removeStyles( styleLabels )
        if( isinstance( self, reactionModule.Reaction ) ) :
            self.availableEnergy.removeStyles( styleLabels )
            self.availableMomentum.removeStyles( styleLabels )

    def toString( self, indent = '' ) :

        return( self.__outputChannel.toString( indent = indent, MT = self.ENDF_MT ) )

    def toXML_strList( self, indent = '', **kwargs ) :

        from . import reaction as reactionModule

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        attributeString = ""
        attributeString += ' ENDF_MT="%s"' % self.ENDF_MT
        if hasattr(self, 'fissionGenre'):
            if self.fissionGenre != enumsModule.FissionGenre.none:
                attributeString += ' fissionGenre="%s"' % self.fissionGenre

        xmlString = [ '%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
        xmlString[-1] += attributeString + '>'

        if self.documentation is not None:
            xmlString += self.documentation.toXML_strList( indent2, **kwargs )

        if hasattr(self, 'doubleDifferentialCrossSection'):
            xmlString += self.doubleDifferentialCrossSection.toXML_strList( indent2, **kwargs )
        xmlString += self.__crossSection.toXML_strList( indent2, **kwargs )
        xmlString += self.__outputChannel.toXML_strList( indent2, **kwargs )

        if( isinstance( self, reactionModule.Reaction ) ) :
            xmlString += self.availableEnergy.toXML_strList( indent2, **kwargs )
            xmlString += self.availableMomentum.toXML_strList( indent2, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append('%s[@label="%s"]' % (node.tag, node.get('label')))

        self.__label = node.get('label')
        self.ENDF_MT = int(node.get('ENDF_MT'))
        if hasattr(self, 'fissionGenre'):
            self.fissionGenre = node.get('fissionGenre', enumsModule.FissionGenre.none)

        childNodesNotParse, membersNotFoundInNode = self.parseAncestryMembers(node, xPath, linkData, **kwargs)
        if len(childNodesNotParse) > 0: raise Exception("Encountered unexpected child nodes '%s' in %s!" % (self.moniker, ', '.join(list(childNodesNotParse.keys()))))

        if hasattr(self, 'fissionGenre'):
            if self.fissionGenre is enumsModule.FissionGenre.none:      # See note in outputChannel.parseNode.
                self.fissionGenre = self.outputChannel.fissionGenre
                if self.fissionGenre is enumsModule.FissionGenre.none:  # Additional kludge needed for some legacy files.
                    self.fissionGenre = {18: enumsModule.FissionGenre.total, 19: enumsModule.FissionGenre.firstChance,
                                         20: enumsModule.FissionGenre.secondChance, 21: enumsModule.FissionGenre.thirdChance,
                                         38: enumsModule.FissionGenre.fourthChance}.get(self.ENDF_MT, enumsModule.FissionGenre.none)

        xPath.pop( )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append('%s[@label="%s"]' % (node.tag, node.get('label')))

        if issubclass(cls, Base_reaction1):
            reaction = cls(node.get('label'), enumsModule.Genre.NBody, node.get('ENDF_MT'), fissionGenre=node.get('fissionGenre', enumsModule.FissionGenre.none))
        else:
            reaction = cls(node.get('label'), enumsModule.Genre.NBody, node.get('ENDF_MT'))

        xPath.pop( )

        reaction.parseNode(node, xPath, linkData, **kwargs)

        return reaction

class Base_reaction1(Base_reaction) :
    """Base class for all types of reaction."""

    def __init__(self, label, genre, ENDF_MT, fissionGenre=enumsModule.FissionGenre.none):

        Base_reaction.__init__(self, label, genre, ENDF_MT)

        self.fissionGenre = fissionGenre

    @property
    def fissionGenre( self ) :

        return( self.__fissionGenre )

    @fissionGenre.setter
    def fissionGenre(self, value):
        """Sets the fission genre to *value* for *self*.

        :param value:   One of the allowed values in **FissionGenre**.
        """

        self.__fissionGenre = enumsModule.FissionGenre.checkEnumOrString(value)

class Base_reaction2(Base_reaction1):
    """Base class for all types of reaction."""

    ancestryMembers = ('doubleDifferentialCrossSection', ) + Base_reaction1.ancestryMembers
    keyName = 'label'

    def __init__(self, label, genre, ENDF_MT, fissionGenre=enumsModule.FissionGenre.none):

        Base_reaction1.__init__(self, label, genre, ENDF_MT, fissionGenre=fissionGenre)

        self.__doubleDifferentialCrossSection = doubleDifferentialCrossSectionModule.Component()
        self.__doubleDifferentialCrossSection.setAncestor(self)

    @property
    def doubleDifferentialCrossSection( self ) :

        return( self.__doubleDifferentialCrossSection )

def isGNDSReaction( o ) :
    """Returns True if o is an instance of base_reaction or of a subclass thereof. """

    return isinstance(o, Base_reaction)
