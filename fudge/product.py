# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the product class.

GNDS classes store reaction data in a hierarchical format. At the top is the product class. Next is the
channel, representing the reaction output channels. An output channel contains the Q-value and a list
of outgoing products (some of which may decay into other particles).
The channel is stored inside the reaction class, which consists of an input and output channel together,
(e.g. n_1 + O_16  --> H_1 + N_16). Finally, the reactionSuite class is used to store a list of reaction
for a given input channel (that is, for the incoming particles that come together in the reaction).
The following outline summarizes the classes.

:product:
    Outgoing particle, including multiplicity and distribution information.

:channel:
    A Q-value plus a list of products.

:reaction:
    Contains a cross section (defined in the gnds.reactionData.crossSection module), along with
    an outgoing channel which must be a outputChannel instance (e.g., 2n + Pu239).
    
:reactionSuite:
    Contains all reactions sharing the same input channel (e.g. n + Pu239)
    
        (e.g., n + Pu239 --> n + Pu239, n + Pu239 --> n + Pu239_e1,
        n + Pu239 --> n + Pu239_e2, n + Pu239 --> 2n + Pu239 ... )
"""

from PoPs import IDs as IDsPoPsModule
from PoPs import misc as miscPoPsModule
from PoPs.decays import misc as miscDecaysPoPsModule
from PoPs.families import particle as PoPsParticleModule

from fudge.core.utilities import brb

from xData import enums as xDataEnumsModule
from xData import vector as vectorModule
from xData import matrix as matrixModule
from xData import productArray as productArrayModule
from LUPY import ancestry as ancestryModule

from fudge import enums as enumsModule
from fudge import reactionProducts as reactionProductsModule
from fudge import outputChannel as outputChannelModule
from fudge import suites as suitesModule

from .productData.distributions import distribution as distributionModule
from .productData.distributions import unspecified as unspecifiedModule
from .productData.distributions import miscellaneous as miscellaneousModule
from .productData import multiplicity as multiplicityModule
from .productData import averageProductEnergy as averageProductEnergyModule
from .productData import averageProductMomentum as averageProductMomentumModule


class Product( ancestryModule.AncestryIO ) :
    """
    This is the class for a gnds product. If the product can decay (e.g. for breakup reactions),
    the resulting decay information is defined in the product outputChannel
    """

    moniker = 'product'
    ancestryMembers = ( 'multiplicity', 'distribution', 'outputChannel', 'averageProductEnergy', 'averageProductMomentum' )
    keyName = 'label'

    def __init__( self, pid, label, outputChannel = None ) :
        """Creates a new product object."""

        ancestryModule.AncestryIO.__init__( self )

        self.__pid = pid
        self.label = label
        self.__particle = None  # lazy evaluation

        self.__outputChannel = None
        if( outputChannel is not None ) : self.addOutputChannel( outputChannel )

        self.__multiplicity = multiplicityModule.Component( )
        self.__multiplicity.setAncestor( self )

        self.__averageProductEnergy = averageProductEnergyModule.Component( )
        self.__averageProductEnergy.setAncestor( self )

        self.__averageProductMomentum = averageProductMomentumModule.Component( )
        self.__averageProductMomentum.setAncestor( self )

        self.__distribution = distributionModule.Component( )
        self.__distribution.setAncestor( self )

    def  __str__( self ) :
        """Converts product object to a string representation."""

        return( self.toString( simpleString = False ) )

    @property
    def pid( self ) :

        return( self.__pid )

    @property
    def outputChannel( self ) :

        return( self.__outputChannel )
    
    @property
    def particle( self ) :

        if( self.__particle is None ) :
            pops = self.findAttributeInAncestry( 'PoPs' )
            try :
                self.__particle = pops[self.pid]
            except :
                baseName, anti, qualifier = miscPoPsModule.baseAntiQualifierFromID( self.pid, qualifierAllowed = True )
                name = baseName + anti
                self.__particle = pops.chemicalElements.getSymbol( name )
        return( self.__particle )

    @property
    def label( self ) :
        """Returns self's label."""

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = value

    @property
    def parentProduct( self ) :
        """Returns the parent product of self or None, if self does not have a parent product."""

        try :
            return( self.ancestor.findClassInAncestry( self.__class__ ) )
        except :
            return( None )

    @property
    def multiplicity( self ) :

        return( self.__multiplicity )

    @property
    def averageProductEnergy(self):

        return self.__averageProductEnergy

    @property
    def averageProductMomentum(self):

        return self.__averageProductMomentum

    @property
    def distribution( self ) :

        return( self.__distribution )

    def addOutputChannel( self, outputChannel ) :
        """Adds outputChannel to particle."""

        if isinstance(outputChannel, outputChannelModule.OutputChannel):
            self.__outputChannel = outputChannel
            self.__outputChannel.setAncestor( self )
        else :
            raise TypeError( 'Invalid decay channel = %s' % brb.getType( outputChannel ) )

    def amendForPatch( self, fromLabel, toLabel ) :

        self.__multiplicity.amendForPatch( fromLabel, toLabel )
        self.__averageProductEnergy.amendForPatch( fromLabel, toLabel )
        self.__averageProductMomentum.amendForPatch( fromLabel, toLabel )
        self.__distribution.amendForPatch( fromLabel, toLabel )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.multiplicity.convertUnits( unitMap )
        self.averageProductEnergy.convertUnits( unitMap )
        self.averageProductMomentum.convertUnits( unitMap )
        self.distribution.convertUnits( unitMap )
        if( self.__outputChannel is not None ) : self.__outputChannel.convertUnits( unitMap )

    def checkProductFrame( self ) :
        """
        Calls checkProductFrame for self's distributions and if present for its outputChannel.
        """

        self.distribution.checkProductFrame( )
        if( self.__outputChannel is not None ) : self.__outputChannel.checkProductFrame( )

    def cullStyles( self, styleList ) :

        self.__multiplicity.cullStyles( styleList )
        self.__distribution.cullStyles( styleList )
        self.__averageProductEnergy.cullStyles( styleList )
        self.__averageProductMomentum.cullStyles( styleList )

    @property
    def domainMin( self ) :

        return( self.multiplicity.domainMin )

    @property
    def domainMax( self ) :

        return( self.multiplicity.domainMax )

    @property
    def domainUnit( self ) :

        return( self.multiplicity.domainUnit )

    def diff( self, other, diffResults ) :

        self.distribution.diff( other.distribution, diffResults )

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers :
            if( getattr( self, ancestryMember ) is None ) : continue
            getattr( self, ancestryMember ).findLinks( links )

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** for the *multiplicity*, *distribution* and *outputChannel* members.
        """

        numberOfFixes  = self.__multiplicity.fixDomains(labels, energyMin, energyMax)
        numberOfFixes += self.__distribution.fixDomains(labels, energyMin, energyMax)
        if self.__outputChannel is not None: numberOfFixes += self.__outputChannel.fixDomains(labels, energyMin, energyMax)

        return numberOfFixes

    def thresholdQAs(self, unit, final=True, pops=None):

        if self.__outputChannel is not None:
            return self.__outputChannel.thresholdQAs(unit, final=final, pops=pops)

        return 0.0

    def getLevelAsFloat( self, unit, default = 0. ) :

        if( hasattr( self.particle, 'getLevelAsFloat' ) ) : return( self.particle.getLevelAsFloat( unit, default = default ) )
        return( default )

    def getMass( self, unit ) :
        """Returns the mass of the particle if possible, otherwise None is returned."""

        particle = self.particle
        PoPs = self.findAttributeInAncestry('PoPs')
        if particle.id in PoPs.aliases:
            particle = PoPs[particle.id]
        return( particle.getMass( unit ) )

    def isSpecified(self):
        '''
        Returns **True** if both the multiplicity and distribution have specified data. That is, each have at least one unspecified data set.
        '''

        return self.multiplicity.isSpecified() and self.distribution.isSpecified()

    def listOfProducts(self, finalOnly=False, includeQualifier=True):
        '''
        Returns, as a python set, the list of PoPs ids for all products (i.e., outgoing particles) for *self*.

        :param finalOnly:           If `True`, only final product ids are returned, otherwise, all are returned..
        :param includeQualifier:    If `True`, particle qualifiers are include in ids, otherwise, they are stripped from ids.

        :return:                    A python set instance.
        '''

        products = set()
        if not (finalOnly and self.__outputChannel is not None):
            pid = self.__pid
            if not includeQualifier:
                pid = miscPoPsModule.idWithQualifierRemoved(pid)
            products.add(pid)
        if self.__outputChannel is not None:
            products.update(self.__outputChannel.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))

        return products

    def calculateAverageProductData( self, style, indent, **kwargs ) :

        verbosity = kwargs['verbosity']
        indent2 = indent  + kwargs['incrementalIndent']
        indent3 = indent2 + kwargs['incrementalIndent']

        reactionSuite = kwargs['reactionSuite']
        energyUnit = kwargs['incidentEnergyUnit']
        momentumUnit = kwargs['momentumUnit']
        massUnit = kwargs['massUnit']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']

        if 'MT504_data' in kwargs and self.pid == IDsPoPsModule.photon:
            MT504_data = kwargs['MT504_data']

            axes = averageProductEnergyModule.defaultAxes(energyUnit=energyUnit)
            averageEnergy = averageProductEnergyModule.XYs1d(data=MT504_data['averageProductEnergy'], axes=axes, label=style.label)
            self.averageProductEnergy.add(averageEnergy)

            axes = averageProductMomentumModule.defaultAxes(energyUnit=energyUnit, momentumUnit=momentumUnit)
            averageMomentum = averageProductMomentumModule.XYs1d(data=MT504_data['averageProductMomentum'], axes=axes, label=style.label) 
            self.averageProductMomentum.add(averageMomentum)

            return

        kwargs['product'] = self
        try :
            kwargs['productMass'] = reactionSuite.PoPs[self.pid].getMass( massUnit )
        except :                        # Can happend when mass is not needed for evaluation and hence not stored.
            kwargs['productMass'] = None

        if( verbosity > 1 ) :
            print('%s%s: label = %s: calculating average product data' % (indent, self.pid, self.label))
        if( len( self.distribution ) > 0 ) :
            multiplicity = style.findFormMatchingDerivedStyle( self.multiplicity )
            if( not( isinstance( multiplicity, (multiplicityModule.Branching1d, multiplicityModule.Unspecified) ) ) ) :
                kwargs['multiplicity'] = multiplicity.toPointwise_withLinearXYs( accuracy = 1e-5, lowerEps = 1e-6, upperEps = 1e-6 )

            energyData, momentumData = self.distribution.calculateAverageProductData( style, indent = indent2, **kwargs )
            if( energyData is not None ) :
                axes = averageProductEnergyModule.defaultAxes( energyUnit = energyUnit )
                if( len( energyData ) == 1 ) :
                    averageEnergy = averageProductEnergyModule.XYs1d( data = energyData[0], axes = axes, label = style.label ) 
                else :
                    averageEnergy = averageProductEnergyModule.Regions1d( axes = axes, label = style.label )
                    for energyDataRegion in energyData :
                        averageEnergyRegion = averageProductEnergyModule.XYs1d( data = energyDataRegion, axes = axes )
                        averageEnergy.append( averageEnergyRegion )
                self.averageProductEnergy.add( averageEnergy )
            if( momentumData is not None ) :
                axes = averageProductMomentumModule.defaultAxes( energyUnit = energyUnit, momentumUnit = momentumUnit )
                if( len( momentumData ) == 1 ) :
                    averageMomentum = averageProductMomentumModule.XYs1d( data = momentumData[0], axes = axes, label = style.label ) 
                else :
                    averageMomentum = averageProductMomentumModule.Regions1d( axes = axes, label = style.label ) 
                    for momentumDataRegion in momentumData :
                        averageMomentumRegion = averageProductMomentumModule.XYs1d( data = momentumDataRegion, axes = axes ) 
                        averageMomentum.append( averageMomentumRegion )
                self.averageProductMomentum.add( averageMomentum )

        if( self.__outputChannel is not None ) :
            self.__outputChannel.calculateAverageProductData( style, indent = indent3, **kwargs )

    def partialProductionIntegral(self, reaction_suite, productID, energyIn, energyOut=None, muOut=None, phiOut=None, 
            frame=xDataEnumsModule.Frame.product, LegendreOrder=0, **kwargs):

        def branchingGammas(initialState, photonBranchingData, probability, LegendreOrder):

            partialProductionIntegralSum = 0.0
            if LegendreOrder != 0:
                return partialProductionIntegralSum             # Currently, assume isotropic scattering in lab frame.
            if initialState in photonBranchingData:
                gammas = photonBranchingData[initialState]['photons']
                for branchingRatio, gammaEnergy, finalState, photonEmissionProbability in gammas :
                    gammaEnergy = float(gammaEnergy)
                    energyOutMin, energyOutMax = miscellaneousModule.domainLimits(energyOut, gammaEnergy, gammaEnergy)
                    if energyOutMin <= gammaEnergy <= energyOutMax:
                        partialProductionIntegralSum += branchingRatio * probability * photonEmissionProbability
                    partialProductionIntegralSum += branchingGammas(finalState, photonBranchingData, branchingRatio * probability, LegendreOrder)

            return partialProductionIntegralSum

        partialProductionIntegralSum = 0.0

        if( self.pid == productID ) :
            if( self.distribution.hasData( ) ) :
                multiplicity = self.multiplicity[0]
                if( isinstance( multiplicity, multiplicityModule.Branching1d ) ) :
                    photonBranchingData = miscDecaysPoPsModule.photonBranchingData( reaction_suite.PoPs, self.parentProduct.pid )
                    factor = miscellaneousModule.muPhiEvaluate( muOut, phiOut )
                    partialProductionIntegralSum += factor * branchingGammas( multiplicity.product.parentProduct.pid, photonBranchingData, 1.0, LegendreOrder )
                else :
                    domainMin = self.multiplicity[0].domainMin
                    domainMax = self.multiplicity[0].domainMax
                    if( domainMin <= energyIn <= domainMax ) :
                        multiplicity = self.multiplicity.evaluate( energyIn )
                        if( multiplicity != 0.0 ) :
                            partialProductionIntegralSum = self.distribution.integrate( reaction_suite, energyIn, energyOut = energyOut, 
                                    muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder )
                            partialProductionIntegralSum *= multiplicity

        if( self.__outputChannel is not None ) :
            partialProductionIntegralSum += self.__outputChannel.partialProductionIntegral( reaction_suite, productID, energyIn, energyOut = energyOut, 
                    muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder )

        return( partialProductionIntegralSum )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']
        if( verbosity > 1 ) : print('%s%s: label = %s: MonteCarlo_cdf processing' % ( indent, self.pid, self.label ))

        self.distribution.processMC_cdf( style, tempInfo, indent )
        if( self.__outputChannel is not None ) : self.__outputChannel.processMC_cdf( style, tempInfo, indent2 )

    def processMultiGroup(self, style, tempInfo, indent):

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        tempInfo['workFile'].append(self.label)

        doIt = not isinstance(self.distribution[0], unspecifiedModule.Form)
        if doIt and self.pid in style.transportables:
            if verbosity > 1:
                print('%s%s: label = %s: multiGroup processing' % (indent, self.pid, self.label))

            productMass = tempInfo['masses']['Product']             # Save to restore later
            tempInfo['masses']['Product'] = self.getMass(tempInfo['massUnit'])
            tempInfo['product'] = self
            tempInfo['multiplicity'] = self.multiplicity

            self.multiplicity.processMultiGroup( style, tempInfo, indent )
            self.averageProductEnergy.processMultiGroup( style, tempInfo, indent )
            self.averageProductMomentum.processMultiGroup( style, tempInfo, indent )

            reprocessingOpts = tempInfo['partialReprocessingOptions']
            if reprocessingOpts is not None and not (tempInfo['crossSectionModified'] or
                        reprocessingOpts['preprocessingStyleLabel'] in self.multiplicity or
                        reprocessingOpts['preprocessingStyleLabel'] in self.distribution):
                # no changes: copy processed data from original processedDataSource
                if verbosity > 1:
                    print('%s%s: label = %s: copying transfer matrix from template file' % (indent, self.pid, self.label))

                transferMatrix = reprocessingOpts['processedDataSource'].followXPath(
                    self.toXLink()).distribution[reprocessingOpts['otherStyleLabel']]

                # FIXME: should the copy method automatically skip copying self.ancestor?
                ancestor = transferMatrix.ancestor
                transferMatrix.setAncestor(None)
                newMatrix = transferMatrix.copy()
                transferMatrix.setAncestor(ancestor)

                newMatrix.label = style.label
                self.distribution.add(newMatrix)

            else:
                try:
                    self.distribution.processMultiGroup( style, tempInfo, indent )
                except:
                    if tempInfo['logFile'] is None:
                        raise
                    else :
                        import traceback
                        tempInfo['logFile'].write( '\n' + self.toXLink() + ':\n' + traceback.format_exc( ) + '\n' )
                        tempInfo['failures'] += 1
            tempInfo['masses']['Product'] = productMass

        if self.__outputChannel is not None: self.__outputChannel.processMultiGroup(style, tempInfo, indent2)

        del tempInfo['workFile'][-1]

    def multiGroupQ(self, multiGroupSettings, temperatureInfo, includeFinalProducts):
        """
        Returns the sum of the multi-group, Q for the requested label for the this product and any sub-product.

        This is a cross section weighted Q. If includeFinalProducts is False, only the Q for the products output channel is returned, 
        otherwise, the Q for all output channels summed, including output channels for each products.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param includeFinalProducts: If true, the Q is calculated for all output channels, including those for products.
        """
        Q = vectorModule.Vector()
        if self.outputChannel is not None:
            Q += self.outputChannel.multiGroupQ(multiGroupSettings, temperatureInfo, includeFinalProducts)

        return Q

    def multiGroupMultiplicity(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, multiplicity for the requested label for the this product and any sub-product.
        
        This is a cross section weighted multiplicity.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        if self.outputChannel is None:
            if self.pid == productID:
                multiplicity = multiGroupSettings.formAsVector(self.multiplicity, temperatureInfo)
            else:
                multiplicity = vectorModule.Vector()
        else:
            multiplicity = self.outputChannel.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, productID)

        return multiplicity

    def multiGroupAverageEnergy(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average energy for the requested label for the requested product.
        
        This is a cross section weighted average energy summed over this and all sub-products.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        if self.outputChannel is None:
            if self.pid == productID:
                averageEnergy = multiGroupSettings.formAsVector(self.averageProductEnergy, temperatureInfo)
            else:
                averageEnergy = vectorModule.Vector()
        else:
            averageEnergy = self.outputChannel.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, productID)

        return averageEnergy

    def multiGroupAverageMomentum(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average momentum for the requested label for the requested product.
        
        This is a cross section weighted average momentum summed over this and all sub-products.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        if self.outputChannel is None:
            if self.pid == productID:
                averageMomentum = multiGroupSettings.formAsVector(self.averageProductMomentum, temperatureInfo)
            else:
                averageMomentum = vectorModule.Vector()
        else:
            averageMomentum = self.outputChannel.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, productID)

        return averageMomentum

    def multiGroupProductMatrix(self, multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder):
        """
        Returns the multi-group, product matrix for the requested label for the requested product index for the requested Legendre order.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        :param legendreOrder: Requested Legendre order.
        """

        if self.outputChannel is None:
            if self.pid == productID:
                productMatrix = multiGroupSettings.formAsMatrix(self.distribution, temperatureInfo, legendreOrder)
            else:
                productMatrix = matrixModule.Matrix()

        else:
            productMatrix = self.outputChannel.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder)

        return productMatrix

    def multiGroupProductArray(self, multiGroupSettings, temperatureInfo, particleIDs, productID):
        """
        Returns the full multi-group, total product array for the requested label for the requested product id.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        """

        productArray = productArrayModule.ProductArray()

        if self.outputChannel is None:
            if self.pid == productID:
                form = multiGroupSettings.form(self.distribution, temperatureInfo)
                if form is not None:
                    productArray = productArrayModule.ProductArray(form.multiGroupSubform.array.constructArray())
        else:
            productArray = self.outputChannel.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

        return productArray

    def check(self, info):
        """ check product multiplicity, distribution and breakup products (if applicable) """
        from fudge import warning
        warnings = []

        multWarnings = self.multiplicity.check(info)

        if self.pid == IDsPoPsModule.photon:
            from fudge.productData.distributions import uncorrelated, energy as energyModule
            distribution = self.distribution.evaluated
            if isinstance(distribution, uncorrelated.Form) and isinstance(distribution.energySubform.data,
                           (energyModule.DiscreteGamma, energyModule.PrimaryGamma)):
                if self.multiplicity.evaluated.rangeMax > 1:
                    multWarnings.append(warning.UnphysicalDiscreteOrPrimaryPhotonMultiplicity(self))

        if multWarnings:
            warnings.append(warning.Context("Multiplicity:", multWarnings))

        if self.label in info['transportables'] and not self.distribution.hasData():
            warnings.append(warning.MissingDistribution(self.label, self))

        distributionWarnings = self.distribution.check(info)
        if distributionWarnings:
            warnings.append(warning.Context("Distribution:", distributionWarnings))

        if self.__outputChannel is not None:
            parentIsTwoBody = info['isTwoBody']
            info['isTwoBody'] = self.__outputChannel.genre == enumsModule.Genre.twoBody
            for decayProduct in self.__outputChannel:
                decayWarnings = decayProduct.check(info)
                if decayWarnings:
                    warnings.append(warning.Context("Decay product: %s" % decayProduct.label, decayWarnings))
            info['isTwoBody'] = parentIsTwoBody # reset to parent channel

        return warnings

    def reactionProducts(self, _reactionProducts):

        category = reactionProductsModule.Category.particle
        if self.__outputChannel is not None: category = reactionProductsModule.Category.intermediate

        multiplicityForm = self.multiplicity.evaluated
        multiplicity = 0
        if isinstance(multiplicityForm, multiplicityModule.Constant1d): multiplicity = multiplicityForm.value
        if int(multiplicity) == multiplicity: multiplicity = int(multiplicity)

        pid = self.pid
        pops = self.findAttributeInAncestry('PoPs')
        pid = pops.final(pid).id                # Some legacy files have meta-stable as a product instead of a particle (e.g., Sc46_m1 instead of Sc46_e2).

        if pid not in _reactionProducts: _reactionProducts[pid] = reactionProductsModule.ReactionProduct(category, 0)
        _reactionProducts[pid].multiplicity += multiplicity

        if self.__outputChannel is not None: self.__outputChannel.reactionProducts(_reactionProducts)

    def removeStyles( self, styleLabels ) :

        self.__multiplicity.removeStyles( styleLabels )
        self.__distribution.removeStyles( styleLabels )
        self.__averageProductEnergy.removeStyles( styleLabels )
        self.__averageProductMomentum.removeStyles( styleLabels )
        if( self.__outputChannel is not None ) : self.__outputChannel.removeStyles( styleLabels )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s pid="%s" label="%s">' % ( indent, self.moniker, self.pid, self.label ) ]

        xmlString += self.multiplicity.toXML_strList( indent2, **kwargs )
        xmlString += self.distribution.toXML_strList( indent2, **kwargs )
        xmlString += self.averageProductEnergy.toXML_strList( indent2, **kwargs )
        xmlString += self.averageProductMomentum.toXML_strList( indent2, **kwargs )

        if( self.__outputChannel is not None ) :
            xmlString += self.__outputChannel.toXML_strList( indent2, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toString( self, simpleString = False, exposeGammaMultiplicity = False ) :
        """
        Returns a string representation of self. If simpleString is True, the string contains only the final 
        particles, and not any intermediate particles.
        """

        def multiplicityString( ) :

            multiplicity = ''
            if( ( self.pid != IDsPoPsModule.photon ) or exposeGammaMultiplicity ) :
                _multiplicity = self.multiplicity.evaluated
                if( isinstance( _multiplicity, multiplicityModule.Constant1d ) ) :
                    iValue = int( _multiplicity.value )
                    if( iValue == _multiplicity.value ) :
                        if( iValue != 1 ) : multiplicity = '%s' % iValue
                    else :
                        multiplicity = '%s ' % _multiplicity.value
                else :
                    multiplicity = 'm(E)*'
            return( multiplicity )

        multiplicity = multiplicityString( )
        if( simpleString == True ) :
            productName = '%s%s' % ( multiplicity, self.pid )
        else :
            productName = self.pid
            productName = '%s%s' % ( multiplicity, productName )
            if( self.__outputChannel is not None ) : productName = '(%s -> %s)' % ( productName, self.__outputChannel )
        return( productName )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Translate a <product> element from xml."""

        xPath.append( '%s[@label="%s"]' % ( node.tag, node.get( 'label' ) ) )

        attrs = dict( node.items( ) )
        prod = cls(attrs.pop('pid'), label=attrs.pop('label'))

        kwargs['membersToSkip'] = [outputChannelModule.OutputChannel.moniker]
        childNodesNotParse, membersNotFoundInNode = prod.parseAncestryMembers(node, xPath, linkData, **kwargs)

        outputChannel = childNodesNotParse.pop(outputChannelModule.OutputChannel.moniker, None)
        if outputChannel is not None:
            prod.addOutputChannel(outputChannelModule.OutputChannel.parseNodeUsingClass(outputChannel, xPath, linkData, **kwargs))

        if len(childNodesNotParse) > 0: raise Exception("Encountered unexpected child nodes '%s' in %s!" % (prod.moniker, ', '.join(list(childNodesNotParse.keys()))))

        xPath.pop( )

        return prod

class Products(suitesModule.ExclusiveSuite):

    moniker = 'products'

    def __init__( self ) :

        suitesModule.ExclusiveSuite.__init__( self, ( Product, ) )

    def uniqueLabel( self, product ) :
        """
        If product's label is the same as another product's label in self, construct a new unique label
        based on product's name appended with '__' and one or more lower case letters (i.e., 'a' to 'z').
        """

        if( product.label is None ) : product.label = product.pid
        return( suitesModule.ExclusiveSuite.uniqueLabel( self, product ) )
