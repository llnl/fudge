# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
In GNDS, an output channel (dubbed outputChannel) contains all outgoing products from either a reaction (also called an interaction)
or since a product may subsequently decay, the decay products for a for. This means that a reaction node has an outputChannel and a product may
have an outputChannel.

For example, consider the following nuclear reaction:

        n + O16 --> n + ( O16_e6 -> He4 + C12 )

This interaction represents a neutron hitting an Oxygen-16 target producing a n and Oxygen-16 in the 6th excited. The
excited Oxygen-16 subsequently decays into He4 and C12. The output channel for this interaction is "n + O16_e6" and
the "O16_e6" contains the decay channel "He4 + C12".

Output and decay channels can be further divided into twoBody or NBody channels.

    :twoBody:       The channel initially has two products whose angular and energy correlates are known. Only the 
                    center-of-mass angular data for one of the products is needed. The products' outgoing energies can 
                    be calculated from their center-of-mass angular information and other data like the masses of the 
                    products involved in a reaction or decay.

    :NBody:         The channel has two or more products whose angular and energy correlates are not known. The full double
                    differential data for each product must be given.
"""

from PoPs import IDs as IDsPoPsModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import vector as vectorModule
from xData import matrix as matrixModule
from xData import productArray as productArrayModule

from . import enums as enumsModule
from . import reactionProducts as reactionProductsModule
from .outputChannelData import Q as QModule

class Processes :

    continuum = 'continuum'
    inclusive = 'inclusive'


class OutputChannel( ancestryModule.AncestryIO ) :
    """This is the GNDS outputChannel which is a Q value and a list of products (i.e., product objects)."""

    moniker = 'outputChannel'
    ancestryMembers = ('Q', 'products', 'fissionFragmentData')

    def __init__( self, genre, process = None ) :
        """Constructor for outputChannel."""

        from . import product as productModule
        from .outputChannelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule

        ancestryModule.AncestryIO.__init__( self )

        self.genre = genre                              # Need to check for valid genre.
        self.__process = process

        self.__Q = QModule.Component( )
        self.__Q.setAncestor( self )

        self.__products = productModule.Products()
        self.__products.setAncestor( self )

        self.__fissionFragmentData = fissionFragmentDataModule.FissionFragmentData( )
        self.__fissionFragmentData.setAncestor( self )

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers :
            if( ancestryMember in ( 'fissionFragmentData' ) ) : continue
            getattr( self, ancestryMember ).findLinks( links )

    def __getitem__(self, i):
        """Returns the i^th product in the channel."""

        return self.__products[i]

    def __len__(self):
        """Returns the number of products for this channel."""

        return len(self.__products)

    def __str__(self):
        """Converts channel object to a string representation."""

        return self.toString(simpleString=False)

    @property
    def Q(self):

        return self.__Q

    @property
    def products(self):

        return self.__products

    @property
    def fissionFragmentData(self):

        return self.__fissionFragmentData

    @property
    def genre(self):

        return self.__genre

    @genre.setter
    def genre(self, value):

        self.__genre = enumsModule.Genre.checkEnumOrString(value)

    @property
    def process(self):

        return self.__process

    @process.setter
    def process(self, value):

        self.__process = value

    def amendForPatch( self, fromLabel, toLabel ) :

        self.__Q.amendForPatch( fromLabel, toLabel )
        self.__products.amendForPatch( fromLabel, toLabel )
        self.__fissionFragmentData.amendForPatch( fromLabel, toLabel )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__Q.convertUnits( unitMap )
        self.__products.convertUnits( unitMap )
        self.__fissionFragmentData.convertUnits( unitMap )

    def checkProductFrame( self ) :
        """Calls checkProductFrame for self's products."""

        for product in self : product.checkProductFrame( )

    def diff(self, other, diffResults):

        self.Q.diff( other.Q, diffResults )

        selfQ = self.Q[0].evaluate(self.Q[0].domainMin)
        otherQ = other.Q[0].evaluate(self.Q[0].domainMin)
        if abs(selfQ - otherQ) > 0.1:
            diffResults.append('Q values differ', '%s vs %s' % (selfQ, otherQ), self.toXLink(), other.toXLink())

        aliases = self.rootAncestor.PoPs.aliases
        productIDs = []
        productIDs1 = {}
        for product in self.products :
            pid = product.pid
            if( pid in aliases ) : pid = aliases[pid].pid             # Handle old GNDS files which used meta-stable ids instead of nuclear level ids.
            if( pid not in productIDs ) : productIDs.append( pid )
            if( pid not in productIDs1 ) : productIDs1[pid] = []
            productIDs1[pid].append( product )

        aliases = other.rootAncestor.PoPs.aliases
        productIDs2 = {}
        for product in other.products : 
            pid = product.pid
            if( pid in aliases ) : pid = aliases[pid].pid             # Handle old GNDS files which used meta-stable ids instead of nuclear level ids.
            if( pid not in productIDs ) : productIDs.append( pid )
            if( pid not in productIDs2 ) : productIDs2[pid] = []
            productIDs2[pid].append( product )

        for productID in productIDs :
            if( productID not in productIDs1 ) :        # This should only happen for photons and for reaction 'sumOfRemainingOutputChannels'.
                diffResults.append( 'Product missing from channel - 1', 'product id = "%s"' % productID, self.toXLink( ), other.toXLink( ) )
            elif( productID not in productIDs2 ) :      # This should only happen for photons and for reaction 'sumOfRemainingOutputChannels'.
                diffResults.append( 'Product missing from channel - 2', 'product id = "%s"' % productID, self.toXLink( ), other.toXLink( ) )
            else :
                if len(productIDs1[productID]) != len(productIDs2[productID]):
                    isspecified1 = True
                    for i1, product in enumerate(productIDs1[productID]):
                        if not product.isSpecified():
                            isspecified1 = False

                    isspecified2 = True
                    for i2, product in enumerate(productIDs2[productID]):
                        if not product.isSpecified():
                            isspecified2 = False

                    if isspecified1 and isspecified2:
                        diffResults.append('More than one product with same id', 'product id = "%s"' % productID, self.toXLink(), other.toXLink())
                    elif isspecified1:
                        diffResults.append('More than one product with same id and unspecified - 2', 'product id = "%s"' % productID, self.toXLink(), other.toXLink())
                    elif isspecified2:
                        diffResults.append('More than one product with same id and unspecified - 1', 'product id = "%s"' % productID, self.toXLink(), other.toXLink())
                    else:
                        diffResults.append('More than one product with same id and unspecified - both', 'product id = "%s"' % productID, self.toXLink(), other.toXLink())
                else :
                    for i1, product1 in enumerate(productIDs1[productID]):
                        product2 = productIDs2[productID][i1]
                        product1.diff(product2, diffResults)

    def getFinalProductList( self ) :
        """Returns a list of [ multiplicity, name ] in order for the final products of this channel."""

        mP = []
        for product in self :
            if( product.outputChannel is None ) :
                try :
                    multiplicity = product.multiplicity.getConstant( )
                except :
                    multiplicity = "(?)"                        # ????? This needs work.
                mP += [ [ multiplicity, product.pid ] ]
            else :
                mP += product.outputChannel.getFinalProductList( )
        return( mP )

    def getProductsWithName(self, name, final=True, delayedNeutrons=True):
        """Returns a list of all the channel's products with given name ('photon','n', etc)."""

        return [product for product in self.iterateProducts(final=final, delayedNeutrons=delayedNeutrons) if product.pid == name]

    def getProductWithName( self, name ) :
        """
        Return the product with given name. If no such product is found, or if more than one are found, raise an Exception.
        """

        prods = self.getProductsWithName( name )
        if( len( prods ) != 1 ) : raise Exception( "Unique product named '%s' not found in %s" % ( name, self ) )
        return( prods[0] )

    def cullStyles( self, styleList ) :

        self.__Q.cullStyles( styleList )
        self.__products.cullStyles( styleList )
        self.__fissionFragmentData.cullStyles( styleList )

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** for the *Q*, *products* and *fissionFragmentDat* members.
        """

        numberOfFixes  = self.__Q.fixDomains(labels, energyMin, energyMax)
        numberOfFixes += self.__products.fixDomains(labels, energyMin, energyMax)
        numberOfFixes += self.__fissionFragmentData.fixDomains(labels, energyMin, energyMax)

        return numberOfFixes

    def iterateProducts(self, final=True, delayedNeutrons=True):
        """
        This is an iterator that returns each product of *self*. If *final* is **True**, products with an
        output channel are not returned but, instead, their nested products are. Ergo, only final products are returned
        if *final* is **True**.

        :param final:               If **True**, only final products are returned; otherwise, only the products directly within *self* are returned.
        :param delayedNeutrons:     If **True**, delayed neutrons are also iterated over if present.
        """

        for product in self.products:
            if product.outputChannel is not None:
                if final:
                    yield from product.outputChannel.iterateProducts(final=final, delayedNeutrons=delayedNeutrons)
                else:
                    yield product
            else:
                yield product

        if delayedNeutrons:
            yield from self.__fissionFragmentData.iterateProducts()

    def listOfProducts(self, finalOnly=False, includeQualifier=True):
        '''
        Returns, as a python set, the list of PoPs ids for all products (i.e., outgoing particles) for *self*.

        :param finalOnly:           If `True`, only final product ids are returned, otherwise, all are returned..
        :param includeQualifier:    If `True`, particle qualifiers are include in ids, otherwise, they are stripped from ids.

        :return:                    A python set instance.
        '''

        products = set()
        for product in self.__products:
            products.update(product.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))
        products.update(self.__fissionFragmentData.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))

        return products

    def reactionProducts( self, _reactionProducts ) :

        for product in self.products : product.reactionProducts( _reactionProducts )
        if( self.process is not None ) : _reactionProducts[self.process] = reactionProductsModule.ReactionProduct( reactionProductsModule.Category.process, 0 )

        return( _reactionProducts )

    def removeStyles( self, styleLabels ) :

        self.__Q.removeStyles( styleLabels )
        self.__products.removeStyles( styleLabels )
        self.__fissionFragmentData.removeStyles( styleLabels )

    def thresholdQAs(self, unit, final=True, pops=None):

        if len(self.__Q) > 0 and pops is None:
            Q = self.__Q.thresholdQAs(unit)
        else:                       # Calculate from particle masses.
            massUnit = unit + '/c**2'
            reactionSuite = self.rootAncestor
            projectile = pops.final(reactionSuite.projectile)
            target = pops.final(reactionSuite.target)
            Q = target.getMass(massUnit) + projectile.getMass(massUnit)
            for product in self:
                try:
                    Q -= pops[product.pid].getMass(massUnit) * product.multiplicity.getConstant()
                except:
                    if product.pid == IDsPoPsModule.photon:
                        continue
                    raise ValueError('Non-constant Q-value must be explicitly listed in GNDS! Maybe missing mass for particle "%s".' % product.pid)

        if final:
            for product in self.__products:
                Q += product.thresholdQAs(unit, final=final, pops=pops)

        return Q

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        kwargs['outputChannel'] = self

        for productIndex, product in enumerate( self.__products ) :
            kwargs['productIndex'] = str( productIndex )
            product.calculateAverageProductData( style, indent = indent, **kwargs )

        self.__fissionFragmentData.calculateAverageProductData( style, indent, **kwargs )

    def partialProductionIntegral( self, reaction_suite, productID, energyIn, energyOut = None, muOut = None, phiOut = None, 
                frame = xDataEnumsModule.Frame.lab, LegendreOrder = 0, **kwargs ) :

        partialProductionIntegralSum = 0.0
        for product in self.__products :
            partialProductionIntegralSum += product.partialProductionIntegral( reaction_suite, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                    phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )

        return( partialProductionIntegralSum )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']

        for productIndex, product in enumerate( self ) :
            tempInfo['productIndex'] = str( productIndex )
            product.processMC_cdf( style, tempInfo, indent2 )

        self.__fissionFragmentData.processMC_cdf( style, tempInfo, indent2 )

    def processMultiGroup( self, style, tempInfo, indent ) :

        indent2 = indent + tempInfo['incrementalIndent']

        self.__Q.processMultiGroup( style, tempInfo, indent )

        for productIndex, product in enumerate( self ) :
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productName'] = product.pid
            tempInfo['productLabel'] = product.label
            product.processMultiGroup( style, tempInfo, indent2 )

        self.__fissionFragmentData.processMultiGroup( style, tempInfo, indent2 )

    def multiGroupQ(self, multiGroupSettings, temperatureInfo, includeFinalProducts):
        """
        Returns the sum of the multi-group, Q for the requested label for the this output channel. This is a cross section weighted Q.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param includeFinalProducts: Boolean value indicating whether to include the contriibution from the final fission product data.
        """

        Q = multiGroupSettings.formAsVector(self.__Q, temperatureInfo)

        if includeFinalProducts:
            for product in self.products:
                Q += product.multiGroupQ(multiGroupSettings, temperatureInfo, includeFinalProducts)

        Q += self.__fissionFragmentData.multiGroupQ(multiGroupSettings, temperatureInfo)

        return Q

    def multiGroupMultiplicity(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group multiplicity for the requested label for the request product of this output channel. This is a cross section weighted multiplicity.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        multiplicity = vectorModule.Vector()
        for product in self.products:
            multiplicity += product.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, productID)

        multiplicity += self.fissionFragmentData.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, productID)

        return multiplicity

    def multiGroupAverageEnergy(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average energy for the requested label for the requested product. This is a cross section weighted average energy.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        averageEnergy = vectorModule.Vector()
        for product in self.products:
            averageEnergy += product.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, productID)

        averageEnergy += self.fissionFragmentData.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, productID)

        return averageEnergy

    def multiGroupAverageMomentum(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average momentum for the requested label for the requested product. This is a cross section weighted average momentum.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        averageMomentum = vectorModule.Vector()
        for product in self.products:
            averageMomentum += product.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, productID)

        averageMomentum += self.fissionFragmentData.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, productID)

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

        productMatrix = matrixModule.Matrix()
        for product in self.products:
            productMatrix += product.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder)

        productMatrix += self.fissionFragmentData.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder)

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

        for product in self.products:
            productArray += product.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

        productArray += self.fissionFragmentData.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

        return productArray

    def QToPointwiseLinear( self, final = True, **kwargs ) :

        linearQ = self.__Q.toPointwise_withLinearXYs( **kwargs )
        if( final ) :
            for product in self.__products :
                if( product.outputChannel is not None ) : linearQ += product.outputChannel.Q.toPointwise_withLinearXYs( final = final, **kwargs )
        return( linearQ )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        qualifiers = ''
        if( self.process is not None ) : qualifiers += ' process="%s"' % str( self.process )
        xmlStringList = [ '%s<%s genre="%s"%s>' % ( indent, self.moniker, self.genre, qualifiers ) ]
        xmlStringList += self.__Q.toXML_strList( indent2, **kwargs )
        xmlStringList += self.__products.toXML_strList( indent2, **kwargs )
        xmlStringList += self.__fissionFragmentData.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def channelParticlesToString( self, prefix = "", suffix = "", MT = 0 ) :

        if self.genre == enumsModule.Genre.sumOfRemainingOutputChannels:
            return str(enumsModule.Genre.sumOfRemainingOutputChannels)

        return( prefix + self.toString( simpleString = True ) + suffix )

    def toStandardString( self ) :
        """
        Returns a more standard nuclear, but less informative, string representation for the channel.
        For output channel it is list of "products)Residual" (e.g., "2n)Pu_238")."""

        mP = self.getFinalProductList( )
        s = ''
        for m, p in mP :
            if( m == 1 ) :
                sp = p
            else :
                sm = m
                if( type( m ) != type( '' ) ) : sm = str( m )
                sp = '%s%s' % ( sm, p )
            if( s == '' ) :
                s = sp
            else :
                s = s + ' + %s' % sp
        return( s )

    def toString( self, indent = '', simpleString = False, exposeGammaMultiplicity = False, MT = 0 ) :
        """Returns a string representation of self. If simpleString is True, the string contains only the initial
        products, and not any decay products."""

        if self.genre == enumsModule.Genre.sumOfRemainingOutputChannels:
            if self.process is None:
                return str(enumsModule.Genre.sumOfRemainingOutputChannels)

        s = ''
        p = ''
        gammaString = None
        for product in self.__products :
            if( product.pid == IDsPoPsModule.photon ) :
                gammaStringP = product.toString( simpleString = True, exposeGammaMultiplicity = exposeGammaMultiplicity )
                if( gammaString is None ) : gammaString = gammaStringP
                if( 'energyDependent' in gammaStringP ) : gammaString = gammaStringP
            else :
                s += '%s%s' % ( p, product.toString( simpleString = simpleString, exposeGammaMultiplicity = exposeGammaMultiplicity ) )
                p = ' + '
        if( gammaString is not None ) : s += '%s%s' % ( p, gammaString )
        if( self.process is not None ) : s += " [%s]" % self.process
        return( indent + s )

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        self.genre = node.get('genre')
        self.process = node.get('process')
        self.fissionGenre = node.get('fissionGenre', enumsModule.FissionGenre.none)      # This is needed as some legacy files have the fissionGenre on the outputChannel 
                                                                # instead of the reaction. Removing this will also require a change to parseNode in reactions/base.py.

        childNodesNotParse, membersNotFoundInNode = self.parseAncestryMembers(node, xPath, linkData, **kwargs)

        if len(childNodesNotParse) != 0: raise Exception("Encountered unexpected child nodes '%s' in %s!" % (self.moniker, ', '.join(list(childNodesNotParse.keys()))))

        xPath.pop( )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Translate '<outputChannel>' from xml."""

        xPath.append( element.tag )

        outputChannel = cls(enumsModule.Genre.NBody)

        xPath.pop( )

        outputChannel.parseNode(element, xPath, linkData, **kwargs)

        return outputChannel
