# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from PoPs import IDs as IDsPoPsModule
from PoPs.fissionFragmentData import rate as rateModule

from fudge import suites as suitesModule
from fudge import product as productModule

class Product(productModule.Product):

    keyName = None

class DelayedNeutron( ancestryModule.AncestryIO ) :

    moniker = 'delayedNeutron'
    keyName = 'label'
    ancestryMembers = ('rate', 'product')

    def __init__(self, label, product):

        ancestryModule.AncestryIO.__init__( self )

        self.__label = label

        self.__rate = rateModule.Suite( )
        self.__rate.setAncestor( self )

        if not isinstance(product, Product): raise TypeError('Invalid product instance "%s".' % type(product))
        self.__product = product
        self.__product.setAncestor( self )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def rate( self ) :

        return( self.__rate )

    @property
    def product( self ) :

        return( self.__product )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        self.__product.calculateAverageProductData( style, indent = indent, **kwargs )

    def check( self, info ):
        """
        Check delayed neutron rate, multiplicity and distribution

        :param info: dict
        :return: list of warnings
        """
        from fudge import warning
        warnings = []

        if self.rate[0].value < 0:
            warnings.append( warning.NegativeDelayedRate(self.rate.time, self) )

        warnings += self.product.check( info )
        return warnings

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__rate.convertUnits( unitMap )
        self.__product.convertUnits( unitMap )

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** for the **product** member.
        """

        return self.__product.fixDomains(labels, energyMin, energyMax)

    def listOfProducts(self, finalOnly=False, includeQualifier=True):
        '''
        Returns, as a python set, the list of PoPs ids for all products (i.e., outgoing particles) for *self*.

        :param includeQualifier:    If `True`, particle qualifiers are include in ids, otherwise, they are stripped from ids.

        :return:                    A python set instance.
        '''

        return self.__product.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier)

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        self.__product.processMC_cdf( style, tempInfo, indent = indent )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productName'] = self.__product.pid
        tempInfo['productLabel'] = self.__product.label
        self.__product.processMultiGroup( style, tempInfo, indent )

    def multiGroupMultiplicity(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group multiplicity for the requested label for the request product of this output channel.
        
        This is a cross section weighted multiplicity.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        return self.product.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, productID)

    def multiGroupAverageEnergy(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average energy for the requested label for the requested product.
        
        This is a cross section weighted average energy.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        return self.product.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, productID)

    def multiGroupAverageMomentum(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average momentum for the requested label for the requested product. 
        
        This is a cross section weighted average momentum.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        return self.product.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, productID)

    def multiGroupProductMatrix(self, multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder):
        """
        Returns the multi-group, product matrix for the requested label for the requested product index for the requested Legendre order.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        :param legendreOrder: Requested Legendre order.
        """
        
        return self.product.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder)

    def multiGroupProductArray(self, multiGroupSettings, temperatureInfo, particleIDs, productID):
        """
        Returns the full multi-group, total product array for the requested label for the requested product id.

        :param multiGroupSettings: MG instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        """

        return self.product.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

#        self.__rate.removeStyles( styleLabels )
        self.__product.removeStyles( styleLabels )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        XMLStringList += self.__rate.toXML_strList( indent2, **kwargs )
        if( self.__product is not None ) : XMLStringList += self.__product.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 2 ) : return( [] )
        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        product = None
        for child in element :
            if( child.tag == Product.moniker ) :
                product = Product.parseNodeUsingClass(child, xPath, linkData, **kwargs)
                break

        delayedNeutron1 = cls(element.get('label'), product)

        for child in element :
            if( child.tag == delayedNeutron1.__rate.moniker ) :
                delayedNeutron1.__rate.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == product.moniker ) :
                pass
            else :
                raise TypeError( "Encountered unknown node '%s' in %s" % ( child.tag, element.tag ) )

        xPath.pop( )

        return( delayedNeutron1 )

class DelayedNeutrons(suitesModule.ExclusiveSuite):

    moniker = 'delayedNeutrons'

    def __init__( self ) :

        suitesModule.ExclusiveSuite.__init__( self, allowedClasses = ( DelayedNeutron, ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        for delayedNeutron1 in self : delayedNeutron1.calculateAverageProductData( style, indent = indent, **kwargs )

    def iterateProducts(self):
        """
        This is an iterator that returns each delayed neutron.
        """

        for delayedNeutron in self:
            yield delayedNeutron.product

    def listOfProducts(self, finalOnly=False, includeQualifier=True):
        '''
        Returns, as a python set, the list of PoPs ids for all products (i.e., outgoing particles) for *self*.

        :param finalOnly:           If `True`, only final product ids are returned, otherwise, all are returned..
        :param includeQualifier:    If `True`, particle qualifiers are include in ids, otherwise, they are stripped from ids.

        :return:                    A python set instance.
        '''

        products = set()
        for delayedNeutron1 in self:
            products.update(delayedNeutron1.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier))

        return products

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        for delayedNeutron1 in self : delayedNeutron1.processMC_cdf( style, tempInfo, indent )

    def processMultiGroup( self, style, tempInfo, indent ) :

        for delayedNeutron1 in self : delayedNeutron1.processMultiGroup( style, tempInfo, indent )
