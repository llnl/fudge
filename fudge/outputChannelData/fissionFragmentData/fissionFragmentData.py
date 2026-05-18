# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from xData import vector as vectorModule
from xData import matrix as matrixModule
from xData import productArray as productArrayModule
from LUPY import ancestry as ancestryModule

from fudge.outputChannelData.fissionFragmentData import fissionEnergyRelease as fissionEnergyReleaseModule

from fudge.processing import transporting as transportingModule

from . import delayedNeutron as delayedNeutronModule
from . import productYield as productYieldModule

class FissionFragmentData(ancestryModule.AncestryIO_base):

    moniker = 'fissionFragmentData'
    ancestryMembers = ('delayedNeutrons', 'fissionEnergyReleases', 'productYields')

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

        self.__delayedNeutrons = delayedNeutronModule.DelayedNeutrons( )
        self.__delayedNeutrons.setAncestor( self )

        self.__fissionEnergyReleases = fissionEnergyReleaseModule.FissionEnergyReleases( )
        self.__fissionEnergyReleases.setAncestor( self )

        self.__productYields = productYieldModule.Suite( )
        self.__productYields.setAncestor( self )

    @property
    def delayedNeutrons( self ) :

        return( self.__delayedNeutrons )

    @property
    def fissionEnergyReleases( self ) :

        return( self.__fissionEnergyReleases )

    @property
    def productYields( self ) :

        return( self.__productYields )

    def amendForPatch( self, fromLabel, toLabel ) :

        self.__delayedNeutrons.amendForPatch( fromLabel, toLabel )
        self.__fissionEnergyReleases.amendForPatch( fromLabel, toLabel )
        self.__productYields.amendForPatch( fromLabel, toLabel )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        self.__delayedNeutrons.calculateAverageProductData( style, indent = indent, **kwargs )

    def check( self, info ):
        """
        Check for problems in delayed fission products, product yields and energy release data

        :param info:
        :return:
        """
        from fudge import warning
        warnings = []

        for term in ('productYields', 'delayedNeutrons'):
            for data in getattr(self, term):
                dwarnings = data.check( info )
                if dwarnings:
                    warnings.append( warning.Context('%s: %s' % (term, data.label), dwarnings) )

        if self.fissionEnergyReleases:
            energyReleaseForm = self.fissionEnergyReleases.evaluated
            fer_warnings = energyReleaseForm.check(info)
            if fer_warnings:
                warnings.append( warning.Context(
                    '%s: %s' % (self.fissionEnergyReleases.moniker, energyReleaseForm.label), fer_warnings) )
        return warnings

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__delayedNeutrons.convertUnits( unitMap )
        self.__fissionEnergyReleases.convertUnits( unitMap )
        self.__productYields.convertUnits( unitMap )

    def cullStyles(self, styleList):
        """ See documentation for reactionSuite.cullStyles. """

        self.__delayedNeutrons.cullStyles(styleList)
        self.__fissionEnergyReleases.cullStyles(styleList)
        self.__productYields.cullStyles(styleList)

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** for the **delayedNeutrons** member.
        """

        return self.__delayedNeutrons.fixDomains(labels, energyMin, energyMax)

    def iterateProducts(self):
        """
        This is an iterator that returns each delayed neutron.
        """

        yield from self.__delayedNeutrons.iterateProducts()

    def listOfProducts(self, finalOnly=False, includeQualifier=True):
        '''
        Returns, as a python set, the list of PoPs ids for all products (i.e., outgoing particles) for *self*.

        :param finalOnly:           If `True`, only final product ids are returned, otherwise, all are returned..
        :param includeQualifier:    If `True`, particle qualifiers are include in ids, otherwise, they are stripped from ids.

        :return:                    A python set instance.
        '''

        return self.__delayedNeutrons.listOfProducts(finalOnly=finalOnly, includeQualifier=includeQualifier)

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        self.__delayedNeutrons.processMC_cdf( style, tempInfo, indent = indent )

    def processMultiGroup( self, style, tempInfo, indent ) :

        self.__delayedNeutrons.processMultiGroup( style, tempInfo, indent )
        self.__fissionEnergyReleases.processMultiGroup( style, tempInfo, indent )

    def multiGroupQ(self, multiGroupSettings, temperatureInfo):
        """
        Returns the sum of the multi-group, Q for the requested label for the this output channel. This is a cross section weighted Q.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        """

        Q = vectorModule.Vector()
        if multiGroupSettings.delayedNeutrons == transportingModule.DelayedNeutrons.on and len(self.fissionEnergyReleases)>0:
            fissionEnergyRelease = multiGroupSettings.form(self.__fissionEnergyReleases, temperatureInfo)
            if fissionEnergyRelease is not None:
                Q += fissionEnergyRelease.multiGroupQ(multiGroupSettings)

        return Q

    def multiGroupMultiplicity(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group multiplicity for the requested label for the request product of this output channel. This is a cross section weighted multiplicity.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        multiplicity = vectorModule.Vector()
        if multiGroupSettings.delayedNeutrons == transportingModule.DelayedNeutrons.on:
            for delayedNeutrons in self.delayedNeutrons:
                multiplicity += delayedNeutrons.multiGroupMultiplicity(multiGroupSettings, temperatureInfo, productID)

        return multiplicity

    def multiGroupAverageEnergy(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average energy for the requested label for the requested product. This is a cross section weighted average energy.
        
        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """        
        averageEnergy = vectorModule.Vector()
        if multiGroupSettings.delayedNeutrons == transportingModule.DelayedNeutrons.on:
            for delayedNeutrons in self.delayedNeutrons:
                averageEnergy += delayedNeutrons.multiGroupAverageEnergy(multiGroupSettings, temperatureInfo, productID)

        return averageEnergy

    def multiGroupAverageMomentum(self, multiGroupSettings, temperatureInfo, productID):
        """
        Returns the sum of the multi-group, average momentum for the requested label for the requested product. This is a cross section weighted average momentum.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param productID: Particle id for the requested product.
        """

        averageMomentum = vectorModule.Vector()
        if multiGroupSettings.delayedNeutrons == transportingModule.DelayedNeutrons.on:
            for delayedNeutrons in self.delayedNeutrons:
                averageMomentum += delayedNeutrons.multiGroupAverageMomentum(multiGroupSettings, temperatureInfo, productID)

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
        if multiGroupSettings.delayedNeutrons == transportingModule.DelayedNeutrons.on:
            for delayedNeutrons in self.delayedNeutrons:
                productMatrix += delayedNeutrons.multiGroupProductMatrix(multiGroupSettings, temperatureInfo, particleIDs, productID, legendreOrder)

        return productMatrix

    def multiGroupProductArray(self, multiGroupSettings, temperatureInfo, particleIDs, productID):
        """
        Returns the full multi-group, total product array for the requested label for the requested product id.

        :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param particleIDs: The list of particles to be transported.
        :param productID: Particle id for the requested product.
        """

        productArray = productArrayModule.ProductArray()

        if multiGroupSettings.delayedNeutrons == transportingModule.DelayedNeutrons.on:
            for delayedNeutrons in self.delayedNeutrons:
                productArray += delayedNeutrons.multiGroupProductArray(multiGroupSettings, temperatureInfo, particleIDs, productID)

        return productArray

    def replicate( self, other ) :

        self.__delayedNeutrons = other.delayedNeutrons
        self.__delayedNeutrons.setAncestor( self )

        self.__fissionEnergyReleases = other.fissionEnergyReleases
        self.__fissionEnergyReleases.setAncestor( self )

        self.__productYields = other.productYields
        self.__productYields.setAncestor( self )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__delayedNeutrons.removeStyles( styleLabels )
        self.__fissionEnergyReleases.removeStyles( styleLabels )
        self.__productYields.removeStyles( styleLabels )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__delayedNeutrons.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__fissionEnergyReleases.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__productYields.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 1 ) : XMLStringList = []

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element :
            if( child.tag == self.__delayedNeutrons.moniker ) :
                self.__delayedNeutrons.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == self.__fissionEnergyReleases.moniker ) :
                self.__fissionEnergyReleases.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == self.__productYields.moniker ) :
                self.__productYields.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError( "Encountered unknown node '%s' in %s" % ( child.tag, element.tag ) )

        xPath.pop( )
