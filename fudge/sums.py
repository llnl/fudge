# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the 'sum', 'summand' and 'Summands' classes
"""

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from LUPY import ancestry as ancestryModule

from xData import link as linkModule

from fudge import warning as warningModule
from fudge import suites as suitesModule

from fudge.reactionData import crossSection as crossSectionModule
from fudge.outputChannelData import Q as QModule
from fudge.productData import multiplicity as multiplicityModule


class Sums(ancestryModule.AncestryIO_bare):
    """
    Contains all summed quantities. Currently supports summed cross sections and multiplicities,
    could extend to other types of sums later.
    """

    moniker = 'sums'
    ancestryMembers = ('crossSectionSums', 'multiplicitySums')
    legacyMemberNameMapping = {'crossSections': 'crossSectionSums', 'multiplicities': 'multiplicitySums'}

    def __init__(self):

        ancestryModule.AncestryIO.__init__(self)

        self.__crossSectionSums = suitesModule.CrossSectionSums()
        self.__crossSectionSums.setAncestor(self)

        self.__multiplicitySums = suitesModule.MultiplicitySums()
        self.__multiplicitySums.setAncestor(self)

    @property
    def crossSectionSums(self):

        return self.__crossSectionSums

    @property
    def multiplicitySums(self):

        return self.__multiplicitySums

    def convertUnits(self, unitMap):

        for child in self.ancestryMembers :
            getattr(self, child).convertUnits(unitMap)

    def cullStyles(self, styleList):

        self.__crossSectionSums.cullStyles(styleList)
        self.__multiplicitySums.cullStyles(styleList)

    def findLinks(self, links):

        for ancestryMember in self.ancestryMembers: item = getattr(self, ancestryMember).findLinks(links)

    def fixDomains(self, labels, energyMax):
        """
        Calls the **fixDomains** for the **crossSectionSums** and **multiplicitySums** members.
        """

        numberOfFixes  = self.__crossSectionSums.fixDomains(labels, 0.0, energyMax)
        numberOfFixes += self.__multiplicitySums.fixDomains(labels, 0.0, energyMax)

        return numberOfFixes

    def removeStyles(self, styleLabels):
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__crossSectionSums.removeStyles(styleLabels)
        self.__multiplicitySums.removeStyles(styleLabels)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent','  ')

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        for child in self.ancestryMembers:
            val = getattr(self, child)
            xmlStringList += val.toXML_strList(indent2, **kwargs)

        if len(xmlStringList) == 1:
            return []

        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        for child in element:
            if child.tag in self.ancestryMembers:
                getattr(self, child.tag).parseNode(child, xPath, linkData, **kwargs)
            elif linkData['format'] in [ GNDS_formatVersionModule.version_1_10, GNDS_formatVersionModule.version_2_0_LLNL_3 ] and child.tag in self.legacyMemberNameMapping:
                getattr(self, self.legacyMemberNameMapping[child.tag]).parseNode(child, xPath, linkData, **kwargs)
            else:
                print("Warning: unexpected child '%s' in sums" % child.tag)

        xPath.pop()

class CrossSectionSum(ancestryModule.AncestryIO):
    """
    Stores a summed quantity (cross section, multiplicity, etc.) along with the list of what's being summed over
    """

    moniker = 'crossSectionSum'
    ancestryMembers = ('summands', 'Q', 'crossSection')
    keyName = 'label'

    def __init__(self, label, ENDF_MT):
        """
        Create a 'sum' instance
        :param label: unique string identifying this sum
        :param ENDF_MT: for converting back to ENDF
        :return:
        """

        ancestryModule.AncestryIO.__init__(self)

        self.label = label
        self.ENDF_MT = int(ENDF_MT)

        self.__summands = Summands()
        self.__summands.setAncestor(self)

        self.__Q = QModule.Component()
        self.__Q.setAncestor(self)

        self.__crossSection = crossSectionModule.Component()
        self.__crossSection.setAncestor(self)

    def __str__(self):
        return self.label

    @property
    def summands(self):

        return self.__summands

    @property
    def crossSection(self):

        return self.__crossSection

    @property
    def Q(self):

        return self.__Q

    def convertUnits(self, unitMap):
        "See documentation for reactionSuite.convertUnits."

        self.__Q.convertUnits(unitMap)
        self.__crossSection.convertUnits(unitMap)

    def findLinks(self, links):

        for ancestryMember in self.ancestryMembers: item = getattr(self, ancestryMember).findLinks(links)

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** method on the **Q** and **crossSection** members.
        """

        numberOfFixes  = self.__Q.fixDomains(labels, energyMin, energyMax)
        numberOfFixes += self.__crossSection.fixDomains(labels, energyMin, energyMax)

        return numberOfFixes

    def heatCrossSection(self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.001,
                         heatAllPoints = False, doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True,
                         setThresholdToZero = False, verbose = 0):

        return self.__crossSection.heat(temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints,
                                        doNotThin, heatBelowThreshold, heatAllEDomain, setThresholdToZero = setThresholdToZero,
                                        addToSuite = True)

    def cullStyles(self, styleList):

        self.__Q.cullStyles(styleList)
        self.__crossSection.cullStyles(styleList)

    def removeStyles(self, styleLabels):
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__Q.removeStyles(styleLabels)
        self.__crossSection.removeStyles(styleLabels)

    def processGriddedCrossSections(self, style, verbosity = 0, indent = '', incrementalIndent = '  '):

        self.__crossSection.processGriddedCrossSections(style, verbosity = verbosity, indent = indent,
                                                        incrementalIndent = incrementalIndent)

    def processMultiGroup(self, style, tempInfo, indent):

        self.__crossSection.processMultiGroup(style, tempInfo, indent)

    def sumSummands(self, label=None):
        """
        Calculate the sum by looping over summands. If a style is specified, use that style if available.
        Otherwise sum the evaluated style.
        :param style: optional string identifying what style of data to sum ('eval', 'recon', etc.)
        :return:
        """
        def getPointwiseLinearForm(component):
            if label is not None and label in component: return component[label]
            else: return component.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)

        sums = getPointwiseLinearForm(self.__summands[0].link)
        for idx in range(1, len(self.__summands)):
            pointwise = getPointwiseLinearForm(self.__summands[idx].link)
            sums, current = sums.mutualify(1e-8, 1e-8, 0, pointwise, 1e-8, 1e-8, 0)
            sums += current
        return sums

    def scaleSummandsToMatchSum(self, label=None):
        """
        rescale the cross sections in the summands list so that their sum matches the cross section in self
        :return:
        """

        if label is not None and label in self.__crossSection:
            correctXsec = self.__crossSection[label]
        else:
            correctXsec = self.__crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        summedXsec = self.sumSummands(label=label)
        if correctXsec.domain() != summedXsec.domain():
            raise warningModule.SummedCrossSectionDomainMismatch(obj=self)
        scaleFactor=correctXsec/summedXsec

        for s in self.__summands:
            label = s.link.evaluated.label
            cutHere=s.link.evaluated.domainMin
            try:
                cutHere=max(cutHere,s.link.evaluated.ancestor.ancestor.getThreshold(unit=s.link.evaluated.domainUnit))
            except: # if getting threshold fails, cut at the domainMin
                pass
            answer = s.link.evaluated.returnAsClass(s.link.evaluated,s.link.evaluated*scaleFactor).domainSlice(domainMin=cutHere)
            answer.label=label
            s.link.remove(label)
            s.link.add(answer)


    def scaleSumToMatchSummands(self, label=None):
        """
        Sum the summands cross section list, then set the sum to the cross section in self.
        If a style is specified, the sum will be stored under that style. Otherwise it replaces the
        evaluated style.
        :return:
        """
        resummed = self.sumSummands()
        if label is None: label = self.__crossSection.evaluated.label
        resummed.label = resummed.style = label
        if label == self.__crossSection.evaluated.label:
            self.__crossSection.remove( label )
        self.__crossSection.add(resummed)

    def check(self, info):

        warnings = self.__crossSection.check(info)

        # does self.__crossSection equal the sum over summand reaction cross sections?
        sums = self.sumSummands(label=info['reconstructedStyleName'])
        if info['reconstructedStyleName'] in self.__crossSection:
            quotedXsec = self.__crossSection[info['reconstructedStyleName']]
        else:
            quotedXsec = self.__crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        if sums.domain() != quotedXsec.domain():
            warnings.append(warningModule.SummedCrossSectionDomainMismatch(obj=self))
        else:
            diff = quotedXsec - sums
            diff.nf_pointwiseXY.setSafeDivide(True)
            relativeDiff = abs(diff) / quotedXsec
            try:    # FIXME instead of try/catch, should really compare relativeDiff.rangeMax with 'inf'
                if relativeDiff.rangeMax > info['crossSectionMaxDiff']:
                    warnings.append(warningModule.SummedCrossSectionMismatch(relativeDiff.rangeMax * 100, obj=self))
            except TypeError:
                warnings.append(warningModule.SummedCrossSectionZeroDivision(obj=self))

        return warnings

    def getLabel(self):

        return self.label

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s label="%s" ENDF_MT="%s">' % (indent, self.moniker, self.label, self.ENDF_MT)]
        xmlString += self.__summands.toXML_strList(indent2, **kwargs)
        xmlString += self.__Q.toXML_strList(indent2, **kwargs)
        xmlString += self.__crossSection.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append('%s[@label="%s"]' % (node.tag, node.get('label')))

        sums = cls(node.get('label'), int(node.get('ENDF_MT')))

        for child in node:
            if child.tag not in (Summands.moniker, QModule.Component.moniker, crossSectionModule.Component.moniker):
                raise TypeError('Unexpected child node "%s" encountered' % child.tag)

        sums.__summands.parseNode(node.find(Summands.moniker), xPath, linkData, **kwargs)
        sums.__crossSection.parseNode(node.find(crossSectionModule.Component.moniker), xPath, linkData, **kwargs)
        sums.__Q.parseNode(node.find(QModule.Component.moniker), xPath, linkData, **kwargs)

        xPath.pop()

        return sums

class MultiplicitySum(ancestryModule.AncestryIO):
    """
    Stores a summed multiplicity along with the list of what's being summed over
    """

    moniker = 'multiplicitySum'
    ancestryMembers = ('summands', 'multiplicity')
    keyName = 'label'

    def __init__(self, label, ENDF_MT):
        """
        Create a 'sum' instance
        :param name: descriptive string
        :param label: unique string identifying this sum
        :param ENDF_MT: for converting back to ENDF
        :return:
        """

        ancestryModule.AncestryIO.__init__(self)
        self.label = label

        self.ENDF_MT = int(ENDF_MT)

        self.__summands = Summands()
        self.__summands.setAncestor(self)

        self.__multiplicity = multiplicityModule.Component()
        self.__multiplicity.setAncestor(self)

    def __str__(self):
        return self.label

    @property
    def summands(self):

        return self.__summands

    @property
    def multiplicity(self):

        return self.__multiplicity

    def convertUnits(self, unitMap):
        "See documentation for reactionSuite.convertUnits."

        self.__multiplicity.convertUnits(unitMap)

    def findLinks(self, links):

        for ancestryMember in self.ancestryMembers: item = getattr(self, ancestryMember).findLinks(links)

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** for the **multiplicity** member.
        """

        return self.__multiplicity.fixDomains(labels, energyMin, energyMax)

    def check(self, info):

        warnings = []

        if self.__multiplicity.isConstant() and self.__multiplicity.getConstant() < 1:
            warnings.append(warningModule.NegativeMultiplicity(self.getConstant(), self))
        else:   # energy-dependent mult.
            for form in self.__multiplicity:
                if hasattr(form, 'rangeMin') and form.rangeMin < 0:
                    warnings.append(warningModule.NegativeMultiplicity(form.rangeMin, obj=form))

        # does multiplicity equal the sum over its summand multiplicities?
        sums = self.__summands[0].link.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        for idx in range(1, len(self.__summands)):
            sums, current = sums.mutualify(1e-8, 1e-8, 0,
                    self.__summands[idx].link.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8), 1e-8, 1e-8, 0)
            sums += current
        quotedXsec = self.__multiplicity.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        if sums.domain() != quotedXsec.domain():
            warnings.append(warningModule.SummedMultiplicityDomainMismatch(obj=self))
        else:
            diff = quotedXsec - sums
            diff.nf_pointwiseXY.setSafeDivide(True)
            relativeDiff = abs(diff) / quotedXsec
            try:    # FIXME instead of try/catch, should really compare relativeDiff.rangeMax with 'inf'
                if relativeDiff.rangeMax > info['multiplicityMaxDiff']:
                    warnings.append(warningModule.SummedMultiplicityMismatch(relativeDiff.rangeMax * 100, obj=self))
            except TypeError:
                warnings.append(warningModule.SummedMultiplicityZeroDivision(obj=self))

        return warnings

    def getLabel(self):

        return self.label

    def cullStyles(self, styleList):

        self.__multiplicity.cullStyles(styleList)

    def removeStyles(self, styleLabels):
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__multiplicity.removeStyles(styleLabels)

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attrStr = ""
        xmlString = ['%s<%s label="%s" ENDF_MT="%s"%s>' % (indent, self.moniker, self.label, self.ENDF_MT, attrStr)]
        xmlString += self.__summands.toXML_strList(indent2, **kwargs)
        xmlString += self.__multiplicity.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append('%s[@label="%s"]' % (node.tag, node.get('label')))

        for key in node.keys():
            if key == 'ENDF_MT':
                MT = int(node.get(key))
            elif key == 'label':
                label = node.get(key)
            else:
                raise ValueError('Unsupported attribute "%s"' % key)

        sums = MultiplicitySum(label, MT)

        for child in node:
            if child.tag not in (Summands.moniker, multiplicityModule.Component.moniker):
                raise TypeError('Unexpected child node "%s" encountered' % child.tag)

        sums.__summands.parseNode(node.find(Summands.moniker), xPath, linkData, **kwargs)
        sums.__multiplicity.parseNode(node.find(multiplicityModule.Component.moniker), xPath, linkData, **kwargs)

        xPath.pop()

        return sums

class Summands(ancestryModule.AncestryIO_base):

    moniker = 'summands'

    def __init__(self):

        ancestryModule.AncestryIO.__init__(self)
        self.__summands = []

    def __getitem__(self, index):

        return self.__summands[index]

    def __len__(self):

        return len(self.__summands)

    @property
    def summands(self):

        return self.__summands

    def append(self, item):

#        if( not( isinstance( iter, ( Add, ) ) ) ) : raise TypeError( 'Invalid item to add: type = %s.' % type( item ) )

        self.__summands.append(item)
        item.setAncestor(self)

    def findLinks(self, links):

        for item in self:
            if isinstance(item, (linkModule.Link, linkModule.Link2)):
                links.append([item, item.link, item.path])

    def index(self, item):
        """Returns the index of item in self."""

        return self.__summands.index(item)

    def pop(self, index = -1):
        """Removes and return item at index. Default is the last item."""

        return self.__summands.pop(index)

    def toXML_strList(self, indent = '', **kwargs):

        xmlList = ['%s<%s>' % (indent, self.moniker)]
        for summand in self.__summands:
            xmlList.append('  %s%s' % (indent, summand.toXML()))
        xmlList[-1] += '</%s>' % self.moniker
        return xmlList

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        for child in element :
            if child.tag == Add.moniker:
                self.append(Add.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            else:
                raise TypeError('Unsupported child node "%s" of node "%s".' % (child.tag, self.moniker))

        xPath.pop()

class Add(linkModule.Link):
    """Link representing one of the quantities that is added to the sum."""

    moniker = 'add'
