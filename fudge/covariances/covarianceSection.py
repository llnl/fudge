# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""A covarianceSuite is organized into sections, where each section contains either
    - a covariance matrix for a single reaction quantity (cross section, multiplicity, etc), or
    - a covariance matrix between two different quantities (off-diagonal block)
    """
import abc

from LUPY import ancestry as ancestryModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import link as linkModule

from fudge import suites as suitesModule

from .base import Covariance
from .covarianceMatrix import CovarianceMatrix
from .mixed import MixedForm
from .summed import SummedCovariance


class CovarianceSection(suitesModule.Suite):
    """
    A covarianceSuite contains sections, where each section represents either a self-covariance for one quantity,
    or a cross-covariance between two quantities

    More generally, the covarianceSuite can be thought of as a single covariance matrix with all covariance data
    for a target/projectile. It is broken into sections, where each section holds a chunk of the full matrix.

    Within each section, covariance data can take multiple forms: :py:class:`covarianceMatrix` is the most common,
    but 'summed', 'mixed' are also possible. 

    section inherits from 'Suite', and can contain anything inheriting from the base covariance class.
    """

    moniker = 'covarianceSection'
    monikerByFormat = {GNDS_formatVersionModule.version_1_10: 'section'}
    keyName = 'label'
    ancestryMembers = ('rowData', 'columnData')

    def __init__(self, label, rowData=None, columnData=None):
        """
        Defines one rectangular block of the full covariance matrix. Contains one or more forms
        (each form inherits from base.covarianceMatrix)

        :param label: str label that identifies the section.
        :param rowData: xData.link.link pointing to data corresponding to rows of the covariance matrix
        :param columnData: xData.link.link pointing to data corresponding to columns of the covariance matrix
        """

        suitesModule.Suite.__init__(self, [Covariance])
        self.label = label

        self.rowData = rowData
        self.rowData.setAncestor(self)

        self.columnData = columnData
        if self.columnData is not None:
            self.columnData.setAncestor(self)

    @property
    def crossTerm(self):
        return self.columnData is not None and self.columnData != self.rowData

    @property
    def evaluated(self):
        """
        Helper method to grab evaluated style
        FIXME method should be inherited, but abstractClasses.Component defines methods that don't make sense for covariances
        """

        if not hasattr(self.rootAncestor, 'styles'):
            return self[0]  # For covarianceSection that is not part of a full covarianceSuite.

        evaluated = self.rootAncestor.styles.getEvaluatedStyle()
        try:
            return self[evaluated.label]
        except IndexError:
            return self[0]

    def check(self, info):
        """ check each section """

        from fudge import warning
        warnings = []
        for form in self:
            formWarnings = form.check(info)
            if formWarnings:
                warnings.append(warning.Context("Form '%s':" % form.label, formWarnings))

        return warnings

    def findInstancesOfClassInChildren(self, cls, level=9999):

        foundInstances = ancestryModule.Ancestry.findInstancesOfClassInChildren(self, cls, level)
        foundInstances += suitesModule.Suite.findInstancesOfClassInChildren(self, cls, level)

        return foundInstances

    def fix(self, **kw):
        """assemble some useful info, to be handed down to children's check() functions"""
        info = {}
        warnings = []
        if self.rowData is None:
            info['rowENDF_MFMT'] = None
        else:
            info['rowENDF_MFMT'] = self.rowData['ENDF_MFMT']
        if self.columnData is None:
            info['columnENDF_MFMT'] = None
        else:
            info['columnENDF_MFMT'] = self.columnData['ENDF_MFMT']
        info.update(kw)
        for form in self: warnings += form.fix(**info)
        return warnings

    def convertUnits(self, unitMap):
        """
        :param unitMap:     A dictionary where each key is a unit that will be replaced by its value. Units must be compatible.
        """
        for member in self.ancestryMembers:
            if getattr(self, member):
                getattr(self, member).convertUnits(unitMap)
        suitesModule.Suite.convertUnits(self, unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        moniker = self.monikerByFormat.get(formatVersion, self.moniker)
        xmlString = [indent + '<%s label="%s"' % (moniker, self.label)]
        if self.crossTerm: xmlString[0] += ' crossTerm="true"'
        xmlString[0] += '>'
        for dataPointer in ('rowData', 'columnData'):
            if getattr(self, dataPointer) is not None:
                xmlString.append(getattr(self, dataPointer).toXML(indent2, **kwargs))
        for form in self:
            xmlString += form.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Translate <section> element from xml."""

        xPath.append('%s[@label="%s"]' % (element.tag, element.get('label')))

        linkData['typeConversion'] = {'domainMin': float, 'domainMax': float}
        rowData_ = RowData.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        columnData_ = None
        if element[1].tag == "columnData":
            columnData_ = ColumnData.parseNodeUsingClass(element[1], xPath, linkData, **kwargs)
        del linkData['typeConversion']
        section_ = cls(element.get('label'), rowData_, columnData_)
        start = 2 if (columnData_ is not None) else 1
        for form in element[start:]:
            formClass = {
                CovarianceMatrix.moniker: CovarianceMatrix,
                MixedForm.moniker: MixedForm,
                SummedCovariance.moniker: SummedCovariance,
            }.get(form.tag)
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % form.tag)
            section_.add(formClass.parseNodeUsingClass(form, xPath, linkData, **kwargs))

        xPath.pop()

        return section_


class DataLink(linkModule.Link, abc.ABC):
    """
    Base class for RowData and ColumnData. Both are links but with some additional attributes.
    """
    ancestryMembers = ('slices',)

    def __init__(self, link=None, root=None, path=None, label=None, relative=False, ENDF_MFMT=None, dimension=None):
        linkModule.Link.__init__(self, link=link, root=root, path=path, label=label, relative=relative)
        self.ENDF_MFMT = ENDF_MFMT
        if dimension is not None:
            dimension = int(dimension)
        self.dimension = dimension

        self.slices = Slices()

    def __eq__(self, other):
        for attr in ('linkWithoutUpdating', 'path', 'root', 'ENDF_MFMT', 'dimension'):
            if getattr(self, attr) != getattr(other, attr):
                return False
        if len(self.slices) != len(other.slices):
            return False
        for s1, s2 in zip(self.slices, other.slices):
            if s1 != s2:
                return False
        return True

    def convertUnits(self, unitMap):
        for slice in self.slices:
            slice.convertUnits(unitMap)

    """ # not sure we need this anymore...
    def __deepcopy__( self, memodict={} ):

        copy_ = super().__deepcopy__()
        for attr in ('ENDF_MFMT', 'L', 'domainMin', 'domainMax', 'domainUnit'):
            setattr(copy_, attr, getattr(self, attr))
        return copy_

    copy = __deepcopy__
    """

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)

        if formatVersion == GNDS_formatVersionModule.version_1_10:
            attributeStr = ''
            if self.ENDF_MFMT is not None:
                attributeStr = ' ENDF_MFMT="%s"' % self.ENDF_MFMT
            if self.slices:
                if len(self.slices) > 1:
                    raise Exception("GNDS-1.10 did not support covariances for functions of 3 dimensions or more")
                slice_ = self.slices[0]
                if slice_.domainValue is not None:
                    attributeStr += ' L="%s"' % slice_.domainValue
                for attr in ('domainMin', 'domainMax', 'domainUnit'):
                    if getattr(slice_, attr) is not None:
                        attributeStr += ' %s="%s"' % (attr, getattr(slice_, attr))
            XMLList = ['%s<%s%s href="%s"/>' % (indent, self.moniker, attributeStr, self.build_href(**kwargs))]
        else:  # GNDS-2.0 or later
            attrs = ""
            for attr in ('ENDF_MFMT', 'dimension'):
                if getattr(self, attr) is not None:
                    attrs += ' %s="%s"' % (attr, getattr(self, attr))
            closeTag = ">" if self.slices else "/>"

            XMLList = ['%s<%s%s href="%s"%s' % (indent, self.moniker, attrs, self.build_href(**kwargs), closeTag)]
            if self.slices:
                XMLList += self.slices.toXML_strList(indent2, **kwargs)
                XMLList[-1] += "</%s>" % self.moniker

        return XMLList

    @classmethod
    def parseNodeUsingClass(cls, linkElement, xPath, linkData, **kwargs):

        attrs1_10 = None
        if (linkElement.get('L') or linkElement.get('domainMin')
                and linkData['formatVersion'] == GNDS_formatVersionModule.version_1_10):
            attrs1_10 = {attr: linkElement.attrib.pop(attr) for attr in
                         ('domainMin', 'domainMax', 'L', 'domainUnit') if attr in linkElement.attrib}
            if 'L' in attrs1_10:  # angular distribution covariance
                attrs1_10['domainValue'] = attrs1_10['L']
                attrs1_10['dimension'] = 1
                attrs1_10['outerDimension'] = 2
                del attrs1_10['L']
            else:  # energy distribution covariance
                attrs1_10['dimension'] = 2
                attrs1_10['outerDimension'] = 1

        instance = super(DataLink, cls).parseNodeUsingClass(linkElement, xPath, linkData, **kwargs)
        if linkElement.find(Slices.moniker):
            instance.slices.parseNode(linkElement.find(Slices.moniker), xPath, linkData, **kwargs)
        elif attrs1_10:
            instance.dimension = attrs1_10.pop('outerDimension')
            instance.slices.add(Slice(**attrs1_10))
        return instance


class RowData(DataLink):
    moniker = 'rowData'


class ColumnData(DataLink):
    moniker = 'columnData'


class Slice(ancestryModule.AncestryIO):
    """
    Used inside covariances for multi-dimensional functions.
    Each Slice fixes a value or range along one dimension of the multi-dimensional function.
    """

    moniker = "slice"

    def __init__(self, dimension: int, domainUnit=None, domainMin=None, domainMax=None,
                 domainValue=None):

        if domainValue is not None:
            assert domainMin is None and domainMax is None, "domainValue must not be supplied along with domainMin/Max"
        else:
            assert None not in (domainMin, domainMax), "domainMin and domainMax must be supplied together"

        super().__init__()
        self.__dimension = int(dimension)
        self.__domainUnit = domainUnit
        self.__domainMin = float(domainMin) if domainMin is not None else None
        self.__domainMax = float(domainMax) if domainMax is not None else None
        if domainValue is not None:
            domainValue = float(domainValue)
            if domainValue.is_integer(): domainValue = int(domainValue)
        self.__domainValue = domainValue

    @property
    def dimension(self):
        return self.__dimension

    @property
    def label(self):
        return str(self.__dimension)

    @property
    def domainUnit(self):
        return self.__domainUnit

    @property
    def domainMin(self):
        return self.__domainMin

    @property
    def domainMax(self):
        return self.__domainMax

    @property
    def domainValue(self):
        return self.__domainValue

    def convertUnits(self, unitMap):
        if self.__domainUnit in unitMap:
            newUnit = unitMap[self.__domainUnit]
            from pqu import PQU
            factor = PQU.PQU(1, self.__domainUnit).getValueAs(newUnit)

            if self.__domainMin is not None:
                self.__domainMin *= factor
            if self.__domainMax is not None:
                self.__domainMax *= factor
            if self.__domainValue is not None:
                self.__domainValue *= factor
            self.__domainUnit = newUnit

    def toXML_strList(self, indent='', **kwargs):

        attributesStr = ""
        for attr in ('domainMin', 'domainMax', 'domainValue', 'domainUnit'):
            if getattr(self, attr) is not None:
                attributesStr += ' %s="%s"' % (attr, getattr(self, attr))
        xmlStringList = ['%s<%s dimension="%d"%s/>' % (indent, self.moniker, self.__dimension, attributesStr)]
        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)
        slice_ = cls(node.get('dimension'), node.get('domainUnit'),
                     node.get('domainMin'), node.get('domainMax'), node.get('domainValue'))
        xPath.pop()
        return slice_


class Slices(suitesModule.Suite):
    moniker = "slices"

    def __init__(self):
        suitesModule.Suite.__init__(self, [Slice])
