# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from xData import enums as enumsModule
from . import values as valuesModule
from . import base as baseModule
from . import XYs1d as XYs1dModule

"""
This module containes all the classes for handling a GNDS xs_pdf_cdf1d data container and its children.

This module contains the following classes:

    +-------------------+---------------------------------------------------------------------------------------+
    | Class             | Description                                                                           |
    +===================+=======================================================================================+
    | Data              | This is the base class for the classes :py:class:`Xs`, :py:class:`Pdf` and            |
    |                   | :py:class:`Cdf`.                                                                      |
    +-------------------+---------------------------------------------------------------------------------------+
    | Xs                | This class stores the domain values of the pdf.                                       |
    +-------------------+---------------------------------------------------------------------------------------+
    | Pdf               | This class stores the probability density function (pdf) values.                      |
    +-------------------+---------------------------------------------------------------------------------------+
    | Cdf               | This class stores the cumulative distribution function (cdf) values of the pdf.       |
    +-------------------+---------------------------------------------------------------------------------------+
    | Xs_pdf_cdf1d      | This class represents a GNDS xs_pdf_cdf1d data container.                             |
    +-------------------+---------------------------------------------------------------------------------------+
"""


class Data(ancestryModule.AncestryIO):
    """
    This is the base class for the classes :py:class:`Xs`, :py:class:`Pdf` and :py:class:`Cdf`.

    The following table list the primary members of this class:

    +-----------+-------------------------------------------------------------------+
    | Member    | Description                                                       |
    +===========+===================================================================+
    | values    | This is the list of values (i.e., data) for the derived class.    |
    +-----------+-------------------------------------------------------------------+
    """

    def __init__(self, values):
        """
        :param values:      This is the list of values (i.e., data) for the derived class.
        """

        if not isinstance(values, valuesModule.Values): raise TypeError('Not a values type')
        ancestryModule.AncestryIO.__init__(self)

        self.values = values

    def __len__(self):
        """
        This method returns the number of values of the domain which is the same as for the pdf and cdf.

        :returns:       A python int.
        """

        return len(self.values)

    def copy(self):
        """
        This method returns a copy of *self*.

        :returns:       An instance of self.
        """

        return self.__class__(self.values)

    def offsetScaleData(self, offset, scale):
        """
        This method modifies every value in *self* by multiply it by *scale* and adding *offset*.

        :param offset:      The offset to apply to each value.
        :param scale:       The multiplicative scale factor to apply to each value.
        """

        self.values.offsetScaleValues(offset, scale)

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        XMLStringList = '%s<%s>' % (indent, self.moniker)
        XMLStringList += ''.join(self.values.toXML_strList('', **kwargs))
        XMLStringList += '</%s>' % self.moniker

        return [XMLStringList]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        data1 = cls(valuesModule.Values.parseNodeUsingClass(node[0], xPath, linkData, **kwargs))

        xPath.pop()

        return data1


class Xs(Data):
    """
    This class stores the domain values of the pdf. 
    """

    moniker = 'xs'


class Pdf(Data):
    """
    This class stores the probability density function (pdf) values.
    """

    moniker = 'pdf'


class Cdf(Data):
    """
    This class stores the cumulative distribution function (cdf) value of the pdf. 
    """

    moniker = 'cdf'


class Xs_pdf_cdf1d(baseModule.XDataFunctional):
    """
    This class represents a GNDS xs_pdf_cdf1d data container.

    The following table list the primary members of this class:

    +-------------------+-----------------------------------------------------------------------+
    | Member            | Description                                                           |
    +===========+===============================================================================+
    | xs                | This is the domain values of the pdf.                                 |
    +-------------------+-----------------------------------------------------------------------+
    | pdf               | This is the probability density function (pdf) values.                |
    +-------------------+-----------------------------------------------------------------------+
    | cdf               | This is the cumulative distribution function (cdf) values of the pdf. |
    +-------------------+-----------------------------------------------------------------------+
    | interpolation     | This is the interpolation for the xs and pdf.                         |
    +-------------------+-----------------------------------------------------------------------+
    | axes              | This is the axes member.                                              |
    +-------------------+-----------------------------------------------------------------------+
    | outerDomainValue  | This is the domain value for the next higher dimension for            |
    |                   | a function that is embedded in a high dimensional functions.          |
    +-------------------+-----------------------------------------------------------------------+
    | index             | This is the index member use by some xData classes.                   |
    +-------------------+-----------------------------------------------------------------------+
    | label             | This is the label member use by some xData classes.                   |
    +-------------------+-----------------------------------------------------------------------+

    """

    moniker = 'xs_pdf_cdf1d'
    dimension = 1
    ancestryMembers = ('xs', 'pdf', 'cdf')

    def __init__(self, _xs, _pdf, _cdf, axes=None, label=None, index=None, outerDomainValue=None,
                 interpolation=enumsModule.Interpolation.linlin):
        """
        :param xs:                  This is the domain values of the pdf.
        :param pdf:                 This is the probability density function (pdf) values.
        :param cdf:                 This is the cumulative distribution function (cdf) values of the pdf.
        :param interpolation:       This is the interpolation for the xs and pdf.
        :param axes:                This is the axes member. 
        :param outerDomainValue:    This is the domain value for the next higher dimension for a function that is embedded in a high dimensional functions.
        :param index:               This is the index member use by some xData classes.
        :param label:               This is the label member use by some xData classes.
        """

        baseModule.XDataFunctional.__init__(self, axes=axes, label=label, index=index,
                                            outerDomainValue=outerDomainValue)

        self.interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)

        if not isinstance(_xs, Xs): raise TypeError('Invalid Xs: "%s"' % type(_xs))
        self.__xs = _xs
        self.__xs.setAncestor(self)

        if not isinstance(_pdf, Pdf): raise TypeError('Invalid pdf: "%s"' % type(_pdf))
        self.__pdf = _pdf
        self.__pdf.setAncestor(self)

        if not isinstance(_cdf, Cdf): raise TypeError('Invalid cdf: "%s"' % type(_cdf))
        self.__cdf = _cdf
        self.__cdf.setAncestor(self)

        if len(_xs) != len(_pdf) != len(_cdf):
            raise ValueError('lengths not the same: %s %s %s' % (len(_xs), len(_pdf), len(_cdf)))

    def __len__(self):
        """
        This method returns the number of values of the domain which is the same as for the pdf and cdf.

        :returns:       A python int.
        """

        return len(self.__xs)

    @property
    def xs(self):
        """
        This method returns a reference to *self* xs child.
        """

        return self.__xs

    @property
    def pdf(self):
        """
        This method returns a reference to *self* pdf child.
        """

        return self.__pdf

    @property
    def cdf(self):
        """
        This method returns a reference to *self* cdf child.
        """

        return self.__cdf

    def asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):
        """
        This method returns a representation of the data in *self* as an :py:class:`XYs1dModule.XYs1d` instance. 

        :param asLinLin:    If **True**, the data have lin-lin interpolation.
        :param accuracy:    Used to determine the accuracy if converting data to lin-lin interpolated data.
        :param lowerEps     Used to dull the lower point for "flat" interpolation.
        :param upperEps     Used to dull the upper point for "flat" interpolation.

        :returns:           A :py:class:`XYs1dModule.XYs1d` instance.
        """

        xys1d = self.toPointwise_withLinearXYs(accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps)

        return xys1d

    def convertUnits(self, unitMap):
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        axes = self.axes
        if axes is None: axes = self.ancestor.axes.copy()  # BRB FIXME, the prior line should have worked.
        factors = axes.convertUnits(unitMap)
        if factors[0] != 1.0: self.pdf.offsetScaleData(0.0, factors[0])
        if factors[1] != 1.0: self.xs.offsetScaleData(0.0, factors[1])
        self.fixValuePerUnitChange(factors)

    def copy(self):
        """
        This method returns a copy of *self*.

        :returns:       An instance of self.
        """

        return self.__class__(self.xs.copy(), self.pdf.copy(), self.cdf.copy(),
                              outerDomainValue=self.outerDomainValue, axes=self.axes,
                              interpolation=self.interpolation)

    def toPointwise_withLinearXYs(self, **kwargs):
        """
        This method returns an :py:class:`XYs1dModule.XYs1d` representation of *self*'s pdf.

        :param kwargs:      Not used but present to be compatible with other similar methods.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        _pdf, _cdf = self.to_pdf_and_cdf()
        return _pdf.toPointwise_withLinearXYs(**kwargs)

    def to_pdf_and_cdf(self):
        """
        This method returns two :py:class:`XYs1dModule.XYs1d` representation of *self*. One representing the pdf and one 
        representing the cdf.  The interpolation is the same as self's interpolation.  Ergo, the interpolation for the cdf 
        is not correct as, in general, the correct cdf interpolation is not defined in GNDS. The units for the cdf will be 
        correct but not the label.

        :returns:           Tow instances of :py:class:`XYs1dModule.XYs1d`.
        """

        _pdf = []
        _cdf = []
        for index, x1 in enumerate(self.__xs.values):
            _pdf.append([x1, self.__pdf.values[index]])
            _cdf.append([x1, self.__cdf.values[index]])

        linear = self.toLinearXYsClass()

        axes = self.ancestor.axes
        _pdf = linear(data=_pdf, axes=axes, interpolation=self.interpolation, outerDomainValue=self.outerDomainValue)
        if axes is not None:
            axes = axes.copy()
            axes[0].unit = axes[0].multiplyUnit(axes[1])
        _cdf = linear(data=_cdf, axes=axes, interpolation=self.interpolation, outerDomainValue=self.outerDomainValue)

        return _pdf, _cdf

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        xs_pdf_cdf1d_singleLine = kwargs.get('xs_pdf_cdf1d_singleLine', True)
        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent
        if xs_pdf_cdf1d_singleLine: indent2 = ''

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        interpolationStr = ''
        if self.interpolation != enumsModule.Interpolation.linlin:
            interpolationStr = ' interpolation="%s"' % self.interpolation
        XML_strList = ['%s<%s%s%s>' % (indent, self.moniker, attributeStr, interpolationStr)]
        if self.isPrimaryXData():
            if self.axes is not None: XML_strList += self.axes.toXML_strList(indent=indent2, **kwargs)
        XML_strList += self.xs.toXML_strList(indent2, **kwargs)
        XML_strList += self.pdf.toXML_strList(indent2, **kwargs)
        XML_strList += self.cdf.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        if xs_pdf_cdf1d_singleLine: XML_strList = [''.join(XML_strList)]

        return XML_strList

    def toLinearXYsClass(self):
        """
        This method always returns the class :py:class:`XYs1dModule.XYs1d`.

        :returns:       The class :py:class:`XYs1dModule.XYs1d`.
        """

        return XYs1dModule.XYs1d

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(
            node, xPath, allowInterapolation=True)  # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0:
            raise Exception('Invalid attributes: %s.' % (', '.join(list(extraAttributes.keys()))))

        datas = {}
        children = {'xs': Xs, 'pdf': Pdf, 'cdf': Cdf}
        for child in node:
            if child.tag in children:
                datas[child.tag] = children[child.tag].parseNodeUsingClass(child, xPath, linkData, **kwargs)

        axes = kwargs.get('axes', None)
        if axes is not None: attributes['axes'] = axes

        instance = cls(datas['xs'], datas['pdf'], datas['cdf'], **attributes)

        extraNodesList = baseModule.XDataFunctional.parseNodeStandardChildren(instance, node, xPath, linkData, **kwargs)
        extraNodes = {}
        for extraNode in extraNodesList: extraNodes[extraNode.tag] = extraNode
        extraNodes.pop('xs')
        extraNodes.pop('pdf')
        extraNodes.pop('cdf')

        if len(extraNodes) > 0:
            raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()

        return instance

    @classmethod
    def fromXYs(cls, xys, outerDomainValue=None, thinEpsilon=None):
        """
        This method creates an :py:class:`Xs_pdf_cdf1d` instance for an :py:class:`XYs1dModule.XYs1d` instance.
        If *thinEpsilon* is not None, then points are thinned so that the returned cdf has cdf[i+1] - cdf[i] > *thinEpsilon*.

        :param cls:                 The :py:class:`Xs_pdf_cdf1d` class of the returned instance.
        :param xys:                 The :py:class:`XYs1dModule.XYs1d` to convert to an `Xs_pdf_cdf1d` instance.
        :param outerDomainValue:    The value of the *outerDomainValue* for the returned instance. If None, taken from *xys*.
        :param thinEpsilon:         The thinning parameter.
        """

        if xys.interpolation not in (enumsModule.Interpolation.linlin, enumsModule.Interpolation.flat):
            xys = xys.toPointwise_withLinearXYs(accuracy=1e-3, lowerEps=0, upperEps=1e-8)

        if outerDomainValue is None:
            outerDomainValue = xys.outerDomainValue

        try:
            norm = xys.normalize()
        except:
            norm = xys.copy()
            for i1 in range(len(norm)):
                norm[i1] = [norm[i1][0], 1.0]
            if len(norm) > 2:
                norm[0] = [norm[0][0], 0.0]
            norm[-1] = [norm[-1][0], 0.0]
            norm = norm.normalize()

        xs, pdf = norm.copyDataToXsAndYs()
        cdf = norm.runningIntegral()
        cdf[-1] = 1.0

        if thinEpsilon is not None and xys.interpolation is not enumsModule.Interpolation.flat:
            indicesToRemove = []
            for index, cdf2 in enumerate(cdf):
                if cdf2 == 1:
                    break
                if index != 0:
                    if cdf2 - cdf1 <= thinEpsilon * cdf2:
                        indicesToRemove.append(index)
                    else:
                        lastKept = index
                cdf1 = cdf2
            for index in range(index, len(cdf) - 1):
                indicesToRemove.append(index)
            if len(indicesToRemove) > 0:
                for index in reversed(indicesToRemove):
                    xs.pop(index)
                    pdf.pop(index)
                    cdf.pop(index)

        xs = Xs(valuesModule.Values(xs))
        pdf = Pdf(valuesModule.Values(pdf))
        cdf = Cdf(valuesModule.Values(cdf))

        return cls(xs, pdf, cdf, outerDomainValue=outerDomainValue, interpolation=xys.interpolation)
