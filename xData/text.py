# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains all the classes for handling GNDS text container.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Encoding                          | This class is an enum of the supported text encodings.                |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Markup                            | This class is an enum of the supported text markups.                  |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Text                              | This class represents a GNDS text node.                               |
    +-----------------------------------+-----------------------------------------------------------------------+

This module contains the following functions:

    +-----------------------------------+---------------------------------------------------------------------------+
    | Class                             | Description                                                               |
    +===================================+===========================================================================+
    | isString                          | This function test if an object is a valid python string.                 |
    +-----------------------------------+---------------------------------------------------------------------------+
    | raiseIfNotString                  | This function does a raise if its argument is not a valid python string.  |
    +-----------------------------------+---------------------------------------------------------------------------+
"""

from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule


def isString(string):
    """
    This function executes a TypeError raise if *string* is not a valid python2 str or unicode instance,
    or if is not a valid python3 str instance.
    This function returns the *string* argument unchanged.

    :param string:  The python instance to check.

    :returns:       The argument *string*.
    """

    if not isinstance(string, str): raise TypeError('Text must be a str')

    return string


class Encoding(enumsModule.Enum):
    """
    This class is an enum of the supported text encodings.
    """

    ascii = enumsModule.auto()
    utf8 = enumsModule.auto()


class Markup(enumsModule.Enum):
    """
    This class is an enum of the supported text markups.
    """

    none = enumsModule.auto()
    xml = enumsModule.auto()
    html = enumsModule.auto()
    latex = enumsModule.auto()


class Text(ancestryModule.AncestryIO):
    """
    This class represents a GNDS text container. GNDS text is a string of characters with an encoding and markup.
    Encoding describes how bits of the text are to be interpreted (e.g., ascii, utf8). Markup describes
    what character sequences have special meaning (e.g., xml, latex).

    If the body is set to None, the Text instance is considered to be empty and its toXML_strList returns an empty list.

    The following table list the primary members of this class:

    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | label     | The label for the text.                                       |
    +-----------+---------------------------------------------------------------+
    | encoding  | The type of encoding for the text.                            |
    +-----------+---------------------------------------------------------------+
    | markup    | The markup for the text.                                      |
    +-----------+---------------------------------------------------------------+
    | body      | The characters representing the text.                         |
    +-----------+---------------------------------------------------------------+
    """

    moniker = 'text'

    # FIXME Should add allowed encodings as an argument with default of all. Like encodings = Encoding.allowed.
    #  Also should check in body setter.

    def __init__(self, text=None, encoding=Encoding.ascii, markup=Markup.none, label=None):
        """
        :param text:        The characters representing the text.
        :param encoding:    The type of encoding for the text.
        :param markup:      The markup for the text.
        :param label:       The label for the text.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.encoding = Encoding.checkEnumOrString(encoding)
        self.markup = Markup.checkEnumOrString(markup)

        if label is not None:
            if not isinstance(label, str): raise TypeError('Invalid label instance.')
        self.__label = label

        self.body = text

    def __len__(self):
        """
        This method returns the number of characters in the body of *self*.
        """

        return len(self.__body)

    def __getitem__(self, index):
        """
        This method returns the character at index in the body of *self*.

        :param index:   This index of the character to return with the body of *self*.
        """

        return self.__body[index]

    @property
    def label(self):
        """
        This method returns *self*'s label.
        """

        return self.__label

    @property
    def encoding(self):
        """
        This method returns *self*'s encoding.

        :returns:       An instance of :py:class:`Encoding`.
        """

        return self.__encoding

    @encoding.setter
    def encoding(self, value):
        """
        This method set the encoding of *self* to *value*.

        :param value:       The new encoding value.
        """

        self.__encoding = Encoding.checkEnumOrString(value)

    @property
    def markup(self):
        """
        This method returns *self*'s markup.

        :returns:       An instance of :py:class:`Markup`.
        """

        return self.__markup

    @markup.setter
    def markup(self, value):
        """
        This method set the markup of *self* to *value*.

        :param value:       The new markup value.
        """

        self.__markup = Markup.checkEnumOrString(value)

    @property
    def body(self):
        """
        This method returns *self*'s body.

        :returns:       A python str.
        """

        return self.__body

    @body.setter
    def body(self, text):
        """
        This method set the body of *self* to *text*, overriding the current text.
        If *text* is None, *self* is considered to be empty.

        :param text:    A python str instance.
        """

        self.__filled = False
        self.__body = ''

        text2 = text
        if text is None: text2 = ''

        isString(text2)

        if self.encoding == Encoding.ascii:  # If encoding is ascii verify text is.
            try:
                text2.encode('ascii')
            except UnicodeEncodeError as ex:
                print("    ERROR: encountered non-ascii character '%s' at index %d of text" %
                      (ex.object[ex.start], ex.start))
                print("    Bad character in context:\n", ex.object[ex.start - 20:ex.end + 20], "\n")
                raise ex

        if text is not None: self.__filled = True
        self.__body = text2

    @property
    def filled(self):
        """
        This method returns False if *self* is empty and True otherwise.
        """

        return self.__filled

    def copy(self):
        """
        This method returns a copy of *self*.

        :returns:           An instance of :py:class:`Text`.
        """

        # No need to copy text has it is immutable.
        return Text(self.body, self.encoding, self.markup)

    __copy__ = copy

    def numberOfLines(self):
        """
        This method returns the number of line-feeds in the body of *self*.
        """

        _numberOfLines = self.__body.count('\n')
        if len(self.__body) > 0:
            if self.__body[-1] != '\n': _numberOfLines += 1

        return _numberOfLines

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns an empty python str.
        Method is designed to be over-ridden by derived classes which add attributes not in the :py:class:`Text` class.

        :returns:       An empty python str instance.
        """

        return ''

    def bodyToXML_CDATA(self, **kwargs):
        """
        This method wraps the body of self into an XML "<![CDATA[" and "]]>" tag.
        """

        return '<![CDATA[%s]]>' % self.__body

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.
        If *self* is empty and kwargs.get('showEmptyText') is False, this method returns an empty list.

        :param indent:     The minimum amount of indentation.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:           List of str instances representing the XML lines of self.
        """

        kwargs['addExtraAttributes'] = self.__filled
        extraAttributes = self.XML_extraAttributes(**kwargs)

        if not self.__filled:
            showEmptyText = kwargs.get('showEmpty', False) or kwargs.get('showEmptyText', False)
            if extraAttributes == '' and not showEmptyText:
                return []
            XML_string = '<%s%s/>' % (self.moniker, extraAttributes)
            if extraAttributes == '':
                XML_string = '%-24s <!-- text -->' % XML_string
            return [indent + XML_string]

        attributeStr = ''
        if self.__label is not None: attributeStr += ' label="%s"' % self.__label
        if self.__encoding != Encoding.ascii: attributeStr += ' encoding="%s"' % self.__encoding
        if self.__markup != Markup.none: attributeStr += ' markup="%s"' % self.__markup
        attributeStr += extraAttributes

        return ['%s<%s%s>%s</%s>' % (indent, self.moniker, attributeStr, self.bodyToXML_CDATA(**kwargs), self.moniker)]

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method sets data in *self* using the contents of *node*.

        :param node:       Node to parse.
        :param xPath:      List containing xPath to current node, useful mostly for debugging.
        :param linkData:   dict that collects unresolved links.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.
        """

        xPath.append(node.tag)

        self.encoding = node.get('encoding', Encoding.ascii)
        self.markup = node.get('markup', Markup.none)

        text = node.find('text')
        if text is None:  # Format >= '2.0.LLNL_4' parsing.
            if node.text is None:
                self.body = None
            elif len(node.text) == 0:
                self.body = None
            else:
                self.body = node.text
        else:  # Pre format '2.0.LLNL_4' parsing where the actual text node is a child node named 'text'.
            self.parseNode(text, xPath, linkData, **kwargs)

        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:        Form class to return.
        :param node:       Node to parse.
        :param xPath:      List containing xPath to current node, useful mostly for debugging.
        :param linkData:   dict that collects unresolved links.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :returns:          An instance of *cls* representing *node*.
        """

        instance = cls()
        instance.parseNode(node, xPath, linkData, **kwargs)

        return instance


def raiseIfNotString(string, name):
    """
    This function does a raise if *string* is not a valid python string.

    :param string:  The python instance to check.
    :param name:    The name of the variable representing *string*.
    """

    if not isinstance(string, str): raise TypeError('Invalid %s instance.' % name)
    return string
