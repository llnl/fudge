# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains all the ancestry classes. These classes allow a node in a hierarchy 
to reference any other node in the hierarchy provided that each node of the hierarchy inherits
from one of the ancestry classes and the method :py:func:`Ancestry.setAncestor` is called to
contect a child node to its parent node.

This module was mainly written to support linking and xpath referencing in a GNDS hierarchy.

    This module contains the following classes:
                             
    +-------------------+-----------------------------------------------------------------------+ 
    | Class             | Description                                                           |
    +===================+=======================================================================+
    | Ancestry          | Base class used by other ancestry classes.                            |
    +-------------------+-----------------------------------------------------------------------+
    | AncestryIO_base   | This class adds reading and writing methods.                          |
    +-------------------+-----------------------------------------------------------------------+
    | AncestryIO        | Base class for :py:class:`FSubform`, :py:class:`FSubform`,            |
    +-------------------+-----------------------------------------------------------------------+
    | AncestryIO_bare   | This is an abstruct class that defined the abstract                   |
    |                   | method :py:class:`AncestryIO_bare.parseNodeUsingClass`.               |
    +-------------------+-----------------------------------------------------------------------+
    | AncestryIO_base   | This is an abstruct class that defined the abstract                   |
    |                   | method :py:class:`AncestryIO_bare.parseNode`.                         |
    +-------------------+-----------------------------------------------------------------------+
"""

import sys
import os
import abc
import copy
import numpy
import pathlib
import inspect

from xml.etree import ElementTree

from LUPY import xmlNode as xmlNodeMode        # Wrapper around the xml parser.
from LUPY import checksums as checksumsModule

XML_declaration = '<?xml version="1.0" encoding="UTF-8"?>\n'

class Ancestry(abc.ABC):
    """
    This class is designed to be a base class for a class (instance) that is a member of another
    class (instance). The function of this class is to aid in tracking a class's ancestors and all
    of the ancestors members. That is, if an instance is part of a hierarchy, this class provides 
    methods, for example, that list the position of the instance within the hierarchy or give its 
    position relative to another member in the hierarchy using xlinks.

    This class defines the following members:

        +---------------------------+-----------------------------------------------------------------------------------------------------------+
        | moniker                   | The name associated with the class type.                                                                  |
        +---------------------------+-----------------------------------------------------------------------------------------------------------+
        | ancestor                  | Instance which *self* is a child of.                                                                      |
        +---------------------------+-----------------------------------------------------------------------------------------------------------+
        | keyName                   | The name of the key. The key's value is used by a suite to indentify an element of a suite.               |
        +---------------------------+-----------------------------------------------------------------------------------------------------------+
        | ancestryMembers           | Tuple of names of the members of *self* that inherit from **Ancestry**.                                   |
        +---------------------------+-----------------------------------------------------------------------------------------------------------+
        | legacyMemberNameMapping   | A map whose keys are legacy member names and whose associated values are the current memeber names.       |
        +---------------------------+-----------------------------------------------------------------------------------------------------------+
        | monikerByFormat           | If the moniker had a different name for in a prior format version, then this is a dict of key/value where |
        +---------------------------+-----------------------------------------------------------------------------------------------------------+

    This class defines the member moniker to be an abstract member, meaning that all classes that inherit from this class must define moniker.
    The main reason for the **keyName** is that is it used by a **Suite** instance to reference (i.e., uniquely identify) its children.

    An instance has the xlink '/the/list/of/ancestors/andSelf' if its attribute member is None or
    '/the/list/of/ancestors/andSelf[@attribute="value"]' if its attribute member is not None where each string
    part of the xpath the value of the moniker for that node. For example,
    for a hierarchy consisting of::

        class A with moniker 'nameA' and 
            member mB of class B, class B with moniker 'nameB' and 
                member mC of class C and class C with moniker 'nameC',

    then the xlink for an instance of a C class in the hierarchy is '/nameA/nameB/nameC'. If the mC instance
    set the attribute member to 'greeting', then for this example with value 'Hi', then the
    xlink for the C class is '/nameA/nameB/nameC[@greeting="Hi"]'.
    """

    keyName = None
    ancestryMembers = tuple()
    legacyMemberNameMapping = {}
    monikerByFormat = {}

    def __init__(self):

        self.__ancestor = None

    def __str__(self):
        """Returns the xpath of *self*."""

        return self.toXLink()

    @property
    @abc.abstractmethod
    def moniker(self):
        """Abstract method that must be defined by all classes that inherit this class. Must return the moniker of the instance."""

        pass

    @property
    def ancestor(self):
        """Returns self's ancestor. That is, its parent."""

        return self.__ancestor

    @property
    def keyValue(self):
        """Returns self's keyValue as determined by the member *keyName*. If keyName is None, then None is returned."""

        if self.keyName is None: return None
        return getattr(self, self.keyName)

    @property
    def rootAncestor(self):
        """
        Traverse up the ancestry tree to the root ancestor and returns a reference to it. 
        The root ancestor is the ancestor instance whose ancestor is None.
        """

        ancestor = self
        while(ancestor.__ancestor is not None): ancestor = ancestor.__ancestor
        return ancestor

    def checkAncestry(self, verbose = 0, level = 0):
        """
        This method checks that all ancestryMembers are properly set up. This method checks *self* and all its child nodes as
        determined by the *ancestryMembers* member.
        """

        def check(child):
            """
            For internal use only.
            """

            if child is None: return
            if self.isChild(child):
                child.checkAncestry(verbose = verbose, level = level)
            else:
                print('WARNING from checkAncestry: member "%s" not a child of %s' % (child, self.toXLink()))
                print('    Its ancestry is: %s' % child.toXLink())

        if len(self.ancestryMembers) == 0: return

        prefix = ( level + 1 ) * '    '
        for member in self.ancestryMembers:
            if member == '': continue
            if verbose > 0: print("%s%s" % (prefix, member))
            doLoop = False
            if member[0] == '[':
                member = member[1:]
                doLoop = True
            if hasattr(self, member):
                m1 = getattr(self, member)
                if doLoop:
                    for child in m1: check(child)
                else:
                    check(m1)
            else:
                print('WARNING from checkAncestry: %s does not have member "%s"' % (self.toXLink(), member))

    def copy(self):
        """
        Uses deepcopy to make a copy of self.

        :return:            A copy of self.
        """

        from xData import link as linkModule

        theCopy = copy.deepcopy(self)
        theCopy._copyPost()

        for link in theCopy.findInstancesOfClassInChildren(linkModule.Link) + theCopy.findInstancesOfClassInChildren(linkModule.Link2): link.updateLink()

        return theCopy

    def _copyPost(self):
        """
        For internal use only. Recursively fixes up the ancestor stuff in self.
        """

        from PoPs import suite as PoPsSuiteModule
        from fudge import suites as suitesModule

        for memberName in self.ancestryMembers:
            name = memberName
            if name[0] == '[': name = name[1:]
            member = getattr(self, name)
            if member is None: continue
            if not isinstance(member, list):
                member.setAncestor(self)
                member._copyPost()
            if isinstance(member, (PoPsSuiteModule.Suite, suitesModule.Suite, list)):
                for child in member:
                    if isinstance(member, list):
                        child.setAncestor(self)
                    else:
                        child.setAncestor(member)
                    child._copyPost()

    def findAttributeInAncestry(self, attributeName):
        """
        This method checks if *self* has an member (i.e., attribute) with name *attributeName*. If so, it returns its value. If
        not, this method recursively searches up the ancestry tree by calling :py:class:`Ancestry.findAttributeInAncestry` on its ancestor.
        If root is reached and no member named *attributeName* is found, a raise is executed.

        :param attributeName:   The name of the member whose value it to be returned.

        :raise:                 If member with name *attributeName* is not found in self or its ancestry.
        """

        if hasattr(self, attributeName): return getattr(self, attributeName)
        if self.__ancestor is None: raise Exception('Could not find attribute name = %s in ancestry' % attributeName)
        return self.__ancestor.findAttributeInAncestry(attributeName)

    def findClassInAncestry(self, class_):
        """
        This method checks if *self* has an member that is an instance of *class_*. If so, a reference to it is returned.
        If not, this method recursively searches up the ancestry tree by calling :py:class:`Ancestry.findClassInAncestry` on its ancestor
        until an instance of *class_* is found or root is reached. If root is reached and no instance of *class_* found,
        a raise is executed.

        :param class_:      The class to search for in self and its ancestors.

        :raises:            If an instance of class *class_* is not found in self or its ancestry.
        """

        if isinstance(self, class_): return self
        if self.__ancestor is None: raise Exception('Could not find class name = %s in ancestry' % class_.__name__)
        return self.__ancestor.findClassInAncestry(class_)

    def findEntity(self, entityName, attribute = None, value = None):
        """
        This method takes one part of an xpath (e.g., the 'pathB' in /pathA/pathB/pathC') and returns a 
        reference to the instance of *self* that match that part of the xpath.
        Thie is a default findEntity method. In general, sub-classes should over-ride this method.
        This method uses the following algorithm to find entity. Firstly, if 'attribute' is None, then self is assumed
        to have an attribute named entityName which is taken to be the desired entity. Otherwise, self is iterated
        over until an item with an attribute named attribute with value *value* is found. In either case, if
        an entity is found, its moniker value must be entityName. If no entity is found, a raise AttributeError is executed.

        :param entityName:      ?
        :param attribute:       ?
        :param value:           ?
        """

        if entityName in self.legacyMemberNameMapping: entityName = self.legacyMemberNameMapping[entityName]

        if entityName in ('.', self.moniker):
            return self
        elif entityName == '..':
            return self.__ancestor

        entity = None
        if attribute is None:
            entity = getattr(self, entityName)
        else:
            try:                       # try needed in case self cannot be iterated.
                for entity1 in iter(self):
                    if str(getattr(entity1, attribute, None)) == value:
                        entity = entity1
                        break
            except TypeError:
                try:
                    entity = getattr(self, entityName)
                except AttributeError:
                    pass
        if entity is None or entityName != getattr(entity, 'moniker'):
            raise AttributeError("Can't find entity %s in %s" % (entityName,self))

        return entity

    def findInstancesOfClassInChildren(self, cls, level = 9999):
        """
        Finds all instances of class *cls* in self's children, grand-children, etc.
        """

        foundInstances = []
        level -= 1
        if level < 0: return foundInstances
        for ancestryMember in self.ancestryMembers:
            instance = getattr(self, ancestryMember)
            if instance is None: continue
            if isinstance(instance, cls): foundInstances.append(instance)
            if isinstance(instance, list):
                for item in instance: foundInstances += item.findInstancesOfClassInChildren(cls, level = level)
            else:
                foundInstances += instance.findInstancesOfClassInChildren(cls, level = level)

        return foundInstances

    def isChild(self, child):

        if isinstance(child, Ancestry): return child.__ancestor == self

        return False

    def isParent(self, parent):

        return self.__ancestor == parent

    def setAncestor(self, ancestor):
        """Sets self's ancestor to ancestor."""

        self.__ancestor = ancestor 

    def toRelativeXLink(self, other = None, formatVersion = None):
        """
        Returns a string that is a relative xlink to another element (using XML xpath syntax).
        Both elements must reside in the same hierarchy.  For a description of xpath, see 
        http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics.
        """

        if other is None:
            if self.__ancestor is None: return ''
            return self.__ancestor.toRelativeXLink() + '../'
        else:
            if self.rootAncestor is not other.rootAncestor:
                raise Exception('Root ancestors not the same ("%s" != "%s")' % (self.toXLink(), other.toXLink()))
            thisPath = self.toXLink(formatVersion = formatVersion).split('/')
            othersPath = other.toXLink(formatVersion = formatVersion).split('/')
            for i1, tag in enumerate(thisPath):
                if i1 >= len(othersPath): break
                if tag != othersPath[i1]: break
            relativePath = ''
            for i2 in range(len(thisPath) - i1): relativePath += '../'
            relativePath += '/'.join(othersPath[i1:])
            return relativePath

    def toXLink(self, formatVersion = None):
        """
        Returns a string that is an xlink to self (using XML xpath syntax).  The resulting 
        xlink starts at the root element. For a description of xpath, see 
        http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics.
        """

        ancestorXLink = ''
        if self.__ancestor is not None: ancestorXLink = self.__ancestor.toXLink(formatVersion = formatVersion)

        attribute = ''
        keyValue = self.keyValue
        if keyValue is not None: attribute = "[@%s='%s']" % (self.keyName, keyValue)

        moniker1 = self.monikerByFormat.get(formatVersion, self.moniker)

        return ancestorXLink + '/%s%s' % ( moniker1, attribute )

    def followXPath(self, xPath):
        """
        This method returns the instance referenced to by the xpath *xPath* starting at *self* if the xpath is relative.
        This method uses Ancestry.findEntity to find each part of the xpath.

        :param xPath:   String xPath, e.g., "/reactionSuite/reactions/reaction[@label='2']".

        :return:        Class instance pointed to by xPath
        """

        def follow2(xPathList, node):
            """
            For internal use. Recursive helper function that descends the xpath starting at *node* to find the instance referenced by *xPathList*.

            :param xPathList:       Python list of xpath parts.
            :param node:            The current node to search for the first element of *xPathList*.

            :returns:               The instance referenced by *xPathList*.
            """

            if len(xPathList) == 0: return node

            xPathNext = xPathList[0]

            try:
                if "[@" in xPathNext:
                    r1,r2 = xPathNext.split("[@",1)
                    r2,r3 = r2.split("=",1)
                    r3 = r3.rstrip(']')[1:-1]
                    nodeNext = node.findEntity(r1, r2, r3)
                else:
                    nodeNext = node.findEntity(xPathNext)
            except:
                raise XPathNotFound()

            return follow2(xPathList[1:], nodeNext)

        xPathList = xPath.split('/')                        # FIXME refactor to use xml.etree.ElementPath.xpath_tokenizer?

        while xPathList[0] == '':
            xPathList = xPathList[1:]                       # Trim empty sections from the beginning.

        while xPathList[-1] == '':
            xPathList = xPathList[:-1]                      # Trim empty sections from the end as can happend with xPath = '../'.

        if '{' in xPath:                                    # Careful, qualifiers may contain '/'.
            xpl2 = []
            for val in xPathList:
                if '}' in val and '{' not in val and '{' in xpl2[-1] and '}' not in xpl2[-1]:
                    xpl2[-1] += '/' + val
                else:
                    xpl2.append(val)
            xPathList = xpl2

        try:
            return follow2(xPathList, self)
        except XPathNotFound:
            raise XPathNotFound("Cannot locate path '%s'" % xPath)

    def printMembers(self, cls2=None, exclude=False, width=3):
        """
        This method prints the type of *self* and then prints all the member added since the base class *cls2*. That is, 
        those members which are not defined in by base class *cls2*.

        :param cls2:        Must be a base class for *self**.
        :param exclude:     If True, any member starting with an underscore (i.e., '_') is also not printed.
        :param width:       The number of members printed per line.
        """

        if cls2 is None:
            cls2 = Ancestry

        if not issubclass(self.__class__, cls2):
            raise TypeError('self is not a derived fom class cls2')

        print(type(self))
        for index, cls in enumerate(inspect.getmro(self.__class__)):
            if index == 0:
                continue
            print('  %s' % cls)

        members = []
        _ancestryMembers = dir(cls2)
        for item in dir(self):
            if item[:1] == '_' and exclude:
                continue
            if item in _ancestryMembers:
                continue
            members.append(item)

        width = max(1, width)
        if len(members) > 0:
            fmt = '%%-%ds' % max(map(len, members))
        counter = 0
        sep = '   '
        print()
        for member in members:
            print(sep, fmt % member, end='')
            sep = ''
            counter += 1
            if counter % width == 0:
                print()
                sep = '   '
        if counter % width != 0:
            print()

    def printNewMembers(self, exclude=False, width=3):
        """
        This method is deprecated, please use the method **printMembers** instead.
        """

        self.printMembers(Ancestry, exclude=exclude, width=width)

    def printNonAncestryMembers(self, exclude=False, width=3):
        """
        This method is deprecated, please use the method **printMembers** instead.
        """

        self.printMembers(Ancestry, exclude=exclude, width=width)

class AncestryIO_base(Ancestry):
    """This class adds methods to read and write *self* to a file. Currently, its supports reading and writing to an XML file."""

    def toXML(self, indent = '', **kwargs):
        """
        Calls self.toXML_strList and joins its returned list with '\n'.

        :param indent:              Amount of starting indentation.
        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :returns:                   XML representation of *self*.
        """

        return '\n'.join(self.toXML_strList(indent=indent, **kwargs))

    @abc.abstractmethod
    def toXML_strList(self, indent = '', **kwargs):
        """
        This methods must be overwritten by the derived class. It must return a Python list of strings that are the XML representation of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.
        """

        pass

    def saveToOpenedFile(self, fOut, **kwargs):
        """
        Writes an XML representation of *self* to the opened file *fOut*. This method uses *self.toXML* to construct XML representation of *self*.

        :param fOut:                File stream to write *self* to.
        :param kwargs:              A dictionary that contains data to control the way this method acts.
        """

        fOut.write(self.toXML(**kwargs))
        fOut.write('\n')

    def saveToFile(self, fileName, **kwargs):
        """
        Writes an XML representation of *self* to path *fileName*. This method opens the file and calls **saveToOpenedFile** 
        and then closes the file.

        :param fileName:            File path to write *self* to.
        :param kwargs:              A dictionary that contains data to control the way this method acts.
        """

        dirname = os.path.dirname(fileName)
        if len(dirname) > 0 and not os.path.exists(dirname): os.makedirs(dirname, exist_ok=True)
        with open(fileName, "w") as fout :
            fout.write(XML_declaration)
            self.saveToOpenedFile(fout, **kwargs)

    def saveToHybrid(self, XML_name, HDF_name=None, HDF_subDir=None, minLength=3, flatten=True, compress=False, **kwargs):
        """
        Saves *self* to as two files, with most of the hierarchy in an XML file but most of the numerical data saved in an associated HDF file.
        To construct the file names, the best is to leave *HDF_name* as None. However, if the HDF file should be put into a sub-directory
        below where the XML file will be placed, then one can leave *HDF_name* as None and set *HDF_subDir*. For example, to
        place the HDF file is the sub-directory 'HDF', set HDF_subDir='HDF'. The following example will write the XML to the file
        'path/to/n-Pu239.xml' and the HDF data to the file 'path/to/HDF/n-Pu239.h5':

        node.saveToHybrid('path/to/n-Pu239.xml', HDF_subDir='HDF')

        :param XML_name:    Name of the XML file written.
        :param HDF_name:    Name of the HDF file written. If None, will be the same as XML_name with extension changed to .h5
        :param HDF_subDir:  If not None, this name is added as a directory to the front of the HDF_name.
        :param minLength:   Minimum number of numerical values in a *Values* node before writing the *Values* node to HDF.
        :param flatten:     If True, GNDS datasets are concatenated into flattened HDF5 datasets.
        :param compress:    Enable gzip + shuffle compression for HDF5 datasets.
        :param kwargs:      A dictionary that contains data to control the way this method acts.
        """

        from fudge import externalFile as externalFileModule
        from .hdf5 import HDF5_present, h5py
        if not HDF5_present:
            raise ImportError("Missing module 'h5py' must be installed to support generating HDF5 files.")

        if not hasattr(self, 'externalFiles'):
            self.saveToFile(XML_name,  **kwargs)
            return

        if HDF_name is None:
            if XML_name.endswith('.xml'):
                HDF_name = XML_name.replace('.xml','.h5')
            else:
                HDF_name = XML_name + '.h5'

        relativeHDF_name = HDF_name
        commonPath = os.path.commonpath([ XML_name, HDF_name ])
        if len(commonPath) > 1:
            relativeHDF_name = HDF_name[len(commonPath)+1:]
            if HDF_subDir is not None: relativeHDF_name = os.path.join(HDF_subDir, relativeHDF_name)

        dirname = os.path.dirname(XML_name)
        if len(dirname) > 0 and not os.path.exists(dirname): os.makedirs(dirname)
        dirname = os.path.dirname(HDF_name)
        if len(dirname) > 0 and not os.path.exists(dirname): os.makedirs(dirname)

        HDF_opts = { 'index': 0, 'minLength': minLength, 'flatten': flatten, 'iData': [], 'dData': [], 'compression': {} }
        if compress: HDF_opts['compression'] = { 'compression': 'gzip', 'shuffle': 'true', 'chunks': True}
        with h5py.File(HDF_name, "w") as h5:
            HDF_opts['h5file'] = h5
            self.externalFiles.add(externalFileModule.ExternalFile("HDF", relativeHDF_name, checksum = "deadbeef", algorithm="sha1"))

            xmlString = self.toXML(HDF_opts=HDF_opts, **kwargs)

            if len(HDF_opts['iData']) > 0:
                iData = numpy.array(HDF_opts['iData'], dtype=numpy.int32)
                h5.create_dataset('iData', data=iData, **HDF_opts['compression'])
            if len(HDF_opts['dData']) > 0:
                dData = numpy.array(HDF_opts['dData'], dtype=numpy.float64)
                h5.create_dataset('dData', data=dData, **HDF_opts['compression'])

        sha1sum = checksumsModule.Sha1sum.from_file(HDF_name)
        for idx, line in enumerate(xmlString):
            if 'externalFile' in line and 'checksum="deadbeef"' in line:
                xmlString[idx] = line.replace("deadbeef", sha1sum)
                break

        with open(XML_name, "w") as fout :
            fout.write(XML_declaration)
            fout.write(xmlString)

        self.externalFiles.pop("HDF")

    def parseAncestryMembers(self, node, xPath, linkData, **kwargs):
        """
        The method parses all members defined in the ancestryMembers of *self* that are  found in *node*.  If *kwargs* contains 
        the key *membersToSkip* then all members of *self* whose moniker is in the list *kwargs['membersToSkip'] are not parsed. 
        This, for example, allows the calling instance to parse the *externalFiles* member and set up needed entries into *kwargs*.

        This method returns a dictionary and a list. The dictionary contains all child nodes of *node* not parsed, including 
        those specified by membersToSkip. The keys in the dictionary are the member monikers and the items are their associated 
        child node of *node*. The list contains all monikers for the *ancestryMembers* not in *node*.

        :param node:        A :py:class:`xmlNodeMode.XML_node` instance.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary that contains data to control the way this method acts.
        """

        membersToSkip = kwargs.get('membersToSkip', [])
        membersNotFoundInNode = []
        childNodes = {}
        for child in node: childNodes[child.tag] = child

        for ancestryMember in self.ancestryMembers:
            if ancestryMember in membersToSkip:
                continue
            elif ancestryMember not in childNodes:
                membersNotFoundInNode.append(ancestryMember)
            else:
                child = childNodes.pop(ancestryMember)
                getattr(self, ancestryMember).parseNode(child, xPath, linkData, **kwargs)

        return childNodes, membersNotFoundInNode

    def parseCleanup(self, node, **kwargs):
        """
        This method does nothing. A place holder incase the derived class needs one.

        :param node:        A :py:class:`xmlNodeMode.XML_node` instance.
        :param kwargs:      A dictionary that contains data to control the way this method acts.
        """

        pass

    @classmethod
    def parseXMLString(cls, string, **kwargs):
        """
        This method parses a XML string to construct an instance of *cls*.

        :param string:      String to parse.
        :param kwargs:      A dictionary that contains data to control the way this method acts.

        :returns:           An instance of *cls*.
        """

        if not isinstance(string, str): raise TypeError('Invalid string.')

        node = ElementTree.fromstring(string)
        node = xmlNodeMode.XML_node(node, xmlNodeMode.XML_node.etree)

        instance = cls.parseNodeUsingClass(node, [], {}, **kwargs)
        instance.parseCleanup(node, **kwargs)

        return instance

    @classmethod
    def readXML_file(cls, fileName, **kwargs):
        """
        This methods reads and parses the data in file *fileName* to construct an instance of *cls*.

        :param cls:         The class to construct from the data in file *fileName*.
        :param fileName:    Path to the file used to construct an instance of *cls*.
        :param kwargs:      A dictionary that contains data to control the way this method acts.

        :returns:           An instance of *cls*.
        """

        if isinstance(fileName, pathlib.Path):
            fileName = str(fileName)
        if not isinstance(fileName, str): raise TypeError('Invalid file name.')

        node = ElementTree.parse(fileName).getroot()
        node = xmlNodeMode.XML_node(node, xmlNodeMode.XML_node.etree)

        if node.tag != cls.moniker:
            raise ValueError('Node name "%s" in XML not the same as requested class moniker "%s".' % (node.tag, cls.moniker))

        kwargs['sourcePath'] = fileName

        instance = cls.parseNodeUsingClass(node, [], {}, **kwargs)
        instance.parseCleanup(node, **kwargs)

        return instance

class AncestryIO(AncestryIO_base):
    """
    This is an abstruct class that defined the abstract method :py:class:`AncestryIO_bare.parseNodeUsingClass`.
    """

    @classmethod
    @abc.abstractmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        This method must be overrwritten by the derived class.
        The overrwritten method must construct a class of instance *cls* using *node* 

        :param node:        A :py:class:`xmlNodeMode.XML_node` instance.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary that contains data to control the way this method acts.
        """

        pass

class AncestryIO_bare(AncestryIO_base):
    """
    This is an abstruct class that defined the abstract method :py:class:`AncestryIO_bare.parseNode`.
    """

    @abc.abstractmethod
    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method must be overrwritten by the derived class.
        The overrwritten method must fill in *self* with data in *node*.

        :param node:        A :py:class:`xmlNodeMode.XML_node` instance.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary that contains data to control the way this method acts.
        """

        pass

class XPathNotFound(Exception):
    """
    Exception class which is raised when the :py:func:`followXPath` cannot follow the path.
    """

    pass
