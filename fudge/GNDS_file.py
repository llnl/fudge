# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the function **type** for determining the type (i.e., *ReactionSuite*, *CovarianceSuite* or *PoPs*) 
of a **GNDS/XML** file and the function **read** for reading into **FUDGE** a *GNDS/XML* file.

If the file is not a *GNDS/XML* file, a raise is executed.
"""

import pathlib
import xml.sax
from xml.etree import ElementTree

from LUPY.hdf5 import HDF5_present, h5py
from LUPY import xmlNode as xmlNodeModule  # wrapper around the xml parser:

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from PoPs import database as databaseModule

from fudge import enums as enumsModule
from fudge import map as mapModule
from fudge import suites as suitesModule
from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule
from fudge.covariances import covarianceSuite as covarianceSuiteModule
from fudge.processing import flux as fluxModule
from fudge.processing import group as groupModule
from fudge.outputChannelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule

HDF5_values = 'HDF5_values'

class GNDSTypeException(Exception):
    """
    For internal use. The raise executed in **GNDSTypeHandler.startElement** if the file is a valid
    *GNDS/XML* file as there is not other way to stop the parser.
    """

    pass

class GNDSTypeHandler(xml.sax.handler.ContentHandler):
    """For internal use. The SAX handler used to determine the *GNDS/XML* file type."""

    def startElement(self, name, attributes):
        """
        The SAX handler's startElement. This method always raises on the first start element to stop the parser.
        If the exception GNDSTypeException is raise, a valid tag name for the first element was found. All
        other exceptions are an error.
        """

        self.name = name
        if name in [reactionSuiteModule.ReactionSuite.moniker, covarianceSuiteModule.CovarianceSuite.moniker]:
            interaction = attributes.get('interaction', None)               # Allow None to support GNDS 1.10. This will be deprecated.
            if interaction == enumsModule.Interaction.legacyTNSL:
                interaction = enumsModule.Interaction.TNSL

            self.data = {'projectile': attributes['projectile'],
                         'target': attributes['target'],
                         'evaluation': attributes['evaluation'],
                         'format': attributes['format'],
                         'interaction': interaction}
        elif name == databaseModule.Database.moniker:
            self.data = {'name': attributes['name'], 'format': attributes['format']}
        elif name == mapModule.Map.moniker:
            self.data = {'library': attributes['library'], 'format': attributes['format']}
        elif name == fissionFragmentDataModule.FissionFragmentData.moniker:
            self.data = {'fissionFragmentData': None}
        elif name == fluxModule.Fluxes.moniker:
            self.data = {'fluxes': None}
        elif name == groupModule.Groups.moniker:
            self.data = {'groups': None}
        else:
            raise TypeError('Invalid XML file with top element = "%s"' % name)
        raise GNDSTypeException()


def type(fileName, show=False, checkForHDF5=True):
    """
    This function determines the type of the *GNDS/XML* file given by *fileName* or raises an error if
    the file is not a *GNDS/XML* file. For a *GNDS/XML* file this function returns a tuple of length 2
    items. The first item is the moniker for the *GNDS/XML* type. The second item is a tuple whose contents
    depend of the *GNDS/XML* file type. If the type is *ReactionSuite* or *covarianceSuite*, then the tuple
    contains 3 items. The first is the projectile's ID, the second is the target's ID and the third is the 
    evaluation string.  If the type is *PoPs*, the tuple contains 1 item which is the database's name.

    If the file is not a *GNDS/XML* file, a raise is executed.
    """

    if pathlib.Path(fileName).is_dir():
        if show:
            print('%-20s ERROR: %s' % ("DIRECTORY", fileName))
            return
        else:
            raise

    if HDF5_present and checkForHDF5:
        try:
            with h5py.File(fileName, 'r') as hdf5:
                if 'iData' not in hdf5.keys():
                    raise Exception('Does not matter as caught by except.')
            return HDF5_values, {'format': '1.10'}
        except:
            pass

    parser = xml.sax.make_parser()
    handler = GNDSTypeHandler()
    parser.setContentHandler(handler)

    try:
        parser.parse(str(fileName))
    except GNDSTypeException:
        if show: print('%-16s %s: %s' % (handler.name, fileName, handler.data))
        return handler.name, handler.data
    except xml.sax._exceptions.SAXParseException:
        if show:
            print('%-20s ERROR: %s' % ("INVALID XML", fileName))
        else:
            raise
    except:
        if show:
            print('%-20s ERROR: %s' % ("Oops", fileName))
        else:
            raise

def read(fileName, reactionSuite=None, warningNoReactionSuite=True, verbosity=1, lazyParsing=True):
    """
    This function uses the function **type** to determine the proper **FUDGE** function to call to read 
    in the file *fileName* into **FUDGE**.  It returns the **FUDGE** instance for the type.
    """

    name, dummy = type(fileName)
    if name == reactionSuiteModule.ReactionSuite.moniker:
        kwargs = {'verbosity': verbosity, 'lazyParsing': lazyParsing}
        return reactionSuiteModule.ReactionSuite.readXML_file(fileName, **kwargs)
    elif name == covarianceSuiteModule.CovarianceSuite.moniker:
        kwargs = {'reactionSuite': reactionSuite, 'warningNoReactionSuite': warningNoReactionSuite}
        return covarianceSuiteModule.CovarianceSuite.readXML_file(fileName, **kwargs)
    elif name == mapModule.Map.moniker:
        return mapModule.Map.readXML_file(fileName)
    elif name == databaseModule.Database.moniker:
        return databaseModule.read(fileName)
    elif name == fluxModule.Fluxes.moniker:
        return fluxModule.read(fileName)
    elif name == groupModule.Groups.moniker:
        return groupModule.read(fileName)
    elif name == fissionFragmentDataModule.FissionFragmentData.moniker:
        element = ElementTree.parse(fileName).getroot()
        element = xmlNodeModule.XML_node(element, xmlNodeModule.XML_node.etree)
        fissionFragmentData = fissionFragmentDataModule.FissionFragmentData()
        fissionFragmentData.parseNode(element, [], {})
        return fissionFragmentData
    else:
        if HDF5_present: return h5py.File(fileName, 'r')
        raise ImportError('Cannot open HDF5 file as h5py module not available.')


class GNDSTypeHandlerPreview(xml.sax.handler.ContentHandler):
    """For internal use. The SAX handler used to determine the *GNDS/XML* file type."""

    def __init__(self, rootMoniker, parser, childMonikers):

        xml.sax.handler.ContentHandler.__init__(self)

        self.parser = parser
        self.childMonikers = childMonikers
        self.GNDS_level = 0
        self.rootMoniker = rootMoniker
        self.haltParsingMoniker = rootMoniker
        self.GNDS_previewLine = -1
        self.GNDS_previewColumn = None

    def startElement(self, name, attributes):

        if self.GNDS_level == 1:
            if name not in self.childMonikers:
                if self.GNDS_previewLine == -1:
                    self.GNDS_previewColumn = self.parser.getColumnNumber()
                    self.GNDS_previewLine = self.parser.getLineNumber()
                raise GNDSTypeException()
        self.GNDS_level += 1

    def endElement(self, name):

        self.haltParsingMoniker = name
        self.GNDS_previewColumn = self.parser.getColumnNumber()
        self.GNDS_previewLine = self.parser.getLineNumber()
        self.GNDS_level -= 1

        if name == self.rootMoniker:
            raise GNDSTypeException()

def preview(fileName, haltParsingMoniker=stylesModule.Styles.moniker):
    """
    Returns a GNDS reactionSuite or covarianceSuite instance that only contains the child nodes up to and including the child
    node whose moniker matchs *haltParsingMoniker*. The default *haltParsingMoniker* is the "styles" node which, by GNDS 2.0 
    ordering, means that the return GNDS reactionSuite will contanin the "externalFiles" and "styles" child nodes of the input 
    file.  If *haltParsingMoniker* is None, then an empty object is returned. If the input file is not a GNDS reactionSuite 
    or covarianceSuite instance a raise is executed.
    """

    name, data = type(fileName)

    if name == reactionSuiteModule.ReactionSuite.moniker:
        childNodeOrder = reactionSuiteModule.ReactionSuite.childNodeOrder[data['format']]
    elif name == covarianceSuiteModule.CovarianceSuite.moniker:
        childNodeOrder = covarianceSuiteModule.CovarianceSuite.childNodeOrder[data['format']]
    else:
        raise ValueError('File "%s" is not a support preview GNDS type.' % fileName)

    if haltParsingMoniker is None:
        childMonikers = []
    else:
        if haltParsingMoniker not in childNodeOrder:
            raise ValueError('Invalid haltParsingMoniker = "%s".' % haltParsingMoniker)
        childMonikers = childNodeOrder[:childNodeOrder.index(haltParsingMoniker) + 1]

    parser = xml.sax.make_parser()
    handler = GNDSTypeHandlerPreview(name, parser, childMonikers)
    parser.setContentHandler(handler)

    try:
        parser.parse(fileName)

    except GNDSTypeException:
        lines = []
        with open(fileName) as fIn:
            for index in range(handler.GNDS_previewLine): lines.append(fIn.readline())

        lines.append(lines.pop(-1)[:handler.GNDS_previewColumn] + '</%s>' % handler.haltParsingMoniker)
        if handler.haltParsingMoniker != name: lines.append('</%s>' % name)

        element = ElementTree.fromstring(''.join(lines))
        element = xmlNodeModule.XML_node(element, xmlNodeModule.XML_node.etree)
        linkData = {'unresolvedLinks': []}
        if name == reactionSuiteModule.ReactionSuite.moniker:
            reactionSuite = reactionSuiteModule.ReactionSuite.parseNodeUsingClass(element, [], linkData,
                                                                                  sourcePath=fileName,
                                                                                  numberOfBrokenLinksToPrint=0)
            if ((reactionSuite.format == GNDS_formatVersionModule.version_1_10) and (
                    haltParsingMoniker == suitesModule.ExternalFiles.moniker)):
                if len(reactionSuite.externalFiles) == 0:
                    # Some older files had the externalFiles node after the styles node.
                    reactionSuite2 = preview(fileName)
                    if len(reactionSuite2.externalFiles) > 0:
                        reactionSuite = reactionSuite2
                        while len(reactionSuite.styles) > 0: reactionSuite.styles.remove(reactionSuite.styles[0])
            return reactionSuite
        else:
            return covarianceSuiteModule.CovarianceSuite.parseNodeUsingClass(element, [], linkData, sourcePath=fileName)
