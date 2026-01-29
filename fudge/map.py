# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__doc__ = """
This module contains classes for dealing with a **GNDS** map file. A map file creates a nuclear data library 
from a list of **GNDS** **reactionSuite** files.

Example of expected usage:

.. code-block::

    from fudge import map as mapModule
    file = '/path/to/GIDI/Test/Data/MG_MC/all_maps.map'
    map = mapModule.read( file )    # Read in a map file.
    m_O16 = map.find( 'n', 'O16' )  # Gets a reference to a map Protare instance for "n + O16".
    O16 = m_O16.read( )             # Creates a FUDGE protare (i.e., reactionSuite instance).

Note, in general, most user should get a protare (i.e., **GNDS** **reactionSuite**) this way.
"""

import os
import abc
import re
import pathlib

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from LUPY import ancestry as ancestryModule
from LUPY import checksums as checksumsModule
from LUPY import misc as LUPY_miscModule

from fudge import enums as enumsModule
from . import styles as stylesModule
from . import GNDS_file as GNDS_fileModule
from fudge.core.utilities import guessInteraction as guessInteractionModule
from fudge import reactionSuite as reactionSuiteModule


__todo = """
    -) Need to make Base and some other classes a meta class.
"""

class FormatVersion:
    """
    This class is only needed because some legacy map files (map files created before the map node was
    defined in **GNDS**) are still used at LLNL. When these map files are no longer used, this class
    will be removed. New map files should only be created with **GNDS** format '2.0'.
    """

    version_0_1 = "0.1"
    version_0_2 = "0.2"
    allowed = (
        version_0_1,
        version_0_2,
        GNDS_formatVersionModule.version_1_10,
        GNDS_formatVersionModule.version_2_0_LLNL_4,
        GNDS_formatVersionModule.version_2_0,
    )
    default = GNDS_formatVersionModule.default

    @staticmethod
    def check(format):
        """
        This methods should be used to check that the format of a map file is valid. Also, the returned value is the format
        that should be used after calling this method.

        :param format:  The format string to check.

        :return:        The format that should be used by other logic in this module.
        """

        if format == FormatVersion.version_0_1 or format == FormatVersion.version_0_2:
            format = GNDS_formatVersionModule.version_1_10
        if format not in FormatVersion.allowed:
            raise ValueError('Unsupported format "%s".' % format)

        return format


class Base(ancestryModule.AncestryIO):
    """
    Base class used by all other map related classes. Mainly inherits class :py:class:`ancestryModule.AncestryIO`
    and adds the members **checksum** and **algorithm**.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | checksum  | The check sum of the file referrence by the path member.  |
    +-----------+-----------------------------------------------------------+
    | algorithm | The algorithm used to calculate *checksum*.               |
    +-----------+-----------------------------------------------------------+
    """

    def __init__(self, checksum=None, algorithm=None):
        """
        :param checksum:        The check sum of the file referrence by the path member.
        :param algorithm:       The algorithm used to calculate *checksum*.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.checksum = checksum
        self.algorithm = algorithm

    @property
    def checksum(self):
        """Return checksum for this map entry."""

        return self.__checksum

    @checksum.setter
    def checksum(self, checksum_):
        """Sets *self*'s checksum to checksum_."""

        if checksum_ is not None and not isinstance(checksum_, str):
            raise TypeError("Checksum must be a string.")
        self.__checksum = checksum_

    @property
    def algorithm(self):
        """Returns the algorithm used to compute the checksum for this entry."""

        if self.__algorithm:
            return self.__algorithm
        if isinstance(self, Map):
            return None
        return self.ancestor.algorithm

    @algorithm.setter
    def algorithm(self, algorithm_):
        """Set *self*'s algorithm to algorithm_."""

        if algorithm_ is not None and algorithm_ not in checksumsModule.supportedAlgorithms:
            raise ValueError("Unsupported checksum algorithm '%s'." % algorithm_)
        self.__algorithm = algorithm_

    @property
    @abc.abstractmethod
    def path(self):
        pass

    def standardXML_attributes(self, checkAncestor=True):
        """
        Returns the XML attribute string for *checksum* and *algorithm*. Mainly for internal use.
        If *checkAncestor* is **True** and *self*'s *ancestor* uses the same *algorithm* as *self*, then
        the algorithm attribute string is **not** included in the returned string.

        :param checkAncestor:   Used to determine if *algorithm* needed to be added to attribute string.

        :return:                XML attribute string for *checksum* and *algorithm*.
        """

        attrs = ''
        if self.checksum:
            attrs += ' checksum="%s"' % self.checksum

        algorithm = self.algorithm
        if self.ancestor is not None:
            if checkAncestor and self.algorithm == self.ancestor.algorithm:
                algorithm = None
        if algorithm is not None:
            attrs += ' algorithm="%s"' % algorithm

        return attrs


class Map(Base):
    """
    This class represents a **GNDS** **map** that is used to create a nuclear data library from a list
    of :py:class:`Protare` :py:class:`TNSL` and :py:class:`Import` instances. Each instance in the list
    is called an entry.

    The following table list the primary members of this class:

    +---------------+-----------------------------------------------------------+
    | Member        | Description                                               |
    +===============+===========================================================+
    | library       | The name of the library.                                  |
    +---------------+-----------------------------------------------------------+
    | path          | The path to the map file read in.                         |
    +---------------+-----------------------------------------------------------+
    | parentsDir    | May be deprecated, needs more investigating.              |
    +---------------+-----------------------------------------------------------+
    | format        | The format version for the map read in.                   |
    +---------------+-----------------------------------------------------------+
    | checksum      | The check sum of the file referrence by the path member.  |
    +---------------+-----------------------------------------------------------+
    | algorithm     | The algorithm used to calculate *checksum*.               |
    +---------------+-----------------------------------------------------------+
    """

    moniker = "map"

    def __init__(self, library, path, parentsDir="", format=FormatVersion.default, checksum=None, algorithm=None):
        """
        Constructor for **Map** class.

        :param library:         The name of the library.
        :param path:            The path to the map file that was used to construct *self*.
        :param parentsDir:      May be deprecated, needs more investigating.
        :param format:          The format version for the map file read in.
        :param checksum:        Hash for the map (computed from the concatenation of all checksums inside the map)
        :param algorithm:       Algorithm used to compute checksums.
        """

        Base.__init__(self, checksum=checksum, algorithm=algorithm)

        self.__library = LUPY_miscModule.isString(library, "Library must be a string.")
        self.__path = LUPY_miscModule.isString(path, "Path must be a string.")
        parentsDir = LUPY_miscModule.isString(
            parentsDir, """Parent's path must be a string."""
        )

        fileName = pathlib.Path(path)
        if not fileName.is_absolute():
            if parentsDir != '':
                fileName = parentsDir / fileName
            fileName = fileName.resolve()
        self.__fileName = str(fileName)

        self.__format = FormatVersion.check(format)

        self.__entries = []

    def __getitem__(self, index):
        """
        Returns the (index-1)^th item of self.

        :param index:       The index of the iter to return.
        """

        return self.__entries[index]

    def __iter__(self):
        """Iterators over each entry in *self*. Does not dive into import entries."""

        for entry in self.__entries:
            yield entry

    def __len__(self):
        """Returns the number of entries in *self*. Does not dive into import entries."""

        return len(self.__entries)

    @property
    def path(self):
        """
        Returns the path this map was read from or may be written to. It may be absolote or relative,
        depending on how it was initialize.
        """

        return self.__path

    @property
    def format(self):
        """Returns to format for *self*."""

        return self.__format

    @property
    def library(self):
        """Returns the library name for *self*."""

        return self.__library

    @property
    def fileName(self):
        """
        Returns the file name for the location of *self*. Unlike the method :py:func:`path`,
        this will always be the absolute path.
        """

        return self.__fileName

    def updateAllChecksums(
        self, algorithm=checksumsModule.Sha1sum.algorithm, mapDirectory=None
    ):
        """Calls *updateChecksum* on all entries and then updates *self*'s *checksum*."""

        import threading
        from LUPY import parallelprocessing

        class computeSum(threading.Thread):
            def __init__(self, entry):
                threading.Thread.__init__(self)
                self.entry = entry

            def run(self):
                self.entry.updateChecksum(algorithm, mapDirectory=mapDirectory)

        entries = list(self)
        parallelprocessing.executeMultiThreaded(entries, runner=computeSum)

        self.updateChecksum(algorithm)

    def updateChecksum(self, algorithm=checksumsModule.Sha1sum.algorithm):
        """
        Computes map file checksum by concatenating checksums for all entries of *self* into a string and compute the
        checksum of that string.

        :param algorithm:           The algorithm to use to calculate the checksum.
        """

        checker = checksumsModule.checkers[
            algorithm
        ]  # Caleb, I think checker and checkers should be called algorithm and algorithms.
        s1 = "".join([entry.checksum for entry in self if entry.checksum is not None])
        self.checksum = checker.from_string(s1)
        self.algorithm = algorithm

    def append(self, entry):
        """Appends *entry* to *self*."""

        if not isinstance(entry, EntryBase):
            TypeError("Invalid entry.")

        self.__entries.append(entry)
        entry.setAncestor(self)

    def insert(self, index, entry):
        """Inserts the entry into *self* at index."""

        if not isinstance(entry, EntryBase):
            TypeError('Invalid entry of type "%s".' % type(entry))

        self.__entries.insert(index, entry)
        entry.setAncestor(self)

    def iterate(self):
        """
        Iterates over all protares and TNSLs in *self* including those in import entries.
        Ergo, unlike :py:class:`__getitem__`, this method dives into import entries.
        That is, no :py:class:`Import` instances is returned.
        """

        for entry in self.__entries:
            if isinstance(entry, Import):
                yield from entry.iterate()
            else:
                yield entry

    def find(self, projectile, target, library=None, evaluation=None):
        """
        Returns the first entry matching *projectile*, *target*, *library* and *evaluation*.
        Searches each entry in the order they were appended and dives into imported maps.
        If *library* is **None**, then all libraries are treated as matched.
        If *evaluation* is **None**, then all evaluations are treated as matched.

        :param projectile:      The **GNDS** **PoPs** id for the projectile to match.
        :param target:          The **GNDS** **PoPs** id for the target to match.
        :param library:         The name of the library to match.
        :param evaluation:      The name of the evaluation to match.

        :return:            Returns the found entry or **None** is no match was found.
        """

        for entry in self.__entries:
            if isinstance(entry, Import):
                foundEntry = entry.find(projectile, target, library, evaluation)
                if foundEntry is not None:
                    return foundEntry
            else:
                if entry.isMatch(projectile, target, library, evaluation):
                    return entry

        return None

    def findAllOf(self, projectile=None, target=None, library=None, evaluation=None):
        """
        Returns a list of all entries matching *projectile*, *target*, *library* and *evaluation*. Searches
        all imported maps. If *library* is **None**, the all libraries are treated as matched.
        If *evaluation* is **None**, the all evaluations are treated as matched.

        :param projectile:      The **GNDS** **PoPs** id for the projectile to match.
        :param target:          The **GNDS** **PoPs** id for the target to match.
        :param library:         The name of the library to match.
        :param evaluation:      The name of the evaluation to match.

        :return:                Returns the list of all entries that match the input arguments.
        """

        allFound = []

        for entry in self.__entries:
            if isinstance(entry, Import):
                allFound += entry.findAllOf(projectile, target, library, evaluation)
            else:
                if entry.isMatch(projectile, target, library, evaluation):
                    allFound.append(entry)

        return allFound

    def unrecurseAndSave(self, newFileName):
        """
        Iterate through all entries recursively and save as a flattened map files to file *newFileName*.
        Ergo, produces a similar :py:class:`Map` instance but with no import entries.

        :param newFileName:  New filename to which flattened mapfile will be saved.

        :return:            Returns the flattened :py:class:`Map` instance.
        """

        newInstance = Map(self.library, newFileName, "", self.format, algorithm=self.algorithm)
        pathList = []
        for entry in self.iterate():
            entry.path = pathlib.Path(entry.path).resolve().as_posix()
            if entry.path not in pathList:
                newInstance.append(entry)
                pathList.append(entry.path)
            else:
                raise RuntimeError(f'Duplicate ({entry.projectile},{entry.target}) entry in map file {self.path}.') 

        newInstance.updateChecksum()
        newInstance.saveToFile(newFileName)

        return newInstance

    def toXML_strList(self, indent="", **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A keyword list.

        :return:                List of str instances representing the XML lines of *self*.
        """

        indent2 = indent + kwargs.get("incrementalIndent", "  ")

        format = FormatVersion.check(kwargs.get('format', kwargs.get('formatVersion', FormatVersion.default)))

        attrs = self.standardXML_attributes(False)

        XML_list = [
            '%s<%s library="%s" format="%s"%s>'
            % (indent, self.moniker, self.library, format, attrs)
        ]
        for entry in self.__entries:
            XML_list += entry.toXML_strList(indent2, **kwargs)
        XML_list[-1] += "</%s>" % self.moniker

        return XML_list

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Creates a Map instance from an map node.

        :param node:        Node to parse.
        :param xPath:       Currently not used.
        :param linkData:    Currently not used.
        :param kwargs:      A keyword list.

        :return:            :py:class:`Map` instance.
        """

        if node.tag != Map.moniker:
            raise TypeError('Invalid node name "%s" for a map file.' % node.tag)

        sourcePath = kwargs["sourcePath"]

        format = node.get("format", None)
        if format is None:
            raise ValueError("Map node does not have 'format' attribute.")

        library = node.get("library", None)
        if library is None:
            raise ValueError("Map node does not have 'library' attribute.")

        map = Map(library, sourcePath, checksum=node.get('checksum'), algorithm=node.get('algorithm'), format=format)

        kwargs["format"] = FormatVersion.check(format)
        for child in node:
            if child.tag == Import.moniker:
                map.append(Import.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            elif child.tag == Protare.moniker:
                map.append(Protare.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            elif child.tag == TNSL.moniker:
                map.append(TNSL.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            else:
                raise ValueError("Invalid child tag '%s' for map file." % child.tag)

        return map

    @staticmethod
    def read(fileName, **kwargs):
        """
        Reads in the file name *fileName* and returns a **Map** instance.

        :param fileName:    The file to read in and parse.
        :param kwargs:      A keyword list.

        :return:            :py:class:`Map` instance.
        """

        return Map.readXML_file(fileName, **kwargs)


class EntryBase(Base):
    """
    Base class for all map entry classes.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | path      | The file path to the file to read in for this entry.      |
    +-----------+-----------------------------------------------------------+
    | checksum  | The check sum of the file referrence by the path member.  |
    +-----------+-----------------------------------------------------------+
    | algorithm | The algorithm used to calculate *checksum*.               |
    +-----------+-----------------------------------------------------------+
    """

    def __init__(self, path, checksum=None, algorithm=None):
        """
        Base constructor for map entries.

        :param path:            Path attribute for the entry.
        :param checksum:        Checksum for this entry, computed using the indicated algorithm.
        :param algorithm:       Algorithm for computing the checksum. Only required if different from parent <map> file.
        """

        Base.__init__(self, checksum=checksum, algorithm=algorithm)
        self.path = path

    @property
    def fileName(self):
        """
        Returns the file name for the location of *self*. Unlike the method :py:func:`path`,
        this will always be the absolute path.
        """

        fileName = pathlib.Path(self.path)
        if not fileName.is_absolute():
            fileName = pathlib.Path(self.ancestor.fileName).parent / fileName
            fileName = fileName.resolve()

        return str(fileName)

    @property
    def path(self):
        """Returns the value of the *__path* member. This may be absolute or relative."""

        return self.__path

    @path.setter
    def path(self, path):
        """
        Sets member *self*.__path to *path*.

        :param path:    The name of the path for *self*.
        """

        if isinstance(path, pathlib.Path):
            path = str(path)

        if not isinstance(path, str):
            raise TypeError("Path must be a string.")
        self.__path = path

    def buildFileName(self, mapDirectory=None):
        """
        This method is designed to aid in building a map file when the map file does not reside in its final resting place.
        The value of *mapDirectory* should be the final resting place (i.e., directory) of the map file.
        Otherwise, this method is the same as the propery methods *fileName*.

        :param mapDirectory:    The final directory where the map file will reside if it is different than "./".
        """

        if mapDirectory is None:
            return self.fileName

        return str(pathlib.Path(pathlib.Path(mapDirectory).resolve()) / self.path)

    def updateChecksum(self, algorithm=checksumsModule.Sha1sum.algorithm, mapDirectory=None):
        """
        Compute the checksum of the file specified by *path* member and store it into the *checksum* member.

        @param algorithm:       The algorithm to use for the checksum.
        """

        self.algorithm = algorithm
        self.checksum = checksumsModule.checkers[algorithm].from_file(self.buildFileName(mapDirectory))

class Import(EntryBase):
    """
    Class representing an 'import' entry in a map. An import entry specifies the location of another map file to import.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | path      | The file path to the map file in for this entry.          |
    +-----------+-----------------------------------------------------------+
    | checksum  | The check sum of the file referrence by the path member.  |
    +-----------+-----------------------------------------------------------+
    | algorithm | The algorithm used to calculate *checksum*.               |
    +-----------+-----------------------------------------------------------+
    """

    moniker = "import"

    def __init__(self, path, checksum=None, algorithm=None):
        """Constructor for the import entry.

        :param path:            Path to the map file in import.
        :param checksum:        Checksum for this entry, computed using the indicated algorithm.
        :param algorithm:       Algorithm for computing the checksum. Only required if different from parent <map> file.
        """

        EntryBase.__init__(self, path, checksum=checksum, algorithm=algorithm)
        self.__map = None

    def __str__(self):
        """Returns a simple string representation of *self*."""

        return '%s with path "%s".' % (self.moniker, self.path)

    @property
    def derivedPath(self):
        """
        Returns the parent part of the ancestor's path. This, together with *self*'s path defines the
        location of the map file to read in.

        :return:        The parent part of the ancestor's path.
        """

        ancestor = self.ancestor
        if isinstance(ancestor, Map):
            return pathlib.Path(ancestor.fileName).parent

        raise TypeError('Import not an entry of a Map instance.')

    @property
    def map(self):
        """Returns a reference to *self*.__map."""

        self.readMap()

        return self.__map

    def find(self, projectile, target, library=None, evaluation=None):
        """
        Calls :py:func:`Map.find` on *self*.__map and returns its results.

        :param projectile:      The name of the projectile to match.
        :param target:          The name of the projectile to match.
        :param library:         The name of the library to match.
        :param evaluation:      The name of the evaluation to match.

        :return:                Returns the found entry or **None** is no match was found.
        """

        self.readMap()

        return self.__map.find(projectile, target, library, evaluation)

    def findAllOf(self, projectile=None, target=None, library=None, evaluation=None):
        """
        Calls :py:func:`Map.findAllOf` on *self*.__map and returns its results.
        Arguments that are **None** match all values of their type. For example, if *projectile* is **None**,
        then all projectiles are a match.

        :param projectile:      The name of the projectile to match.
        :param target:          The name of the projectile to match.
        :param library:         The name of the library to match.
        :param evaluation:      The name of the evaluation to match.

        :return:                Returns the list of all matches found.
        """

        self.readMap()
        return self.__map.findAllOf(projectile, target, library, evaluation)

    def iterate(self):
        """
        Iterates over all protares and TNSLs in *self*'s map including nested :py:class:`Map` instances.
        No :py:class:`Import` instances is returned.
        """

        self.readMap()
        yield from self.map.iterate()

    def readMap(self):
        """Reads in the map file pointed to by *self* if not already read in. An import only reads its map file when needed."""

        if self.__map is None:
            mapFilePath = pathlib.Path(self.path)
            if not mapFilePath.is_absolute():
                mapFilePath = self.derivedPath / mapFilePath
            self.__map = Map.readXML_file(mapFilePath)
            self.__map.setAncestor(self)

        return self.__map

    def toXML_strList(self, indent="", **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A keyword list.

        :return:                List of str instances representing the XML lines of *self*.
        """

        attrs = self.standardXML_attributes()
        return ['%s<%s path="%s"%s/>' % (indent, self.moniker, self.path, attrs)]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Creates an :py:class:`Import` instance from an import node.

        :param node:            Node to parse.
        :param xPath:           Currently not used.
        :param linkData:        Currently not used.
        :param kwargs:          A keyword list.

        :return:                :py:class:`Import` instance.
        """

        if node.tag != Import.moniker:
            raise TypeError("Invalid node name.")

        kwargs = {attr: node.get(attr) for attr in ("path", "checksum", "algorithm")}
        _import = Import(**kwargs)
        return _import


class ProtareBase(EntryBase):
    """
    Base class for all protare-like entry classes (i.e., the :py:class:`Protare` and :py:class:`TNSL` classes).

    The following table list the primary members of this class:

    +---------------+-----------------------------------------------------------+
    | Member        | Description                                               |
    +===============+===========================================================+
    | projectile    | The name of the projectile.                               |
    +---------------+-----------------------------------------------------------+
    | target        | The name of the target.                                   |
    +---------------+-----------------------------------------------------------+
    | evaluation    | The name of the evaluation.                               |
    +---------------+-----------------------------------------------------------+
    | interaction   | The type of interacter for the protare.                   |
    +---------------+-----------------------------------------------------------+
    | path          | The file path to the file to read in for this entry.      |
    +---------------+-----------------------------------------------------------+
    | checksum      | The check sum of the file referrence by the path member.  |
    +---------------+-----------------------------------------------------------+
    | algorithm     | The algorithm used to calculate *checksum*.               |
    +---------------+-----------------------------------------------------------+
    """

    def __init__(self, projectile, target, evaluation, path, interaction, checksum=None, algorithm=None):
        """
        Construtor for ProtareBase instance.

        :param projectile:          Name for the projectile.
        :param target:              Name for the target.
        :param evaluation:          Name for the protare's Evaluation.
        :param path:                Path to the protare file.
        :param interaction:         Type of interaction for the data in the protare file (see class reactionSuite.Interaction).
        :param checksum:            Checksum for this entry, computed using the indicated algorithm.
        :param algorithm:           Algorithm for computing the checksum. Only required if different from parent <map> file.
        """

        EntryBase.__init__(self, path, checksum=checksum, algorithm=algorithm)

        if not isinstance(projectile, str):
            raise TypeError("Projectile must be a string.")
        self.__projectile = projectile

        if not isinstance(target, str):
            raise TypeError("Target must be a string.")
        self.__target = target

        if not isinstance(evaluation, str):
            raise TypeError("Evaluation must be a string.")
        self.__evaluation = evaluation

        if isinstance(self, TNSL):
            interaction = None  # In GNDS 2.0, TNSL no longer as interaction attribution as it is not needed.
        else:
            interaction = enumsModule.Interaction.checkEnumOrString(interaction)
        self.__interaction = interaction

        self.__guessedInteraction = False

    @property
    def projectile(self):
        """Returns the projectile's name."""

        return self.__projectile

    @property
    def target(self):
        """Returns the target's name."""

        return self.__target

    @property
    def evaluation(self):
        """Returns the evaluation string."""

        return self.__evaluation

    @property
    def interaction(self):
        """Returns the interaction attribute."""

        return self.__interaction

    @property
    def guessedInteraction(self):
        """
        Returns the guessedInteraction value. This is for pre GNDS 2.0 use only and should not be used with GNDS 2.0 or higher
        (i.e., it should always return **False** for GNDS 2.0 or higher).
        """

        return self.__guessedInteraction

    @guessedInteraction.setter
    def guessedInteraction(self, value):
        """
        Set *self*'s guessedInteraction member. This is for pre GNDS 2.0 use only and should not be
        used with GNDS 2.0 or higher.
        """

        self.__guessedInteraction = value

    @property
    def library(self):
        """Returns the name of the library the protare resides in."""

        return self.ancestor.library

    @property
    def derivedPath(self):
        """
        Returns the location of the file pointed to by the *path* member relative to the ancestor's locations.
        """

        ancestor = self.ancestor
        if self.path[0] == os.sep or not isinstance(ancestor, Base):
            return self.path
        return os.path.join(os.path.dirname(ancestor.path), self.path)

    def isMatch(self, projectile, target, library=None, evaluation=None):
        """
        Returns **True** is *self* matches projectile, target, library and evaluation, and **False** otherwise.
        If an argument has a value of **None**, that argument is a match.
        For example, isMatch( None, "O16", None, None ) will match any entry with a target of "O16".

        :param projectile:      The requested projectile's **PoPs** id to check for match.
        :param target:          The requested target's **PoPs** id to check for match.
        :param library:         The requested library to check for match.
        :param evaluation:      The requested evaluation to check for match.

        :return:                Returns **True** if a match and **False** otherwise.
        """

        if self.adaptable(self.projectile, projectile):
            if self.adaptable(self.target, target):
                if self.adaptable(self.library, library):
                    return self.adaptable(self.__evaluation, evaluation)

        return False

    def preview(self, haltParsingMoniker=stylesModule.Styles.moniker):
        """
        This method calls :py:func:`GNDS_fileModule.preview with *haltParsingMoniker*` with *self*'s *fileName* and returns
        its retults.  See :py:func:`GNDS_fileModule.preview for more details.

        :param haltParsingMoniker:  The child node of the protare to halt parsing on.

        :return:                    An instance of :py:class:`reactionSuiteModule.ReactionSuite`.
        """

        return GNDS_fileModule.preview(self.fileName, haltParsingMoniker=haltParsingMoniker)

    def protare(self, verbosity=0, lazyParsing=True):
        """Calls *self*.read. This method is deprecated."""

        return self.read(verbosity=verbosity, lazyParsing=lazyParsing)

    def read(self, verbosity=0, lazyParsing=True):
        """
        Reads in the protare using the **read** function in the module ReactionSuite.py and returns the instance.

        :param verbosity:       Passed onto ReactionSuite.read.
        :param lazyParsing:     Passed onto ReactionSuite.read.
        """

        return reactionSuiteModule.read(self.fileName, verbosity=verbosity, lazyParsing=lazyParsing)

    def standardXML_attributes(self, checkAncestor=True):
        """
        Returns the XML attribute string for protare, target, evalution and interaction as well that those
        returned by :py:func:`Base.standardXML_attributes`.

        :param checkAncestor:   Passed to the method :py:func:`Base.standardXML_attributes`.

        :returns:               A XML attribute string.
        """

        attributeString = ' projectile="%s" target="%s" evaluation="%s" path="%s"' % (self.projectile, self.target, self.evaluation, self.path)
        if not self.guessedInteraction:
            if self.interaction is not None:
                attributeString += ' interaction="%s"' % self.interaction

        return attributeString + Base.standardXML_attributes(self, checkAncestor=checkAncestor)

    @staticmethod
    def adaptable(item1, item2):
        """
        Returns **True** whenever item2 is **None** or *item1* and *item2* are the same. For internal use.

        :param item1:       Item to compare to item2.
        :param item2:       Item to compare to item1.
        """

        if item2 is None:
            return True

        if item1 == item2:      # For case when TNSL target has a regular expression special character in its name.
            return True

        return re.fullmatch(item2, item1) is not None


class Protare(ProtareBase):
    """
    The class for a map's protare child node.

    The following table list the primary members of this class:

    +---------------+-----------------------------------------------------------+
    | Member        | Description                                               |
    +===============+===========================================================+
    | projectile    | The name of the projectile.                               |
    +---------------+-----------------------------------------------------------+
    | target        | The name of the target.                                   |
    +---------------+-----------------------------------------------------------+
    | evaluation    | The name of the evaluation.                               |
    +---------------+-----------------------------------------------------------+
    | interaction   | The type of interacter for the protare.                   |
    +---------------+-----------------------------------------------------------+
    | path          | The file path to the file to read in for this entry.      |
    +---------------+-----------------------------------------------------------+
    | checksum      | The check sum of the file referrence by the path member.  |
    +---------------+-----------------------------------------------------------+
    | algorithm     | The algorithm used to calculate *checksum*.               |
    +---------------+-----------------------------------------------------------+
    """

    moniker = "protare"

    def __init__(self, projectile, target, evaluation, path, interaction, checksum=None, algorithm=None):
        """
        Construtor for a Protare instance.

        :param projectile:          Name for the projectile.
        :param target:              Name for the target.
        :param evaluation:          Name for the protare's Evaluation.
        :param path:                Path to the protare file.
        :param interaction:         Type of interaction for the data in the protare file (see class reactionSuite.Interaction).
        :param checksum:            Checksum for this entry, computed using the indicated algorithm.
        :param algorithm:           Algorithm for computing the checksum. Only required if different from parent <map> file.
        """

        ProtareBase.__init__(self, projectile, target, evaluation, path, interaction,
                             checksum=checksum, algorithm=algorithm)

    def __str__(self):
        """Returns a simple string representation of *self*."""

        return ('%s with projectile "%s", target "%s", evaluation "%s", path "%s" and interaction "%s".' %
                (self.moniker, self.projectile, self.target, self.evaluation, self.path, self.interaction))

    def toXML_strList(self, indent="", **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A keyword list.

        :return:                List of str instances representing the XML lines of *self*.
        """

        format = FormatVersion.check(kwargs.get('format', kwargs.get('formatVersion', FormatVersion.default)))

        return ["%s<%s%s/>" % (indent, self.moniker, self.standardXML_attributes())]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Creates a Protare instance from a protare node.

        :param node:        Node to parse.
        :param xPath:       Currently not used.
        :param linkData:    Currently not used.
        :param kwargs:      A keyword list.

        :return:            A :py:class:`Protare` instance.
        """

        if node.tag != Protare.moniker:
            raise TypeError("Invalid node name.")

        kwargs = {attr: node.get(attr) for attr in
                  ('projectile', 'target', 'evaluation', 'path', 'interaction', 'checksum', 'algorithm')}

        guessedInteraction, kwargs['interaction'] = guessInteractionModule.guessInteraction(
            kwargs['interaction'], kwargs['projectile'], kwargs['target'])

        protare = Protare(**kwargs)

        protare.guessedInteraction = guessedInteraction

        return protare


class TNSL(ProtareBase):
    """
    The class for a map's TNSL node.

    The following table list the primary members of this class:

    +-----------------------+-----------------------------------------------------------+
    | Member                | Description                                               |
    +=======================+===========================================================+
    | projectile            | The name of the projectile.                               |
    +-----------------------+-----------------------------------------------------------+
    | target                | The name of the target.                                   |
    +-----------------------+-----------------------------------------------------------+
    | evaluation            | The name of the evaluation.                               |
    +-----------------------+-----------------------------------------------------------+
    | path                  | The file path to the file to read in for this entry.      |
    +-----------------------+-----------------------------------------------------------+
    | standardTarget        | Name for the standard target.                             |
    +-----------------------+-----------------------------------------------------------+
    | standardEvaluation    | Name for the standard evaluation.                         |
    +-----------------------+-----------------------------------------------------------+
    | checksum              | The check sum of the file referrence by the path member.  |
    +-----------------------+-----------------------------------------------------------+
    | algorithm             | The algorithm used to calculate *checksum*.               |
    +-----------------------+-----------------------------------------------------------+
    """

    moniker = "TNSL"

    def __init__(self, projectile, target, evaluation, path, standardTarget, standardEvaluation, 
                interaction=None, checksum=None, algorithm=None):
        """
        Construtor for a TNSL instance.

        :param projectile:              Name for the projectile.
        :param target:                  Name id for the target.
        :param evaluation:              Evaluation for the protare.
        :param path:                    Path to the protare file.
        :param standardTarget:          Name for the standard target.
        :param standardEvaluation:      Name for the standard evaluation.
        :param checksum:                Checksum for this entry, computed using the indicated algorithm.
        :param algorithm:               Algorithm for computing the checksum. Only required if different from parent <map> file.
        """

        ProtareBase.__init__(self, projectile, target, evaluation, path, None,
                             checksum=checksum, algorithm=algorithm)

        if not isinstance(standardTarget, str):
            raise TypeError("Standard target must be a string.")
        self.__standardTarget = standardTarget

        if not isinstance(standardEvaluation, str):
            raise TypeError(
                'Standard evaluation must be a string, is of type "%s".'
                % type(standardEvaluation)
            )
        self.__standardEvaluation = standardEvaluation

    def __str__(self):
        """Returns a simple string representation of *self*."""

        return (
            '%s with projectile "%s", target "%s", evaluation "%s", path "%s", standardTarget "%s" and standardEvaluation "%s".'
            % (self.moniker, self.projectile, self.target, self.evaluation, self.path, self.standardTarget, self.standardEvaluation)
        )

    @property
    def standardTarget(self):
        """Returns the standardTarget."""

        return self.__standardTarget

    @property
    def standardEvaluation(self):
        """Returns the standardEvaluation."""

        return self.__standardEvaluation

    def toXML_strList(self, indent="", **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A keyword list.

        :return:                List of str instances representing the XML lines of *self*.
        """

        format = FormatVersion.check(kwargs.get('format', kwargs.get('formatVersion', FormatVersion.default)))
        attrs = self.standardXML_attributes()

        if format == FormatVersion.version_0_1 or format == GNDS_formatVersionModule.version_1_10:
            indent2 = indent + kwargs.get('incrementalIndent', '  ')
            XML_list = ['%s<%s%s>' % (indent, self.moniker, attrs)]
            XML_list.append('%s<protare projectile="n" target="%s" evaluation="%s"/></%s>' %
                            (indent2, self.standardTarget, self.standardEvaluation, self.moniker))
        else:
            XML_list = [
                '%s<%s%s standardTarget="%s" standardEvaluation="%s"/>'
                % (indent, self.moniker, attrs, self.standardTarget, self.standardEvaluation)
            ]

        return XML_list

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Creates a TNSL instance from a TNSL node.

        :param node:            Node to parse.
        :param xPath:           Current not used.
        :param linkData:        Currently not used.
        :param kwargs:          A keyword list.

        :return:                A :py:class:`TNSL` instance.
        """

        if node.tag != TNSL.moniker:
            raise TypeError("Invalid node name.")

        format = FormatVersion.check(
            kwargs.get("format", kwargs.get("formatVersion", FormatVersion.default))
        )

        kwargs = {attr: node.get(attr) for attr in
                  ('projectile', 'target', 'evaluation', 'path', 'checksum', 'algorithm')}

        if format == FormatVersion.version_0_1 or format == GNDS_formatVersionModule.version_1_10:
            kwargs['standardTarget'] = node[0].get('target', None)
            kwargs['standardEvaluation'] = node[0].get('evaluation', None)
        else:
            kwargs['standardTarget'] = node.get('standardTarget', None)
            kwargs['standardEvaluation'] = node.get('standardEvaluation', None)

        return TNSL(**kwargs)


def read(fileName, **kwargs):
    """
    Reads in the file named *fileName* and returns a :py:class:`Map` instance.
    """

    return Map.read(fileName, **kwargs)
