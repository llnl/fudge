# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes for reading in an RIS file.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | RIS                               | This class stores the contents of an RIS file.                        |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Protare                           | This class represents one protare in an RIS file.                     |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Reaction                          | This class represents one reaction in an RIS file.                    |
    +-----------------------------------+-----------------------------------------------------------------------+

This module contains the following functions:

    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Function          | Description                                                                                           |
    +===================+=======================================================================================================+
    | checkDataLength   | This function checks that the data at line *index* - 1 are consistence with *length*.                 |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | getLine           | This function gets the next line to analyze.                                                          |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | read              | This function reads in a RIS file.                                                                    |
    +-------------------+-------------------------------------------------------------------------------------------------------+
"""

import pathlib

from pqu import PQU as PQUModule
from fudge  import protareProductInfo as protareProductInfoModule

class RIS:
    """
    This class stores the contents of an RIS file. The contents of the file are stored as a list of :py:class:`Protare`
    instances inside a nested directionary. The top nested directionary's keys are the projectiles. Each value of 
    the projectile directionary is a directionary of the targets for that projectile. Each target's value is a directionary
    of the evaluations for that projectile/target combination. Finally, each evaluation's value is a :py:class:`Protare`
    instance.

    The following table list the primary members of this class:

    +---------------+-----------------------------------------------------------------------------------+
    | Member        | Description                                                                       |
    +===============+===================================================================================+
    | format        | The format of the RIS file.                                                       |
    +---------------+-----------------------------------------------------------------------------------+
    | projectiles   | This is a dictionary of the projectiles.                                          |
    +---------------+-----------------------------------------------------------------------------------+
    """

    def __init__(self, format):
        """
        :param format:     The format of the RIS file.
        """

        if format not in protareProductInfoModule.Formats.allowed:
            raise ValueError('Unsupported RIS format "%s"' % format)
        self.__format = format

        self.__projectiles = {}

    def __iter__(self):
        """
        This method iterates all the protares in *self*.

        :returns:                   An instance of :py:class`Protare`.
        """

        for projectileId in self.__projectiles:
            projectile = self.__projectiles[projectileId]
            for targetId in projectile:
                target = projectile[targetId]
                for evaluation in target:
                    yield target[evaluation]

    @property
    def projectiles(self):
        """
        This method returns a reference to *self*'s projectiles member.
        """

        return self.__projectiles

    @property
    def format(self):
        """
        This method returns a reference to *self*'s format member.
        """

        return self.__format

    def addProtare(self, protare):
        """
        This method adds *protare* to *self*'s list of protares.

        :param protare:     The :py:class`Protare` to add to *self*.
        """

        if not isinstance(protare, Protare):
            raise TypeError('Invalid instance: got type "%s".' % type(protare))

        if protare.projectileId not in self.__projectiles:
            self.__projectiles[protare.projectileId] = {}

        projectile = self.__projectiles[protare.projectileId]
        if protare.targetId not in projectile:
            projectile[protare.targetId] = {}

        target = projectile[protare.targetId]
        if protare.evaluation in target:
            print('WARNING: evalution "%s" already in "%s + %s" so skipping.' % (protare.evaluation, protare.projectileId, protare.targetId))
            return

        target[protare.evaluation] = protare

    def protare(self, projectileId, targetId, evaluation=None):
        """
        This method returns the RIS protare match the input arguments.

        :param projectileId:        The GNDs PoPs id of the projectile.
        :param targetId:            The GNDs PoPs id of the target.
        :param evaluation:          The evaluation string. If None, the first protare match the other arguments is returned.

        :returns:                   An instance of :py:class`Protare`.
        """

        if projectileId not in self.__projectiles:
            raise ValueError('Projectile "%s" not in RIS.' % projectileId)
        projectileRIS = self.__projectiles[projectileId]

        if targetId not in projectileRIS:
            raise ValueError('Target "%s" not in RIS.' % targetId)
        targetRIS = projectileRIS[targetId]

        if evaluation is None:                          # Get the first evaluation.
            for evaluation in targetRIS:
                protareRIS = targetRIS[evaluation]
                break
        else:
            protareRIS = targetRIS[evaluation]

        return protareRIS

    def print(self):
        """Prints the contents of *self* to the screen in a format similar to how the data are stored in a file."""

        print('#ris : 1.0')         # For format "1.*", the first six elements of each reaction data must be the sample as for "1.0".
        for projectile in self.__projectiles:
            targets = self.__projectiles[projectile]
            for target in targets:
                evaluations = targets[target]
                for evaluation in evaluations:
                    evaluations[evaluation].print()

    @staticmethod
    def read(path):
        """
        This static method reads in the RIS file specified by *path*.

        :param path:        The path to the RIS file to read.

        :returns:           An instance of :py:class`RIS`.
        """

        with open(path) as fIn:
            lines = fIn.readlines()

        index, command, data = getLine(0, lines, True, False)
        if command != '#ris':
            raise ValueError('File is not an RIS file.')
        format = checkDataLength(index, lines, 1, data)
        specialFormatFlag = ' : ' in lines[0]

        ris = RIS(format)

        numberOfLines = len(lines)
        while index < numberOfLines:
            index, command, data = getLine(index, lines, True, True)
            if command == '#import':
                importPath = checkDataLength(index, lines, 1, data)
                importPath2 = pathlib.Path(path).parent / importPath
                importedRIS = read(importPath2)
                counter = 0
                for protare in importedRIS:
                    counter += 1
                    ris.addProtare(protare)
            elif command == '#protare':
                projectileId, targetId, evaluation, energyUnit = checkDataLength(index, lines, 4, data)
                protare = Protare(projectileId, targetId, evaluation, energyUnit)
                index = protare.parse(index, lines, specialFormatFlag)
                ris.addProtare(protare)
            else:
                raise ValueError('Invalid command "%s".' % command)

        return ris

class Protare:
    """
    This class represents one protare in an RIS file.

    The following table list the primary members of this class:

    +---------------+-----------------------------------------------------------------------------------+
    | Member        | Description                                                                       |
    +===============+===================================================================================+
    | projectileId  | This is the GNDS PoPs id for the projectile of the protare.                       |
    +---------------+-----------------------------------------------------------------------------------+
    | targetId      | This is the GNDS PoPs id for the target of the protare.                           |
    +---------------+-----------------------------------------------------------------------------------+
    | evaluation    | This is the evaluation string for the protare.                                    |
    +---------------+-----------------------------------------------------------------------------------+
    | energyUnit    | This is the energy unit used in the protare.                                      |
    +---------------+-----------------------------------------------------------------------------------+
    | aliases       | This is a dictionary of the aliases defined in the protare.                       |
    +---------------+-----------------------------------------------------------------------------------+
    | reactions     | This is a python list of the reactions as :py:class:`Reaction` instances.         |
    +---------------+-----------------------------------------------------------------------------------+
    """

    def __init__(self, projectileId, targetId, evaluation, energyUnit):
        """
        :param projectileId:        The GNDs PoPs id of the projectile.
        :param targetId:            The GNDs PoPs id of the target.
        :param evaluation:          The evaluation string. If None, the first protare match the other arguments is returned.
        :param energyUnit:          The energy of the threshold for each reaction.
        """

        self.__projectileId = projectileId
        self.__targetId = targetId
        self.__evaluation = evaluation
        self.__energyUnit = energyUnit
        self.__aliases = {}
        self.__reactions = []

    def __str__(self):
        """
        This method returns a simple string representation of *self*.

        :returns:       A python str instance.
        """

        return '%s %s %s %s' % (self.__projectileId, self.__targetId, self.__evaluation, self.__energyUnit)

    def __iter__(self):
        """
        Iterator over the list of reactions.
        """

        for reaction in self.__reactions:
            yield reaction

    @property
    def projectileId(self):
        """
        This method returns a reference to *self*'s projectileId member.
        """

        return self.__projectileId

    @property
    def targetId(self):
        """
        This method returns a reference to *self*'s targetId member.
        """

        return self.__targetId

    @property
    def evaluation(self):
        """
        This method returns a reference to *self*'s evaluation member.
        """

        return self.__evaluation

    @property
    def energyUnit(self):
        """
        This method returns a reference to *self*'s energyUnit member.
        """

        return self.__energyUnit

    @property
    def aliases(self):
        """
        This method returns a reference to *self*'s aliases member.
        """

        return self.__aliases

    @property
    def reactions(self):
        """
        This method returns a reference to *self*'s reactions member.
        """

        return self.__reactions

    def parse(self, index, lines, specialFormatFlag):
        """
        This method parse all the data for a protare in an RIS file after the command line of the protare.

        :param index:               The index of the current line in *lines* to process.
        :param lines:               All the lines in the RIS file.
        :param specialFormatFlag:   If **True**, the data separator is ' : ', otherwise it is ': '.

        :returns:               The index of the next line to processed.
        """

        numberOfLines = len(lines)
        while index < numberOfLines:
            index, command, data = getLine(index, lines, True, True)
            if   command == '#aliases':
                numberOfAliases = int(checkDataLength(index, lines, 1, data))
                for index2 in range(numberOfAliases):
                    index, command, data = getLine(index, lines, False, True)
                    aid, pid = checkDataLength(index, lines, 2, data)
                    self.__aliases[aid] = pid
            elif command == '#reactions':
                numberOfReactions = int(checkDataLength(index, lines, 1, data))
                for index2 in range(numberOfReactions):
                    index, command, data = getLine(index, lines, False, specialFormatFlag)
                    if len(data) == 4:
                        label, threshold, intermediate, process = checkDataLength(index, lines, 4, data)
                        reactionLabel, covarianceFlag = None, None
                    else:
                        label, threshold, intermediate, process, reactionLabel, covarianceFlag = checkDataLength(index, lines, 6, data)
                    threshold = PQUModule.PQU(float(threshold), self.__energyUnit)
                    reaction = Reaction(label, threshold, intermediate, process, reactionLabel, covarianceFlag)
                    self.__reactions.append(reaction)
            else:
                break

        if index < numberOfLines:
            index -= 1

        return index

    def print(self):
        """Prints the contents of *self* to the screen in a format similar to how the data are stored in a file."""

        print('#protare : %s : %s : %s : %s' % (self.__projectileId, self.__targetId, self.__evaluation, self.__energyUnit))

        if len(self.__aliases) > 0:
            print('#aliases : %s' % len(self.__aliases))
            for key, value  in self.__aliases.items():
                print('    %s : %s' % (key, value))

        print('#reactions : %s' % len(self.__reactions))
        processWidth = 9
        for reaction in self.__reactions:
            if reaction.process is not None:
                if reaction.process != 'totalFission':
                    processWidth = max(processWidth, len(reaction.process))
        processFormat = '%%-%ds' % processWidth
        reactionLabelWidth = max([len(reaction.reactionLabel) for reaction in self.__reactions])
        reactionLabelFormat = '%%-%ds' % reactionLabelWidth

        for reaction in self.__reactions:
            reaction.print(processFormat, reactionLabelFormat)

class Reaction:
    """
    This class represents one reaction in an RIS file.

    The following table list the primary members of this class:

    +-------------------+-----------------------------------------------------------------------------------+
    | Member            | Description                                                                       |
    +===================+===================================================================================+
    | label             | This the label for the reaction.                                                  |
    +-------------------+-----------------------------------------------------------------------------------+
    | threshold         | This is the effective threshold for the reaction.                                 |
    +-------------------+-----------------------------------------------------------------------------------+
    | intermediate      | This is the intermediate product for the reaction.                                |
    +-------------------+-----------------------------------------------------------------------------------+
    | process           | This is the process label for the reaction.                                       |
    +-------------------+-----------------------------------------------------------------------------------+
    | reactionLabel     | This is the GNDS label for the reaction.                                          |
    +-------------------+-----------------------------------------------------------------------------------+
    | covarianceFlag    | This is the string 'covariance' if the reaction has covariance and None otherwise.|
    +-------------------+-----------------------------------------------------------------------------------+
    """

    def __init__(self, label, threshold, intermediate, process, reactionLabel=None, covarianceFlag=None):
        """
        :param label:               The label for the reaction.
        :param threshold:           The threshold for the reaction.
        :param intermediate:        The intermediate product for the reaction.
        :param process:             The process string for the reaction.
        """

        self.__label = label
        self.__threshold = threshold
        self.__intermediate = intermediate
        self.__process = process
        self.__reactionLabel = reactionLabel
        self.__covarianceFlag = covarianceFlag

    def __str__(self):
        """
        This method returns a simple string representation of *self*.

        :returns:       A python str instance.
        """

        processStr = ''
        if len(self.__process) > 0:
            processStr = ' [%s]' % self.__process

        return '%s [%s] %s' % (self.__label, processStr, self.__threshold)

    @property
    def label(self):
        """
        This method returns a reference to *self*'s label member.
        """

        return self.__label

    @property
    def threshold(self):
        """
        This method returns a reference to *self*'s threshold member.
        """

        return self.__threshold

    @property
    def intermediate(self):
        """
        This method returns a reference to *self*'s intermediate member.
        """

        return self.__intermediate

    @property
    def process(self):
        """
        This method returns a reference to *self*'s process member.
        """

        return self.__process

    @property
    def reactionLabel(self):
        """
        This method returns a reference to *self*'s reactionLabel member.
        """

        return self.__reactionLabel

    @property
    def covarianceFlag(self):
        """
        This method returns a reference to *self*'s covarianceFlag member.
        """

        return self.__covarianceFlag

    def print(self, processFormat, reactionLabelFormat):
        """Prints the contents of *self* to the screen in a format similar to how the data are stored in a file."""

        threshold = float(PQUModule.floatToShortestString(self.__threshold))
        process = processFormat % self.__process
        reactionLabel = reactionLabelFormat % self.__reactionLabel
        print('    %-31s : %-12s : %-10s : %s : %s : %s' % 
                (self.__label, threshold, self.__intermediate, process, reactionLabel, self.__covarianceFlag))

def checkDataLength(index, lines, length, data):
    """
    For internal use only. This function checks that the data at line *index* - 1 are consistence with *length*.

    :param index:           The index of the current line in *lines* to process.
    :param lines:           All the lines in the RIS file.
    :param length:          The required number of item in *data*.
    :param data:            A python list of entries on a RIS line excluding the command entry if present.

    :returns:               The single datum if *length* is 1 and python list *data* otherwise.
    """

    if len(data) != length:
        raise Exception('Wrong number of data at index (%s): %s vs %s: "%s".' % (index, len(data), length, lines[index-1][:-1]))

    if length == 1:
        data = data[0]

    return data

def getLine(index, lines, commandLine, specialFormatFlag):
    """
    For internal use only. This function gets the next line in *lines* at index *index* and returns the next index,
    the comnnad if *commandLine* is True and the parsed data on the line.

    :param index:               The index of the current line in *lines* to process.
    :param lines:               All the lines in the RIS file.
    :param commandLine:         If True, line at index is a command line. Otherwise it is not a command line.
    :param specialFormatFlag:   If **True**, the data separator is ' : ', otherwise it is ': '.

    :returns:               The tuple index, command, data.
    """

    sep = ' : ' if specialFormatFlag else ': '
    line = lines[index]
    data = line.split(sep)
    data = [datum.strip() for datum in data]
    if commandLine:
        command = data.pop(0)
    else:
        command = None

    index += 1

    return index, command, data

def read(path):
    """
    This function reads in the RIS file specified by *path*.

    :param path:        The path to the RIS file to read.

    :returns:           An instance of :py:class`RIS`.
    """

    return RIS.read(path)
