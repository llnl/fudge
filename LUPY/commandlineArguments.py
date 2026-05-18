# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains a function for getting the list of arguments from a :py:class:`argparse.ArgumentParser` instance
and two functions for reading and writing data to the computerCodes member of a :py:class:`documentationModule.Documentation` instance.
"""

import argparse
import getpass
import socket

import fudge
from xData.Documentation import documentation as documentationModule
from PoPs import IDs as IDsModule

from LUPY import argumentsForScripts as argumentsForScriptsModule

from xData import date as xDataDateModule
from xData.Documentation import computerCode as computerCodeModule
from xData.Documentation import dates as datesModule

def getArgparseArguments(parser, args, excludeArguments=[]):
    """
    This function returns command line arguments of an :py:class:`argparse.ArgumentParser` instance as a list of 
    positional arguments and a list of optional arguments, excluding any arguments listed in *excludeArguments*.

    :param parser:              A :py:class:`argparse.ArgumentsParser` instance.
    :param args:                Populated namespace returned by the :py:func:`argparse.ArgumentParser.parse_args` method.
        :type args:             argparse.Namespace
    :param excludeArguments:    List of names (i.e., arguments) in namespace *args* to ignore.

    :return:                    The tuple (positionalArguments, optionalArguments).
    """

    if not isinstance(parser, argparse.ArgumentParser):
        raise RuntimeError(f'First argument expected to be of type {argparse.ArgumentParser}')

    if not isinstance(args, argparse.Namespace):
        raise RuntimeError(f'Second argument expected to be of type {argparse.Namespace}')

    if not isinstance(excludeArguments, (list, tuple)):
        raise RuntimeError('The "excludeArguments" argument is expected to be a list of command line arguments to ignore')

    parserActions = dict([(x.dest, x) for x in parser._actions])
    positionalArguments = []
    optionalArguments = []

    for arg in vars(args):
        if arg in excludeArguments:
            continue

        instance = parserActions[arg]
        values = getattr(args, instance.dest)
        if values is None:
            continue

        if len(parserActions[arg].option_strings) == 0:                                         # Positional arguments.
            if isinstance(values, (list, tuple)):
                for value in values:
                    positionalArguments.append(value)
            else:
                positionalArguments.append(values)
        elif isinstance(instance, argparse._StoreAction):                                       # store
            optionalArguments.append('%s %s' % (instance.option_strings[0], values))
        elif isinstance(instance, argparse._StoreConstAction):                                  # store_const
            optionalArguments.append(instance.option_strings[0])
        elif isinstance(instance, (argparse._StoreTrueAction, argparse._StoreFalseAction)):     # store_true or store_false
            if value != instance.default:
                optionalArguments.append(instance.option_strings[0])
        elif isinstance(instance, argparse._AppendAction):                                      # append
            for value in values:
                optionalArguments.append('%s %s' % (instance.option_strings[0], value))
        elif isinstance(instance, argparse._AppendConstAction):                                 # append_const
            raise ValueError('Action "append_const" is not supported.')
        elif isinstance(instance, argparse._CountAction):                                       # count
            optionalArguments += values * [instance.option_strings[0]]
        elif isinstance(instance, argparse._VersionAction):                                     # version
            raise ValueError('This if should always be false.')
        else:
            if hasattr(argparse, '_ExtendAction'):                                              # extend: only since python 3.8
                if isinstance(instance, argparse._ExtendAction):
                    for value in values:
                        optionalArguments.append('%s %s' % (instance.option_strings[0], value))
                    continue
            raise ValueError('Unsupported action.')

    return positionalArguments, optionalArguments

def commandLineArgumentsToDocumentation(codeName, commandLineArguments, documentation, date=None):
    """
    This function writes a list of command line arguments and a note to the computerCodes member 
    of a :py:class:`documentationModule.Documentation` instance.

    :param codeName:            String representing the computer code's name.
    :param currentArguments:    List of command line arguments.
    :param documentation:       :py:class:`xData.Documentation.documentation.Documentation` instance whose child nodes are populated.
    :param date:                Date code was run which is added to the dates member of *documentation*.
        :type date:             :py:class:`xData.date.Date`, optional
    """

    if not isinstance(commandLineArguments, (tuple, list)):
        raise TypeError('First argument expected to be a list of command line arguments.')

    if not isinstance(documentation, documentationModule.Documentation):
        raise TypeError(f'Second argument expected expected to be of type {documentationModule.Documentation}.')

    if codeName in documentation.computerCodes:
        raise KeyError(f'Entry "{codeName}" already in documentation.computerCodes.')

    if date is None:
        date = xDataDateModule.Date()
    if not isinstance(date, xDataDateModule.Date):
        raise TypeError('Date must be a {xDataDateModule.Date} instance.')
    date = datesModule.Date(date, dateType=datesModule.DateType.created)
    documentation.dates.add(date)

    computerCode = computerCodeModule.ComputerCode(codeName, codeName, fudge.__version__)
    computerCode.executionArguments.body = ' '.join(commandLineArguments)

    computerCode.note.body = f'Code {codeName} executed by {getpass.getuser()} on platform {socket.gethostname()}.'

    documentation.computerCodes.add(computerCode)

def commandLineArgumentsFromDocumentation(codeName, documentation):
    """
    This function reads command line arguments for computer code named *codeName* from a :py:class:`documentationModule.Documentation` instance.

    :param codeName:            String representing the computer code's name. Must also match the label of a **ComputerCode** 
                                instance in documentation.computerCodes.
    :param documentation:       A :py:class:`documentationModule.Documentation` instance.

    :return:                    The execution arguments for the code as a string. None is returned if *codeName* is not 
                                in documentation.computerCodes.
    """

    try:
        computerCode = documentation.computerCodes[codeName]
    except:
        return None

    return computerCode.executionArguments.body
