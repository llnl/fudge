#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module contains functions for creating an ris file from a GNDS map file. The main user function is **buildReactionInfoSummary**.
This module can also be run from the command line and uses the module **argparse** to parse input arguments.
'''

import sys
import pathlib
import argparse
import multiprocessing

from LUPY import subprocessing as subprocessingModule

from fudge import enums as enumsModule
from fudge import map as mapModule
from fudge import protareProductInfo as protareProductInfoModule
from fudge import reactionSuite as reactionSuiteModule

summaryDocString__FUDGE = """Builds **ris** (reaction infomation summary) data for each protare in the specified map file and writes the results to a **ris** file."""

description = summaryDocString__FUDGE

def worker(inputQueue, outputQueue):

    for function, args in iter(inputQueue.get, 'STOP'):
        result = function(*args)
        outputQueue.put(result)

def processImport(index, path):
    '''
    Processes an import node.

    :param index:           The index used to do the final sorting so that the order in the ris file is the same as in the map file.
    :param path:            The path specified in the import node.
    '''

    return '%s___%s' % (index, '#import : %s\n' % pathlib.Path(path).with_suffix('.ris'))

def processProtare(index, fileName):
    '''
    Processes a protare node which points to a GNDS reactionSuite.

    :param index:           The index used to do the final sorting so that the order in the ris file is the same as in the map file.
    :param fileName:        The path to the GNDS reactionSuite file.
    '''

    lines = []
    try:
        protare = reactionSuiteModule.read(fileName, lazyParsing=True)
    except:
        return '%s___%s in buildReactionInfoSummary.py for file %s' % (index, 'ERROR', fileName)

    lines.append('#protare : %s : %s : %s : %s' % (protare.projectile, protare.target, protare.evaluation, protare.domainUnit))

    aliases = [[alias.id, protare.PoPs.final(alias.id).id] for alias in protare.PoPs.aliases if alias.isMetaStable()]
    if len(aliases) > 0:
        lines.append('#aliases : %s' % len(aliases))
        for pid, alias in aliases:
            lines.append('    %s : %s' % (pid, alias))

    lines.append('#reactions : %s' % len(protare.reactions))
    lines += protareProductInfoModule.protareProductInfo2(protare, '    ') + ['']

    return '%s___%s' % (index, '\n'.join(lines))

def buildReactionInfoSummary(input, recursive=False, numberOfProcesses=None):
    '''
    This function reads the map file specified by *input* and creates it ris file. This function uses the python module
    multiprocessing to parallelize the processing of all the child nodes nodes.

    :param input:               The name of the map file to process.
    :param recursive:           If **True**, each map file specified by an import child node is also processed.
    :param numberOfProcesses:   The number of multiprocessing.Process to use.
    '''

    if numberOfProcesses is None:
        numberOfProcesses = multiprocessing.cpu_count()

    map = mapModule.read(input, lazyParsing=True)

    inputQueue  = multiprocessing.Queue()
    outputQueue = multiprocessing.Queue()

    recursiveMapNames = []
    counter = 0
    for index, entry in enumerate(map):
        if isinstance(entry, mapModule.Import):
            inputQueue.put([processImport, (index, entry.path)])
            recursiveMapNames.append(entry.buildFileName())
        else:
            if entry.interaction in (enumsModule.Interaction.atomic, enumsModule.Interaction.TNSL):
                continue
            inputQueue.put([processProtare, (index, entry.buildFileName())])
        counter += 1

    numberOfProcesses = min(counter, max(1, numberOfProcesses))
    processes = []
    for index in range(numberOfProcesses):
        process = multiprocessing.Process(target=worker, args=(inputQueue, outputQueue, ))
        process.start()
        processes.append(process)

    for index in range(numberOfProcesses):
        inputQueue.put('STOP')

    results = {}
    for queueIndex in range(counter):
        result = outputQueue.get()
        offset = result.find('___')
        index = int(result[:offset])
        data = result[offset+3:]
        if data[:5] == 'ERROR':
            print(data)
            data = ''
        results[index] = data

    input = pathlib.Path(map.path)
    output = args.output
    if output is None:
        output = input.with_suffix('.ris')
    with output.open('w') as fOut:
        fOut.write('#ris : 1.0\n')
        for index in sorted(results):
            fOut.write(results[index])

    if recursive:
        for mapName in recursiveMapNames:
            subprocessingModule.executeCommand([sys.executable, __file__, '--recursive', mapName])

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('input', type=pathlib.Path,                         help='The map file whose protares are analyzed.')
    parser.add_argument('output', type=pathlib.Path, default=None, nargs='?',
                                                                            help='The path of the outputted RIS file. Default is to use input map file with extension ".ris".')
    parser.add_argument('--recursive', action='store_true',                 help='If present, imported map files a also processed; otherwise, import map files are not processed.')
    parser.add_argument('-n', '--numberOfProcesses', action='store', type=int, default=multiprocessing.cpu_count(),
                                                                            help='Limits the number of child process to it value. Default is the node number of CPUs.')

    args = parser.parse_args()

    buildReactionInfoSummary(args.input, recursive=args.recursive, numberOfProcesses=args.numberOfProcesses)
