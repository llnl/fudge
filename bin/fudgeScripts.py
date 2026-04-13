#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = '''
This script examines all "*.py" files in the FUDGE bin directory and prints a summary about each.
'''

import re
import pathlib
import argparse
import textwrap
import shutil

maxSpecialVariableLineIndex = 100       # Special variable definition must come before this line or script is ignored.
numberOfColumns = shutil.get_terminal_size().columns

parser = argparse.ArgumentParser(description=description)
args = parser.parse_args()

def addScriptsInBin(binPath, summaryDocString):

    for scriptPath in binPath.glob('*.py'):
        if scriptPath.name == thisScriptsName:
            continue
        with scriptPath.open('r') as fIn:
            lines = fIn.readlines()
            for index, line in enumerate(lines):
                if index > maxSpecialVariableLineIndex:     # Extract check in case there are non-FUDGE scripts in bin directory.
                    break
                    print('        ', re.match(r'^summaryDocString__\s*=\s*.+$', line))
                if re.match(r'^summaryDocString__[a-zA-Z]*\s*=\s*.+$', line):
                    variable = line.split('=')[0].strip()
                    locals = {}
                    exec(line, None, locals)  # locals is only supported as kwarg starting in 3.13
                    subModule = line.split('=')[0].strip().split('__')[1]
                    if subModule not in summaryDocString:
                        summaryDocString[subModule] = {}
                    summaryDocString[subModule][scriptPath.name] = locals[variable]
                    break

def printInfo(moduleName, scripts):
    r'''Prints the list of all python files in *scripts* that contain the line "'^summaryDocString__%s\s*=\s*.+$' % moduleName".'''

    print()
    print('Scripts in %s:' % moduleName)
    nameWidth = 0
    for name in scripts:
        nameWidth = max(nameWidth, len(name))
    format = '%%-%ds - %%s' % nameWidth

    indent = 4 * ' '
    indent2 = (nameWidth + 3) * ' ' + indent

    for name in sorted(scripts):
        line = format % (name, scripts[name])
        print('\n'.join(textwrap.wrap(line, width=numberOfColumns,  initial_indent=indent, subsequent_indent=indent2)))

binDir = pathlib.Path(__file__).resolve().parent
thisScriptsName = pathlib.Path(__file__).name

scripts = sorted(binDir.glob('*.py'))
subModules = {}
for binPath in binDir.parent.glob('**/bin'):
    addScriptsInBin(binPath, subModules)

for subModule in ['FUDGE', 'PoPs', 'xData', 'Brownies']:
    if subModule not in subModules:
        continue
    printInfo(subModule, subModules.pop(subModule))

for subModule in subModules:
    printInfo(subModule, subModules[subModule])
