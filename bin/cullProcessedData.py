#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib

from LUPY import argumentsForScripts as argumentsForScriptsModule
import fudge.styles as stylesModule

summaryDocString__FUDGE = """Removes all processed data (except optionally reconstructed resonances) from the specified GNDS reactionSuite."""

description = summaryDocString__FUDGE + """ If option "--outputDir" is present, the directory is stripped from the output path and replaced
with the specified path."""

parser = argparse.ArgumentParser(description=description, allow_abbrev=False)
singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)
parser.add_argument('outputPath', nargs='?', default=None, type=pathlib.Path,  help='Output file name.')
parser.add_argument('-o', '--outputDir', default=None, type=pathlib.Path,      help='The output directory to write the output file to.')
parser.add_argument('-s', '--suffix', default='.culled.xml', type=str,         help='Suffix to add to culled file.')
parser.add_argument('-r', '--reconstructedResonances', action='store_true',    help='Also cull reconstructed resonances.')

args = parser.parse_args()

protare = singleProtareArguments.protare(args, lazyParsing=False)
protare.cullProcessedData(args.reconstructedResonances)

outputPath = args.outputPath
if outputPath is None:
    outputPath = pathlib.Path(protare.sourcePath).with_suffix(args.suffix)
if args.outputDir is not None:
    outputPath = args.outputDir / outputPath.name

protare.saveAllToFile(outputPath)
