#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import argparse
import pathlib

from LUPY import times as timesModule
from LUPY import subprocessing as subprocessingModule

from PoPs import IDs as PoPsIDsModule
from PoPs import alias as PoPsAliasModule
from PoPs.chemicalElements import misc as PoPsGroupsMiscModule

from fudge import enums as enumsModule
from fudge import map as mapModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import suites as suitesModule
from fudge import GNDS_file as GNDS_fileModule

summaryDocString__FUDGE = '''Creates a map file from a list of GNDS reactionSuite and map files.'''

description = '''
This module creates a GNDS map file from a list of GNDS files. Required options are '--path', '--library'. The option
'--interaction' is also required if any GNDS file is stored as GNDS format 1.10. It should not be used for GNDS format 
2.0 or higher; in fact, it is ignored for files with format 2.0 or higher.  Currently, interactions "nuclear", "atomic" 
and "TNSL" cannot be mixed for GNDS format 1.10. That is, if interaction is 'nuclear', all GNDS files must be 'nuclear'
protare files, if interaction is 'atomic', all GNDS files must be 'atomic' protare files, and if interaction is 'TNSL',
all GNDS files must be 'TNSL' protare files. For format 2.0 or higher, mixing of interactions is supported.

In addition to building a map file from a list of protare files, the list of input files can be other map files. In this
case, the built map file will contain a list of import nodes.

Suppose one want to create 4 map files with the following contents:

    neutrons.noMetastables.map  : All protares in the **neutrons** sub-directory excluding the ones with meta-stable targets.
    neutrons.map                : Importing **neutrons.noMetastables.map** plus all protares in the **neutrons** sub-directory
                                    that have meta-stables targets.
    all.noMetastables.map       : Importing the neutrons.noMetastables.map, protons.map, deuterons.map, tritons.map, 
                                    helium3s.map, alphas.map, photo_atomic.map and tslmap.map map files.
    all.map                     : Importing the neutrons.map, protons.map, deuterons.map, tritons.map, 

The following show how one might use this python script to create the 4 map files:

    otherMaps="protons.map deuterons.map tritons.map helium3s.map alphas.map photo_atomic.map tsl.map"
    buildMapFile.py -l neutrons -p neutrons.noMetastables.map --ignoreMetaStables                                    neutrons/*.xml
    buildMapFile.py -l neutrons -p neutrons.map               --nonMetaStablesMapFile neutrons.noMetastables.map neutrons/*.xml
    buildMapFile.py -l endl2009.4-rc2 -p all.noMetastables.map neutrons.noMetastables.map $otherMaps
    buildMapFile.py -l endl2009.4-rc2 -p all.map               neutrons.map               $otherMaps

For a TNSL map node, standardTarget and standardEvaluation attributes are needed. If the option '--standards' is specified,
the next argument is a json file that specifies the standardTarget, or the standardTarget and standardEvaluation for
each TNSL target key listed. This file is read in with the json python module and must contain a json dictionary.  The 
dictionary keys are TNSL target names (e.g., 'Be-metal', 'HinH2O'). Each key's associated value is either the standard 
target (e.g., 'Be9', 'H1'), or a list containing the standard target and standard evaluation (e.g., [ 'H1', 'ENDF/B-VII.2' ] )
for the TNSL target. If no file is given or no key exists for a TNSL target, its standardTarget attribute is left empty 
(i.e., '').  If no standard evaluation exists for a TNSL target, the standard evaluation is taken from the TNSL target's
evaluation.

An example of a --standards file is:

    {   "HinCH2"                    : "H1",
        "HinH2O"                    : "H1",
        "HinYH2"                    : [ "H1", "ENDF/B-VII.1" ],
        "DinD2O"                    : "H2",
        "OinD2O"                    : "O16",
        "YinYH2"                    : "Y89" }

Sometime the map file being built cannot be placed (via the *--path* option) in its final location. For example, 
one does not have permission to write the map file to the final location. If one is not putting the map file in 
its final location, and checksums are being added, then one must use the *--mapDirectory* to specify the location 
of the final map file. Here is an example where the final location of the map file will ve in the current directory
but it is initially put into the user's home directory:

    /path/to/buildMapFile.py neutrons/* -o ~/neutron.map -l lib11 --mapDirectory ./
'''

interactionChoices = enumsModule.Interaction.allowedStrings()
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, nargs='*',                               help='A list of GNDS and map files to put into a map file.')
parser.add_argument('-o', '-p', '--path', required=True,                        help='The path/name of the built map file.')
parser.add_argument('-l', '--library', required=True,                           help='The library name to use for the map file.')
parser.add_argument('-i', '--interaction', action='store', choices=interactionChoices, default=None,
                                                                                help='Specifies the interaction type of the GNDS files to map. Currently, only one interaction type can be \nhandled per execution. This options is not needed for GNDS 2.0 and will be deprecated.')
parser.add_argument('--format', action='store', choices=mapModule.FormatVersion.allowed, default=mapModule.FormatVersion.default,
                                                                                help='Specifies the format for the map file. Default is "%s".' % mapModule.FormatVersion.default)
parser.add_argument('--standards', action='store',                              help='Specifies a json file that is a dictionary with TNSL target names as keys, and each with an associated \nvalue of a string "standardTargets" or a list of "standardTargets" and "standardEvaluation".')
parser.add_argument('--ignoreMetaStables', action='store_true',                 help='If present, meta-stable targets in the "files" argument are ignored. Also see option *--nonMetaStablesMapFile*.')
parser.add_argument('--nonMetaStablesMapFile', action='store', default=None,    help='If present, the option "--ignoreMetaStables" is ignored and instead, only meta-stable targets are \nincluded in the map file. The next argument should be the non-metastable map file which will be included \nas an "import" node at the beginning of the generated map file.')
parser.add_argument('--skipChecksums', action='store_true',                     help='If present, *checksums* are not added to the map file.')
parser.add_argument('--skipRIS', action='store_true',                           help='If present, an *ris* is not create for the map file.')
parser.add_argument('--mapDirectory', action='store', default=None,             help='If the outputted map file will not reside in its final location (i.e., directory), then use this option to specify its final location. Typically, one should enter this as "--mapDirectory ./".')
parser.add_argument('-v', '--verbose', action='count', default=0,               help='Determines how verbose %(prog)s will be.')

args = parser.parse_args()

def particleSortIndex(name):
    '''
    Returns a tuple of length 5 integers with the following meanings:

        +-----------+---------------------------------------------------+
        | index     | description                                       |
        +-----------+---------------------------------------------------+
        | 0         | Type of particle.                                 |
        |           |   0 for a nuclide or a neutron                    |
        |           |   1 for a photon                                  |
        |           |   2 for an electron                               |
        +-----------+---------------------------------------------------+
        | 1         | If index 0 is 0 this represents the particle's Z. |
        |           | Otherwise it is 0.                                |
        +-----------+---------------------------------------------------+
        | 2         | If index is 0 this represents the particle's A.   |
        |           | Otherwise it is 0.                                |
        +-----------+---------------------------------------------------+
        | 3         | If index 0 is 0 this represents the nuclide's     |
        |           | nuclear excitation level index. Otherwise it is 0.|
        +-----------+---------------------------------------------------+
        | 4         | If index 0 is 0 this represents the nuclide's     |
        |           | metastable index.  Otherwise it is 0.             |
        +-----------+---------------------------------------------------+
    '''

    nuclideName, metaStableLevel = PoPsAliasModule.MetaStable.nuclideNameAndMetaStableIndexFromName(name)
    if metaStableLevel != 0:
        name = nuclideName

    if name == PoPsIDsModule.neutron:
        return 0, 0, 1, 0, 0
    elif name == PoPsIDsModule.photon:
        return 1, 0, 0, 0, 0
    elif name == PoPsIDsModule.electron:
        return 2, 0, 0, 0, 0

    info = list(PoPsGroupsMiscModule.chemicalElementALevelIDsAndAnti(str(name)))
    if args.verbose > 4:
        print('    ', name, info)

    if info[1] is None:             # This is the chemical symbol.
        return 1, 0, 0, 0, 0
    if info[2] is None:             # This is A.
        info[2] = 0
    if info[3] is None:             # This is nuclear level index.
        info[3] = 0

    return 0, PoPsGroupsMiscModule.ZFromSymbol[info[1]], info[2], info[3], metaStableLevel

def interactionSortIndex(interaction):

    if interaction is None:
        interaction = args.interaction

    interaction = enumsModule.Interaction.checkEnumOrString(interaction)
    if interaction == enumsModule.Interaction.nuclear:
        return 0
    if interaction == enumsModule.Interaction.TNSL:
        return 1
    if interaction == enumsModule.Interaction.atomic:
        return 2
    if interaction == enumsModule.Interaction.LLNL_TNSL:
        return 3

    return 4

def sortFiles(files):

    _files = []

    for file, data in files:
        if args.verbose > 4:
            print(file, data)
        interaction = interactionSortIndex(data['interaction'])
        projectile = particleSortIndex(data['projectile'])
        target = particleSortIndex(data['target'])
        _files.append((interaction, projectile, target, file, data ))

    return sorted(_files)

standards = {}
if args.standards is not None:
    import json
    with open( args.standards, 'r') as f :
        standards = json.load( f )
    for target in standards :
        if isinstance( standards[target], list ): continue
        standards[target] = [ standards[target], None ]

mode = 'all'
if args.nonMetaStablesMapFile is not None:
    mode = 'metastables'
elif args.ignoreMetaStables:
    mode = 'non-metastables'

fileType = None
groups = []
for file in args.files :
    try :
        name, data = GNDS_fileModule.type(file)
    except :
        print('    WARNING: Invalid file "%s".' % file)
        continue
    if name == reactionSuiteModule.ReactionSuite.moniker:
        if data['interaction'] is None:
            if data['projectile'] == PoPsIDsModule.photon:
                raise Exception('File %s does not have interaction attribute.' % file)
            data['interaction'] = enumsModule.Interaction.nuclear
        pass
    elif name == mapModule.Map.moniker:
        pass
    else :
        if args.verbose > 2: print('    WARNING: Ignoring file of type %s.' % name)
        continue
    if fileType != name:
        if fileType is not None: groups.append( [ fileType, files ] )
        fileType = name
        files = []
    files.append( [ file, data ] )
if fileType is not None: groups.append( [ fileType, files ] )

map = mapModule.Map( args.library, args.path )
if args.nonMetaStablesMapFile is not None:                              # ensure file exists and is a map file
    fileType, info = GNDS_fileModule.type(args.nonMetaStablesMapFile)
    assert fileType == 'map', "nonMetaStablesMapFile argument must be a map file! Received %s instead" % fileType
    map.append(mapModule.Import(args.nonMetaStablesMapFile))

TNSL_missingStandardTargets = []
for fileType, files in groups :
    if fileType == mapModule.Map.moniker:
        for file, data in files : map.append( mapModule.Import( file ) )
    else :
        for interaction, projectile, target, file, data in sortFiles( files ) :
            if mode == 'non-metastables' and target[-1] != 0: continue
            if mode == 'metastables'     and target[-1] == 0: continue
            interaction = data['interaction']
            if interaction is None:
                protare = GNDS_fileModule.preview(file, haltParsingMoniker=None)
                interaction = protare.interaction
            interaction = enumsModule.Interaction.checkEnumOrString(interaction)
            if interaction == enumsModule.Interaction.legacyTNSL:
                interaction = enumsModule.Interaction.TNSL
            if interaction == enumsModule.Interaction.TNSL:
                # Check if the file has targetInfo at /reactionSuite/styles/evaluated/targetInfo
                protare = GNDS_fileModule.preview(file, haltParsingMoniker=suitesModule.Reactions.moniker)
                hasTargetInfo = False
                evaluated_style = protare.styles.getEvaluatedStyle()
                if hasattr(evaluated_style, 'targetInfo') and evaluated_style.targetInfo is not None:
                    hasTargetInfo = True

                if hasTargetInfo:
                    # Skip standardTarget and standardEvaluation if targetInfo exists
                    map.append( mapModule.TNSL( str( data['projectile'] ), str( data['target'] ), str( data['evaluation'] ), file,
                            None, None, interaction ) )
                else:
                    try :
                        standardTarget = standards[data['target']][0]
                    except :
                        standardTarget = ''
                    if standardTarget == '': TNSL_missingStandardTargets.append( data['target'] )

                    try :
                        standardEvaluation = standards[data['target']][1]
                        if standardEvaluation is None: standardEvaluation = str( data['evaluation'] )
                    except :
                        standardEvaluation = str( data['evaluation'] )

                    map.append( mapModule.TNSL( str( data['projectile'] ), str( data['target'] ), str( data['evaluation'] ), file,
                            standardTarget, standardEvaluation, interaction ) )
            else :
                map.append( mapModule.Protare( str( data['projectile'] ), str( data['target'] ), str( data['evaluation'] ), file, interaction ) )

if not args.skipChecksums:
    if args.verbose > 0: print('    INFO: Calculating checksums', end='')
    checksumTime = timesModule.Times()
    map.updateAllChecksums(mapDirectory=args.mapDirectory)
    checksumTime = checksumTime.delta()
    if args.verbose > 0: print(': cpu %.2f s, wall %.2f s, sys %.2f s' % (checksumTime['user'], checksumTime['wall'], checksumTime['sys']))

map.saveToFile(args.path, format=args.format)

if not args.skipRIS:
    if args.verbose > 0: print('    INFO: building RIS file with buildReactionInfoSummary.py')
    RIS_file = pathlib.Path(__file__).parent / 'buildReactionInfoSummary.py'
    subprocessingModule.executeCommand((sys.executable, RIS_file, args.path), stdout=sys.stdout, stderr=sys.stderr)

if len( TNSL_missingStandardTargets ) > 0:
    print( '    WARNING: you still need to update the standard target (and maybe the standard evaluation) for %s TNSL nodes.' % len( TNSL_missingStandardTargets ) )
    if args.verbose > 1:
        for TNSL_missingStandardTarget in TNSL_missingStandardTargets : print('        Standard target missing for %s' % TNSL_missingStandardTarget)
