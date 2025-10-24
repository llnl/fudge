# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains functions and/or classes that are used by some of the FUDGE scripts.
"""

import os
import argparse
import pathlib

from PoPs import IDs as IDsModule

from fudge import enums as enumsModule
from fudge import map as mapModule
from fudge import GNDS_file as GNDS_fileModule
from fudge import reactionSuite as reactionSuiteModule

class SingleProtareArguments :
    """
    This class adds a positional argument and two option arguments to an :py:class:`argparse.ArgumentParser` instance 
    to give a common way to specify GNDS protare file to load as need by some of the FUDGE scripts.
    The positional argument can be a GNDS protare (i.e., reactionSuite) or map file path.
    If the positional argument is a map file path, then the optional argument "--pid" ("--tid") is used to specified
    the projectile (target) of the protare using a GNDS id.
    """

    def __init__( self, parser ) :
        """
        :param parser:      An :py:class:`argparse.ArgumentParser` instance which will be updated with 
                            parameters for accessing a profile from a file.
        """

        self.__parser = parser

        parser.add_argument( 'mapOrProtareFileName',                                    help = 'The path to a map or protare file. If it is a map file then option "--tid" (and maybe "--pid") is required.' )
        parser.add_argument( '--pid', action = 'store', default = IDsModule.neutron,    help = 'The PoPs id for the projectile. Only used if the argument "mapOrProtareFileName" is a map file; otherwise, it is ignored. Default is "%s".' % IDsModule.neutron )
        parser.add_argument( '--tid', action = 'store', default = None,                 help = 'The PoPs id for the target. Required if the argument "mapOrProtareFileName" is a map file; otherwise, it is ignored.' )
        parser.add_argument( '--noLazyParse', action = 'store_true',                    help = 'Disable lazy parsing (increases load time but sometimes useful).' )

        self.path = None

    def protare(self, args, verbosity=0, lazyParsing=True):
        """
        Returns a protare (i.e., reactionSuite) instance that has been loaded per the "mapOrProtareFileName", "--pid" and "--tid" arguments.
        """

        lazyParsing = lazyParsing and not args.noLazyParse

        mapOrProtare = GNDS_fileModule.read(args.mapOrProtareFileName, verbosity=verbosity, lazyParsing=lazyParsing)

        if isinstance(mapOrProtare, reactionSuiteModule.ReactionSuite):
            self.path = args.mapOrProtareFileName
            return mapOrProtare
        elif isinstance(mapOrProtare, mapModule.Map):
            if args.tid is None: raise ValueError('--tid option required for the map file.')
            mapProtare = mapOrProtare.find( args.pid, args.tid )
            self.path = mapOrProtare.path if os.path.isabs(mapProtare.path) else os.path.join(os.path.dirname(args.mapOrProtareFileName), mapProtare.path)
            if mapProtare is None: raise Exception('No match in map file for pid = "%s" and tid = "%s".' % (args.pid, args.tid))
            return mapProtare.protare(verbosity=verbosity, lazyParsing=lazyParsing)
        else:
            raise TypeError('File "%s" is not a map or protare file.' % mapOrProtare)

def mapFromMapOrProtarePath(path):
    """
    This function takes as an argument a path that must point to a map file or a protare file or a list of protare file, and returns a :py:class:`mapModule.Map` instance.

        -) if path is a map file, it reads in the map file and returns the :py:class:`mapModule.Map` instance, or
        -) if path is a protare file, creates a map file, adds the protare to it and returns the :py:class:`mapModule.Map` instance.
        -) if path is a list of paths to protare files, creates a map file, adds the protares to it and returns the :py:class:`mapModule.Map` instance.

    This function was written to aid in writing scripts that need to either loop over a map file and perform an action on each
    protare in the map file, or read in a protare and perform the same action on it.

    :param path:            Path to map or protare file, or a list of paths to protare files.

    :return:                A :py:class:`mapModule.Map` instance.
    """

    if not isinstance(path, str):               # Support path being list with one map file (i.e., ['/path/to/map/file']).
        try:
            if len(path) == 1:
                path = [p for p in path][0]
        except:
            pass

    if isinstance(path, (str, pathlib.Path)):
        name, data = GNDS_fileModule.type(path)

        if name == reactionSuiteModule.ReactionSuite.moniker:
            map1 = mapModule.Map('Dummy', 'dummy.map')
            interaction = enumsModule.Interaction.fromString(data['interaction'])
            map1.append(mapModule.Protare(data['projectile'], data['target'], data['evaluation'], path, interaction))
            return map1
        elif name == mapModule.Map.moniker:
            return mapModule.read(path)
        else :
            raise TypeError('File "%s" is not a map or protare file.' % mapOrProtare)
    else:                               # Must be a list of protare files.
        map1 = mapModule.Map('Dummy', 'dummy.map')
        for file in path:
            name, data = GNDS_fileModule.type(file)
            if name == reactionSuiteModule.ReactionSuite.moniker:
                interaction = enumsModule.Interaction.fromString(data['interaction'])
                map1.append(mapModule.Protare(data['projectile'], data['target'], data['evaluation'], file, interaction))
            else:
                raise TypeError('Path is not a protare (i.e., reactionSuite) file: %s' % file)

        return map1
