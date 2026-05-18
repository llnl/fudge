#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from PoPs import IDs as PoPsIDsModule
from PoPs.families import particle as PoPsParticleModule
from fudge import reactionSuite as reactionSuiteModule

description = '''For a protare (i.e., reactionSuite), list each reaction and its outgoing particles (i.e., products).'''

class Formats:
    """
    This classes stores for the RIS formats and the member **allowed** with are all the allowed formats.  Currently, 
    there is only one allowed format which is '1.0'. However, this format can be one of two formats: old and new. 

    The old 1.0 format has only 4 columns for the items under the '#reactions' node. These items are
    label, threshold, intermediate and process. The new version of the 1.0 format has two additional column which are
    reactionLabel and covarianceFlag. While prior versions of FUDGE cannot read the new format, newer versions of FUDGE
    can read either format. Also, prior verions of GIDI+ can read either format but only translate the first 4 columns.
    New versions of GIDI++ read the additional 2 columns.
    """

    v1_0 = '1.0'

    allowed = (v1_0, )

__doc__ = description

def protareProductInfo(fileName, indent, format=Formats.v1_0):
    """
    Opens the protare file *fileName* and calls **protareProductInfo2**.

    :param fileName:    The file name of the protare to read in.
    :param indent:      The amount of indentation for each line.

    :return:            The returned results from **protareProductInfo2**.
    """

    return protareProductInfo2(reactionSuiteModule.read(fileName, lazyParsing=True), indent, format=format)

def protareProductInfo2(protare, indent, format=Formats.v1_0):
    """
    Constructs a line for a **ris** file for each reaction in *protare*. The list of lines constructed is returned.
    Each line contains four items separated by the colon (i.e., ":") character. The four items are:

        -) A formula for all final products.
        -) The effective threshold for the reaction.
        -) An intermediate product if one exists.
        -) The reaction process.

    :param protare:     A ReactionSuite instance whose reactions are analyzed.
    :param indent:      The amount of indentation for each line.

    :return:            The list of lines constructed.
    """

    
    width = max(map(len, list(reaction.label for reaction in protare.reactions)))
    fmtLabel = '%%-%ds' % width

    width = 9
    for reaction in protare.reactions:
        if reaction.outputChannel.process is not None:
            width = max(width, len(str(reaction.outputChannel.process)))
    fmtProcess = '%%-%ds' % width
        
    lines = []
    for reaction in protare.reactions:
        effectiveThreshold = reaction.crossSection.effectiveThreshold()

        reactionProducts = reaction.reactionProducts()
        particles2, intermediates, processes, dummy = reactionProducts.asSortedList(protare.PoPs, excludePhotons=False)
        particles = {}
        for pid, multiplicity in particles2:
            particles[pid] = multiplicity

        sortedParticles = PoPsParticleModule.ParticleSorter()
        for pid in particles:
            sortedParticles.add(protare.PoPs[pid])

        photons = ''
        label = ''
        sep = ''
        for particle in sortedParticles:
            multiplicity = ''
            if particles[particle.id] != 1:
                multiplicity = '%s' % particles[particle.id]
            if particle.id == PoPsIDsModule.photon:
                photons = '%s%s' % (multiplicity, particle.id)
            else:
                label += '%s%s%s' % (sep, multiplicity, particle.id)
                sep = ' + '
        if len(photons) > 0:
            label += '%s%s' % (sep, photons)

        intermediate = ''
        if len(intermediates) > 0:
            intermediate = intermediates[0][0]
        process = ''
        if len(processes) > 0:
            process = processes[0][0]

        hasCovariance = ''
        if len(reaction.crossSection) > 0:
            crossSection = reaction.crossSection[0]
            if hasattr(crossSection, 'uncertainty'):
                if crossSection.uncertainty is not None:
                    if crossSection.uncertainty.data is not None:
                        hasCovariance = 'covariance'

# All outputted ':' must be surrounded by spaces (i.e., ' : ') for parsing as some RIS files contain the process "ENDL:S2". This requirement was added 20/Aug/2025.
        lines.append('%s%-31s : %-12s : %-10s : %s : %s : %s' % (indent, label, effectiveThreshold, intermediate, 
                fmtProcess % process, fmtLabel % reaction.label, hasCovariance))

    return lines

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('input', type=str,                                          help='''The protare (i.e., reactionSuite) whose reactions are analyzed.''')
    parser.add_argument('--format', choices=Formats.allowed, default=Formats.v1_0,  help='''The format version of the printed RIS lines.''')

    args = parser.parse_args()
    print('\n'.join(protareProductInfo(args.input, '    ', args.format)))
