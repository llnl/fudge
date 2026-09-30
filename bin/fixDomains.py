#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """
"""

import sys
import os
import argparse

from fudge import reactionSuite as reactionSuiteModule
from fudge.reactions import base as reactionBaseModule

parser = argparse.ArgumentParser(description=description)
parser.add_argument("inputFile", type=str,                                      help="""The input file whose 'sums' will be fix.""")
parser.add_argument("outputFile", type=str,                                     help="""The output file name of fixed file.""")
parser.add_argument("--energyMax", type=float, default=-1.0,                    help="""Trims any incident energy domains to energyMax.""")
parser.add_argument("--onlyWriteWhenModified", action='store_false',            help="""If this option is entered and there are no domain fixes, then the output file is not written.""")

args = parser.parse_args()

protare = reactionSuiteModule.read(args.inputFile, verbosity=0)

def checkSums(reaction):

    def isReaction(link):

        linkReaction = link.link.findClassInAncestry(reactionBaseModule.Base_reaction)

        return reaction is linkReaction

    crossSectionSumsToRemove = []
    for crossSectionSum in protare.sums.crossSectionSums:
        summandsToRemove = []
        for index, summand in enumerate(crossSectionSum.summands):
            if isReaction(summand):
                summandsToRemove.append(index)
        for index in reversed(summandsToRemove):
            crossSectionSum.summands.summands.pop(index)
        if len(crossSectionSum.summands) == 0:
            crossSectionSumsToRemove.append(crossSectionSum.label)
    for label in crossSectionSumsToRemove:
        print('    INFO: removing empty crossSectionSum "%s".' % label)
        protare.sums.crossSectionSums.pop(label)

    multiplicitySumToRemove = []
    for multiplicitySum in protare.sums.multiplicitySums:
        if len(multiplicitySum.summands) == 0:
            multiplicitySumToRemove.append(multiplicitySum.label)
        else:
            if isReaction(multiplicitySum.summands[0]):
                multiplicitySumToRemove.append(multiplicitySum.label)
    for label in multiplicitySumToRemove:
        print('    INFO: removing empty multiplicitySum "%s".' % label)
        protare.sums.multiplicitySums.pop(label)

def removeReactionsWithHighThresholds(reactions):

    reactionsToRemove = [reaction for reaction in reactions if reaction.crossSection.effectiveThreshold() >= args.energyMax]
    for reaction in reactionsToRemove:
        checkSums(reaction)
        print('    INFO: removing high threshold reaction "%s".' % reaction.label)
        reactions.remove(reaction.label)

    return len(reactionsToRemove)

if args.energyMax == -1.0:
    args.energyMax = protare.styles.projectileEnergyDomain().max

numberOfFixes = 0
numberOfFixes += removeReactionsWithHighThresholds(protare.reactions)
numberOfFixes += removeReactionsWithHighThresholds(protare.orphanProducts)
numberOfFixes += removeReactionsWithHighThresholds(protare.productions)
numberOfFixes += removeReactionsWithHighThresholds(protare.fissionComponents)
numberOfFixes += removeReactionsWithHighThresholds(protare.incompleteReactions)

numberOfFixes += protare.fixDomains(args.energyMax)

if numberOfFixes != 0 or args.onlyWriteWhenModified:
    protare.saveToFile(args.outputFile, formatVersion=protare.format)
