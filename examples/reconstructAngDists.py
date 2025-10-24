# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import GNDS_file
from fudge import reactionSuite as reactionSuiteModule
from fudge.productData.distributions import angular as angularModule
from fudge.processing.resonances import reconstructResonances


def parseArgs():
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser(
        description='Reconstruct angular distributions from a GNDS file, optionally replacing original distribution')
    parser.add_argument('gnds', type=str, help='Initial GNDS file')
    parser.add_argument('--mt', type=int, help='MT of the reaction to replace / write to table.')
    parser.add_argument('--check', action='store_true', help="Check the resulting evaluation")
    parser.add_argument('--format', choices=['gnds', 'endf', 'table'], default='gnds',
                        help="Format to save the results in outFile, with results from --strategy switch")
    parser.add_argument('--strategy', choices=['merge', 'replace', 'dryrun'], default='dryrun',
                        help="How to merge with existing fast angular distribution")
    parser.add_argument('-o', dest='outFile', default=None, type=str, help='Output file')
    return parser.parse_args()


def smallBanner(x, wingsize=10):
    """Small banner, with wings"""
    return wingsize * '*' + ' ' + x.replace('\n', '; ') + ' ' + wingsize * '*'


def getPlotTable(xys2d):
    """Make a table that we can feed to a spreadsheet from an instance of distributions.angular.XYs2d"""
    if not isinstance(xys2d, angularModule.XYs2d): raise TypeError(
        "Must be of type fudge.productData.distributions.angular.XYs2d")
    table = []
    maxLegendre = max([len(legendre) for legendre in xys2d])
    for legendre in xys2d:
        table.append([0.0] * (2 + maxLegendre))
        table[-1][0] = legendre.outerDomainValue
        for i, y in enumerate(legendre): table[-1][i + 1] = y
    return ['\t'.join(map(str, row)) for row in table]


if __name__ == "__main__":

    args = parseArgs()

    print(smallBanner("Reading %s" % args.gnds))
    RS = GNDS_file.read(args.gnds)
    assert isinstance(RS, reactionSuiteModule.ReactionSuite), "Input must be a GNDS reactionSuite file!"

    if args.format == 'table' and args.mt is None:
        raise ValueError("Need to supply an MT number to write output as a table")

    # Reconstruct the angular distributions
    print(smallBanner("Reconstructing resonance angular distributions"))
    newAngDists = reconstructResonances.reconstructAngularDistributions(RS)
    # dictionary with reaction labels as keys and angular.XYs2d (containing Legendre series) as values

    if args.strategy == 'dryrun': exit()

    # If we reconstructed the angular distributions too, put them in the appropriate 
    # reaction/outputChannel/products/product/distribution/angular/XYs2d
    for key in newAngDists:
        newDist = newAngDists[key]
        reaction = RS.getReaction(key)
        if args.mt is not None and reaction.ENDF_MT != args.mt:
            continue
        # FIXME better would be for reconstructAngularDistributions to return ejectile pid along with reaction and xys2d
        ejectile = key.split(' + ')[0]
        product = reaction.outputChannel.getProductWithName(ejectile)

        if args.strategy == 'replace':
            print(smallBanner("Replacing angular distributions for reaction `%s`" % key))
        elif args.strategy == 'merge':
            print(smallBanner("Merging angular distributions for reaction `%s`" % key))
            originalDist = product.distribution.evaluated.angularSubform
            Emax = newDist.domainMax

            if not isinstance(originalDist, (angularModule.XYs2d, angularModule.Regions2d)):
                raise NotImplementedError(
                    "Don't know how to merge reconstructed data into angular distribution of type %s" %
                    type(originalDist))

            newRegions = angularModule.Regions2d(axes=originalDist.axes)
            newRegions.append(newDist)

            if isinstance(originalDist, angularModule.XYs2d):
                originalDist = [originalDist]  # pretend it was Regions2d

            for region in originalDist:
                if region.domainMax <= newDist.domainMax:
                    continue
                if region.domainMin < newDist.domainMax:
                    regionSlice = region.domainSlice(newDist.domainMax, region.domainMax)
                    newRegions.append(regionSlice)
                else:
                    newRegions.append(region)
            newDist = newRegions

        product.distribution.evaluated.subforms.pop()
        product.distribution.evaluated.subforms.append(newDist)

    # Check evaluation
    if args.check:
        print(smallBanner("Checking final evaluation"))
        print(RS.check())

    # Save the results
    if args.outFile is not None:
        print(smallBanner("Saving results as a " + args.format + " in file %s" % args.outFile))
        if args.format == 'table':
            distribution = RS.getReaction(args.mt).outputChannel[0].distribution.evaluated.angularSubform
            with open(args.outFile, mode='w') as fout:
                fout.write('\n'.join(getPlotTable(distribution)))
        elif args.format == 'endf':
            from brownies.legacy.toENDF6 import toENDF6
            CS = None
            covariances = RS.loadCovariances()
            if len(covariances):
                CS = covariances[0]
            with open(args.outFile, mode='w') as fout:
                fout.write(RS.toENDF6(styleLabel=RS.styles.getEvaluatedStyle().label,
                                      flags={'verbosity': 0}, covarianceSuite=CS))
        else:
            # save to GNDS
            RS.saveToFile(args.outFile)
