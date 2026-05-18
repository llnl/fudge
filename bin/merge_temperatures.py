#!/usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
from fudge import GNDS_file
from fudge import reactionSuite
from fudge import styles

summaryDocstring__FUDGE = "Merge processed GNDS files at various temperatures into one multi-temperature file"

description = """Script for merging processed GNDS files at different temperatures.
Running processProtare.py with many temperatures can be expensive. It may be useful to process fewer temperatures at a time,
then merge them using this script.

Sample use:
    processProtare.py <evaluation.xml> processed_1.xml -t 0 -t 293.6 -t 600 --temperatureUnit K --baseTemperatureIndex 0 ...
    processProtare.py <evaluation.xml> processed_2.xml -t 900 -t 1200 -t 2000 --temperatureUnit K --baseTemperatureIndex 3 ...

    merge_temperatures.py processed_1.xml processed_2.xml -o final_file.xml --hybrid
"""

def parse_args():
    parser = argparse.ArgumentParser(description)
    parser.add_argument("gndsFiles", nargs="+", help="Files to merge")
    parser.add_argument("-o", "--outputFile", required=True, help="Final merged file name")
    parser.add_argument("--hybrid", action="store_true", help="Write final file in hybrid XML/HDF5")

    return parser.parse_args()


def copy_component(src_comp, dest_comp, toCopy):
    """ Copy all processed forms in a component like <crossSection>, <multiplicity>, etc. """
    for key in src_comp.keys():
        if key in toCopy:
            dest_comp.add(src_comp[key])


def copy_styles(src_styles, dest_styles, toCopy):
    """ Styles section requires special handling """
    for key in src_styles.keys():
        if key not in toCopy: continue

        style = src_styles[key]
        # don't repeat processing documentation for higher temps:
        style.documentation.dates.clear()
        style.documentation.computerCodes.clear()  

        if isinstance(style, styles.HeatedMultiGroup):
            # replace transportables with a link to the first temperature
            style.transportables.clear()
            dest_lowest_temp = dest_styles.temperatures()[0]
            # FIXME href doesn't support assignment
            style.transportables._Suite__href = dest_styles[dest_lowest_temp.heatedMultiGroup].transportables.toXLink()

        dest_styles.add(style)


def copy_product(src_product, dst_product, toCopy):
    copy_component(src_product.multiplicity, dst_product.multiplicity, toCopy)
    copy_component(src_product.distribution, dst_product.distribution, toCopy)
    copy_component(src_product.averageProductEnergy, dst_product.averageProductEnergy, toCopy)
    copy_component(src_product.averageProductMomentum, dst_product.averageProductMomentum, toCopy)

    if src_product.outputChannel:
        copy_outputChannel(src_product.outputChannel, dst_product.outputChannel, toCopy)


def copy_outputChannel(src_channel, dst_channel, toCopy, fission=False):
    copy_component(src_channel.Q, dst_channel.Q, toCopy)

    if fission:
        for other_dn, base_dn in zip(
                src_channel.fissionFragmentData.delayedNeutrons,
                dst_channel.fissionFragmentData.delayedNeutrons):
            copy_product(other_dn.product, base_dn.product, toCopy)

        copy_component(
            src_channel.fissionFragmentData.fissionEnergyReleases,
            dst_channel.fissionFragmentData.fissionEnergyReleases, toCopy)

    for src_product in src_channel.products:
        dst_product = dst_channel.products[src_product.label]
        copy_product(src_product, dst_product, toCopy)


def copy_reaction(src_reaction, dst_reaction, toCopy):
    """
    Copy all styles from the toCopy list from src_reaction to dst_reaction
    """

    copy_component(src_reaction.crossSection, dst_reaction.crossSection, toCopy)
    copy_component(src_reaction.availableEnergy, dst_reaction.availableEnergy, toCopy)
    copy_component(src_reaction.availableMomentum, dst_reaction.availableMomentum, toCopy)

    copy_outputChannel(src_reaction.outputChannel, dst_reaction.outputChannel, toCopy,
            fission=src_reaction.isFission())


# --------------------------------------------------------
# Merge multiple ReactionSuites (first file is accumulator)
# --------------------------------------------------------
def merge_reaction_suites(files):
    print(f"Loading base ReactionSuite from {files[0]}")
    base = reactionSuite.read(files[0])

    for idx, fn in enumerate(files[1:]):
        other = reactionSuite.read(fn)

        print(f"Merging ReactionSuite from {fn}")
        toCopy = set([getattr(temp, processedStyle)
            for processedStyle in ('heated', 'griddedCrossSection', 'heatedMultiGroup', 'URR_probabilityTables', 'heatedMultiGroup', 'SnElasticUpScatter')
            for temp in other.styles.temperatures()
            ])

        for otherStyle in toCopy:
            if otherStyle in base.styles:
                # FIXME should this script have an option to re-index style labels?
                raise Exception(f"Style {otherStyle} is already present in merged file! May indicate missing or incorrect '--baseTemperatureIndex' when running processProtare!")

        copy_styles(other.styles, base.styles, toCopy)

        for other_rx in other.reactions:
            base_rx = base.reactions[other_rx.label]
            copy_reaction(other_rx, base_rx, toCopy)

        for other_op in other.orphanProducts:
            base_op = base.orphanProducts[other_op.label]
            copy_component(other_op.crossSection, base_op.crossSection, toCopy)
            copy_outputChannel(other_op.outputChannel, base_op.outputChannel, toCopy)

        for other_crossSectionSum in other.sums.crossSectionSums:
            base_crossSectionSum = base.sums.crossSectionSums[other_crossSectionSum.label]
            copy_component(other_crossSectionSum.crossSection, base_crossSectionSum.crossSection, toCopy)

        for other_production in other.productions:
            base_production = base.productions[other_production.label]
            copy_component(other_production.crossSection, base_production.crossSection, toCopy)
            copy_component(other_production.outputChannel.Q, base_production.outputChannel.Q, toCopy)

        # LLNL-specific stuff in applicationData:
        for key, copy_method in (
                ("LLNL::multiGroupReactions", copy_reaction),
                ("LLNL::multiGroupDelayedNeutrons", copy_outputChannel),
                ("LLNL::URR_probability_tables", copy_component)
                ):
            if key in base.applicationData:
                for other_appData, base_appData in zip(
                        other.applicationData[key],
                        base.applicationData[key]):
                    copy_method(other_appData, base_appData, toCopy)

    return base


if __name__ == "__main__":
    args = parse_args()

    previews = [GNDS_file.preview(gndsFile) for gndsFile in args.gndsFiles]
    # sanity check:
    for field in ('projectile', 'target', 'evaluation'):
        field_values = set([getattr(f, field) for f in previews])
        if len(field_values) != 1:
            raise Exception(f"Error! {field} must match for all input files, but found multiple values: {field_values}")

    first_temperatures = [f.styles.temperatures()[0].temperature for f in previews]
    # sort by temperature
    first_temperatures, previews, gndsFiles = zip(*sorted(zip(first_temperatures, previews, args.gndsFiles)))

    # temperatures must be increasing monotonically:
    all_temps = [t.temperature for preview in previews for t in preview.styles.temperatures()]
    if all_temps != sorted(all_temps):
        raise Exception(f"Error! Temperatures were out of order! {all_temps}")

    print(f"Merging {len(all_temps)} temperatures from {len(gndsFiles)} files")
    merged = merge_reaction_suites(gndsFiles)


    print(f"\nWriting merged ReactionSuite to {args.outputFile}")
    merged.saveAllToFile(args.outputFile, hybrid=args.hybrid)

