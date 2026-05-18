#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import numpy
import argparse

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import styles as stylesModule
from fudge import institution as institutionModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData import URR_probabilityTables as URR_pdfsModule
from fudge.resonances import probabilityTables as probabilityTablesModule

from fudge.processing.resonances import makeUnresolvedProbabilityTables

from xData import table as tableModule

from LUPY import times as timesModule

summaryDocString__FUDGE = """This script process unresolved resonances to create probability tables."""

parser = argparse.ArgumentParser(
    "Process unresolved resonances to create probability tables and optional cross section PDFs")
parser.add_argument("gnds", help="Input GNDS file to process. Should have already run through processProtare")
parser.add_argument("nSamples", type=int, help="Number of realizations")
parser.add_argument("-o", "--output", help="Output GNDS file to write")
parser.add_argument("-g", "--generator", default="wigner", help="level generator: 'wigner', 'goe', etc.")
parser.add_argument("-f", "--formatVersion", default=GNDS_formatVersionModule.default,
                    help="GNDS format version to write out")
parser.add_argument("--skipPDFs", action="store_true", help="Generate probability tables but not PDFs. A bit faster.")
parser.add_argument("--hybrid", action="store_true", help="Write hybrid XML/HDF5 files.")
parser.add_argument("--debug", help="Save realization details to specified file (pickle format).")
parser.add_argument("-v", "--verbose", action="store_true", help="enable verbose output")
args = parser.parse_args()

if args.output is None:
    args.output = args.gnds + "_urr.xml"

RS = reactionSuiteModule.read(args.gnds)

if not RS.resonances or not RS.resonances.unresolved:
    if args.verbose:
        print(f"Skipping URR processing since {args.gnds} contains no URR data")
    RS.saveAllToFile(args.output, xs_pdf_cdf1d_singleLine=True, formatVersion=args.formatVersion, hybrid=args.hybrid)
    sys.exit()

energyUnit = RS.styles.getEvaluatedStyle().projectileEnergyDomain.unit
urrStyles = RS.styles.getStylesOfClass(stylesModule.URR_probabilityTables)
assert len(urrStyles) == 0, "Error: file already contains processed URR data"

griddedStyles = RS.styles.getStylesOfClass(stylesModule.GriddedCrossSection)
assert len(griddedStyles) > 0, ("No processed Monte Carlo data found! "
                                "Make sure to include the '-mc' option when running processProtare.py")

temperatures = [style.temperature.value for style in griddedStyles]
temperatureUnits = [style.temperature.unit for style in griddedStyles]
assert len(set(temperatureUnits)) == 1, "Inconsistent temperature units"

print(f"Processing URR tables at {len(temperatures)} temperatures")
tableGenerator = makeUnresolvedProbabilityTables.ProbabilityTableGenerator(RS)

timer = timesModule.Times()
results = tableGenerator.generatePDFs(args.nSamples, temperatures, temperatureUnit=str(temperatureUnits[0]),
                                      verbose=args.verbose, style=args.generator, makePDFs=(not args.skipPDFs),
                                      debugFile=args.debug)
print(timer.toString())

# package results into GNDS:
probabilityTables = probabilityTablesModule.ProbabilityTables()
columnLabels = ['probability', 'total', 'elastic', 'capture', 'fission']
columnHeaders = []
for idx, label in enumerate(columnLabels):
    reaction = RS.getReaction(label)
    if reaction is not None:
        label = reaction.label
    columnHeaders.append(tableModule.ColumnHeader(idx, label, ''))  # FIXME assumes tables are unitless

for idx, gridded in enumerate(griddedStyles):
    urrStyle = stylesModule.URR_probabilityTables(label="URR_%s" % str(idx).zfill(3), derivedFrom=gridded.label)
    RS.styles.add(urrStyle)

    PT = probabilityTablesModule.ProbabilityTable(label=urrStyle.label)
    for energy, probabilities, tableData in results['probabilityTables'][gridded.temperature.value]:
        headers = columnHeaders[:]
        if 'fission' not in tableData:
            headers.pop()

        data = [probabilities]
        for key in columnLabels:
            if key in tableData:
                data.append(tableData[key])
        data = numpy.array(data).T.tolist()

        table_ = tableModule.Table(columns=headers, data=data)

        PT.add(probabilityTablesModule.IncidentEnergy(energy, energyUnit, table_))

    probabilityTables.add(PT)

    if results['pdfs'] is not None:
        # also store cross section pdfs for each reaction
        pdfsNow = results['pdfs'][gridded.temperature.value]
        for key in pdfsNow:
            # FIXME probably need to do some smoothing on each xys1d to avoid ginormous files
            # or maybe use kernel density estimator?
            xys2d = URR_pdfsModule.XYs2d(axes=pdfsNow[key].axes)
            for xys1d in pdfsNow[key]:
                integral = xys1d.integrateWithWeight_x()
                xys1d.scaleOffsetXAndY(xScale=1/integral, yScale=integral, insitu=True)
                xys2d.append(URR_pdfsModule.Xs_pdf_cdf1d.fromXYs(xys1d))
            URR_pdfs = crossSectionModule.URR_probabilityTables1d(urrStyle.label, xys2d)
            reaction = RS.getReaction(key)
            reaction.crossSection.add(URR_pdfs)

institution = institutionModule.Institution(probabilityTablesModule.LLNLProbabilityTablesToken)
institution.append(probabilityTables)
RS.applicationData.add(institution)

RS.saveAllToFile(args.output, xs_pdf_cdf1d_singleLine=True, formatVersion=args.formatVersion, hybrid=args.hybrid)
