# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the URR probability tables if they exists.
"""

import numpy
from xData import XYs1d as XYs1dModule
from fudge import styles as stylesModule
from fudge.reactionData import crossSection as crossSectionModule

inelasticMTs = [4] + list(range(51, 92))

def getILF_IOA(self, skipILF_logic=False):
    '''
    This function determines the ACE ILF and IOA parameters.

    :param skipILF_logic:   If True, the returned ILF and IOA values are both -1 This is for debugging and should be False otherwise.

    :return:                Integer ILF and IOA values.
    '''

    ILF = -1
    ILFs = []
    IOA = -1
    IOAs = []

    unresolved = self.resonances.unresolved
    if unresolved and not skipILF_logic:
        for reaction in self.reactions:
            if isinstance(reaction.crossSection[0], crossSectionModule.ResonancesWithBackground):
                continue
            if reaction.ENDF_MT in inelasticMTs:
                if reaction.effectiveThreshold() < unresolved.domainMax:
                    ILFs.append(reaction.ENDF_MT)
            else:
                if reaction.effectiveThreshold() < unresolved.domainMax:
                    IOAs.append(reaction.ENDF_MT)

    if len(ILFs) == 1:
        ILF = ILFs[0]
    elif len(ILFs) > 1:
        ILF = 4

    if len(IOAs) > 0:
        IOA = 0

    return ILF, IOA

def toACE_tables(reactionCrossSection, GNDS_URR_data, cumulatives):

    ACE_URR_data = {}
    data = GNDS_URR_data.data
    for function in data.functionals:
        pdf = XYs1dModule.XYs1d([function.xs.values, function.pdf.values], dataForm = 'xsandys')
        cdf = XYs1dModule.XYs1d([function.xs.values, function.cdf.values], dataForm = 'xsandys')
        inv_cdf = cdf.inverse()

        ACE_URR_crossSections = []
        xs = [inv_cdf.evaluate(cumulative) for cumulative in cumulatives] + [inv_cdf[-1][1]]
        for i1, x2 in enumerate(xs):
            if i1 > 0:
                ACE_URR_crossSections.append(pdf.integrateWithWeight_x(domainMin = x1, domainMax = x2) / pdf.integrate(domainMin = x1, domainMax = x2))
            x1 = x2
        ACE_URR_data[function.outerDomainValue] = [reactionCrossSection.evaluate(function.outerDomainValue), ACE_URR_crossSections]

    return ACE_URR_data

def URR_probabilityTable(self, styleLabel, fromPDF=False, numberOfProbabilities=20, skipILF_logic=False, verbose=0):
    """
    Dump probability tables to ACE. By default use the probability table
    data from <applicationData>, or if fromPDF=True compute the probability
    tables from cross-section pdfs (not expected to work as well).
    numberOfProbabilities is ignored unless fromPDF=True
    """

    URRPT = []
    IFF = 1
    ILF, IOA = getILF_IOA(self, skipILF_logic=skipILF_logic)
    URR_styles = [style for style in self.styles.getStylesOfClass(stylesModule.URR_probabilityTables)
                  if style.findDerivedFromStyle(stylesModule.GriddedCrossSection).label == styleLabel]
    if len(URR_styles) != 1:
        if verbose > 0:
            print("Skipping URR as no unique URR style corresponding to label %s found" % styleLabel)
        return ILF, IOA, URRPT

    URR_style = URR_styles[0]

    if fromPDF:
        heatedStyle = URR_style
        while True:
            heatedStyle = self.styles[heatedStyle.derivedFrom]
            if heatedStyle is None:
                raise Exception('Could not find heated style for URR data.')
            if isinstance(heatedStyle, stylesModule.Heated):
                break

        cumulatives = [i1 / float(numberOfProbabilities) for i1 in range(numberOfProbabilities)]
        URR_label = URR_style.label
        MTs_URR = {}
        for reaction in self.reactions:
            if URR_label in reaction.crossSection:
                reactionCrossSection = reaction.crossSection[heatedStyle.label]
                ACE_I_index = { 2: 3, 18: 4, 102: 5 }[reaction.ENDF_MT]
                MTs_URR[ACE_I_index] = toACE_tables(reactionCrossSection, reaction.crossSection[URR_label], cumulatives)

        energies = sorted(MTs_URR[3].keys())
        URRPT = [len(energies), numberOfProbabilities, 2, ILF, IOA, IFF] + energies

        ones = numberOfProbabilities * [1.0]
        cumulatives = cumulatives + [1.0]
        probabilities = cumulatives[1:]
        for energy in energies:
            elastic = MTs_URR[3][energy]
            fission = [0.0, ones]
            if 4 in MTs_URR:
                fission = MTs_URR[4][energy]
            capture = MTs_URR[5][energy]
            total = []
            totalCrossSection = elastic[0] + fission[0] + capture[0]
            for i1 in range(len(elastic[1])):
                total.append((elastic[0] * elastic[1][i1] + fission[0] * fission[1][i1] + capture[0] * capture[1][i1]) / totalCrossSection)

            URRPT += probabilities
            URRPT += total
            URRPT += elastic[1]
            URRPT += fission[1]
            URRPT += capture[1]
            URRPT += ones

    else: # use FUDGE probability tables
        from fudge.resonances import probabilityTables
        pts = self.applicationData[probabilityTables.LLNLProbabilityTablesToken]
        ptNow = pts.probabilityTables[URR_style.label]

        energies = [pt.value for pt in ptNow]
        numberOfProbabilities = len(ptNow[0].table)
        URRPT = [len(energies), numberOfProbabilities, 2, ILF, IOA, IFF] + energies

        columnLabels = [
            self.getReaction('total').label,
            self.getReaction('elastic').label,
            None,   # placeholder for fission
            self.getReaction('capture').label,
            None    # enter dummy values for last column (heating)
        ]
        if self.getReaction('fission') is not None:
            fissionLabel = self.getReaction('fission').label
            if fissionLabel in [k.name for k in ptNow[0].table.columns]:
                columnLabels[2] = fissionLabel

        for pt in ptNow:  # loop over incident energies
            probabilities = pt.table.getColumn('probability')
            URRPT += list(numpy.cumsum(probabilities))
            for columnLabel in columnLabels:
                if columnLabel is None:
                    URRPT += [1] * len(pt.table)
                else:
                    URRPT += pt.table.getColumn(columnLabel)

    return ILF, IOA, URRPT
