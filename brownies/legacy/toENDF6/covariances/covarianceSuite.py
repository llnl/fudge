# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.families import nuclide as nuclideModule
from xData import enums as xDataEnumsModule

from brownies.legacy.converting import endf_endl as endf_endlModule

from fudge.covariances import covarianceSuite as covarianceSuiteModule
from fudge.covariances import mixed as covarianceMixedModule

from .. import endfFormats as endfFormatsModule
from .. import gndsToENDF6 as gndsToENDF6Module
from .modelParameters import averageParametersToENDF6

def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
    """Convert to ENDF format."""

    reactionSuite = targetInfo['reactionSuite']

    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN = 1, 1.0
    ZAI = ZAM  # assuming one isotope/file
    MTL = 0 # mtl=1 sections are handled in lumpedCovariance

    if self.parameterCovariances:
        self.parameterCovariances.toENDF6(endfMFList, flags, targetInfo, verbosityIndent)

    sections = self.covarianceSections
    # sort covariances by MF/MT:
    mfmts = []
    for section_ in sections:
        mfmts.append( list( map( int, section_.rowData.ENDF_MFMT.split(',') ) ) )
    if len(mfmts)==0: return

    mfs, mts = list( zip( *mfmts ) )
    zipList = list( zip( mfs, mts, sections ) )
    idx = 0
    while idx<len(zipList):
        mf,mt,covar = zipList[idx]
        thisMFMT = [a[2] for a in zipList if a[:2]==(mf,mt)]
        idx += len(thisMFMT)
        endf = []

        if mf in (31,33):
            endf += [endfFormatsModule.endfHeadLine(ZAM, AWT, 0, MTL, 0, len(thisMFMT))]
        elif mf==34:
            LTT = 1
            if any([slice.domainValue == 0 for section in thisMFMT for slice in section.rowData.slices]):
                LTT = 3
            NMT1 = 1    # cross-reaction terms not yet supported
            endf += [(endfFormatsModule.endfHeadLine( ZAM, AWT, 0, LTT, 0, NMT1))]
            MAT1 = 0
            MT1 = mt
            L1s, L2s = [],[]
            for section_ in thisMFMT:
                L1s.append( section_.rowData.slices[0].domainValue )
                if section_.columnData:
                    L2s.append( section_.columnData.slices[0].domainValue )
                else:
                    L2s.append( L1s[-1] )
            NL = len(set(L1s))
            NL1 = len(set(L2s))
            endf += [endfFormatsModule.endfHeadLine(0.0, 0.0, MAT1, MT1, NL, NL1)]
        elif mf==35:
            endf += [endfFormatsModule.endfHeadLine(ZAM, AWT, 0, MTL, len(thisMFMT), 0)]
        elif mf==40:
            endf += [endfFormatsModule.endfHeadLine(ZAM, AWT, 0, 0, len(thisMFMT), 0)]
        for section_ in thisMFMT:
            MAT1 = 0
            form = gndsToENDF6Module.getForm( targetInfo['style'], section_ )
            conversionFlags = targetInfo['ENDFconversionFlags'].get(form, "")
            if (section_.columnData is not None and
                    section_.columnData.root is not None and
                    section_.columnData.root[1:] in self.externalFiles):
                otherTarget = section_.columnData.root[1:]
                if otherTarget != 'reactions':
                    ZA, MAT1 = endf_endlModule.ZAAndMATFromParticleName( otherTarget )
            if mf==34:
                L1 = int(section_.rowData.slices[0].domainValue)
                if section_.columnData:
                    L2 = int(section_.columnData.slices[0].domainValue)
                else:
                    L2 = L1
                NI = 1
                if isinstance(form, covarianceMixedModule.MixedForm):
                    NI = len(form)
                    frame = form[0].productFrame
                else:
                    frame = form.productFrame
                LCT = {xDataEnumsModule.Frame.lab: 1, xDataEnumsModule.Frame.centerOfMass: 2}[frame]
                if 'LCT=0' in conversionFlags:
                    LCT = 0
                endf += [ endfFormatsModule.endfHeadLine( 0.0, 0.0, L1, L2, LCT, NI ) ]
            if mf==40:
                rowData = section_.rowData
                if isinstance(rowData, str):
                    raise Exception( "Don't string me along!" ) # FIXME
                else:
                    quant = rowData.link
                    if mt == 18:
                        QM, QI = quant.findAttributeInAncestry('getQ')('eV'), 0
                        LFS, level, IZAP = 0,0,-1
                    else:
                        product = quant.findAttributeInAncestry('outputChannel')[0]
                        QI = quant.findAttributeInAncestry('getQ')('eV')
                        LFS, level = 0, 0.
                        particle = reactionSuite.PoPs[product.pid]
                        if( isinstance( particle, nuclideModule.Particle ) ) :
                            LFS = particle.index
                            level = particle.energy[0].float( 'eV' )
                        QM = QI + level
                        IZAP = chemicalElementMiscPoPsModule.ZA( reactionSuite.PoPs[product.pid] )
                    NL = 1
                    endf += [endfFormatsModule.endfHeadLine( QM, QI, IZAP, LFS, 0, NL ) ]
                    XMF1, XLFS1, NC, NI = 10,LFS, 0,1
                    endf += [endfFormatsModule.endfHeadLine( XMF1,XLFS1,MAT1,mt,NC,NI )]
            targetInfo['MAT1'] = MAT1
            targetInfo['dataPointer'] = [section_.rowData,section_.columnData]
            endf += form.toENDF6(flags, targetInfo)
            targetInfo.pop( 'dataPointer' )
            targetInfo.pop( 'MAT1' )
        endf.append( endfFormatsModule.endfSENDLineNumber() )
        if mt not in endfMFList[mf]:
            endfMFList[mf][mt] = []
        endfMFList[mf][mt] += endf
    # also add ENDF-style pointers for lumped covariance data:
    for reactionSum in targetInfo['reactionSuite'].sums.crossSectionSums:
        MT1 = reactionSum.ENDF_MT
        if MT1 not in list( range( 851, 872 ) ) : continue
        for summand in reactionSum.summands:
            mt = summand.link.findAttributeInAncestry('ENDF_MT')
            endfMFList[33][mt] = [endfFormatsModule.endfHeadLine(ZAM,AWT,0,MT1,0,0),
                    endfFormatsModule.endfSENDLineNumber() ]
    return

covarianceSuiteModule.CovarianceSuite.toENDF6 = toENDF6
