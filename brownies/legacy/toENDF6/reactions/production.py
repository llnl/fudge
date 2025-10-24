# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs import IDs as IDsPoPsModule

from fudge.reactions import production as productionModule

#
# production
#
def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):


    if flags['verbosity'] >= 10:
        print('%sproduction: %s' % (verbosityIndent, self.outputChannel.toString(simpleString = True)))
    MT = self.ENDF_MT
    conversionFlag = targetInfo['ENDFconversionFlags'].get(self, None)
    if conversionFlag == 'MF3production':
        # MT203-207 production data stored in MF=3 instead of 8 & 10. Treat like normal reaction:
        from . import base
        base.toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=verbosityIndent)
        return

    MF = None
    for product in self.outputChannel:
        if product.pid == IDsPoPsModule.photon:
            if product.distribution.hasData():
                productMF = 12
                if 'MF13' in targetInfo['ENDFconversionFlags'].get(product,""): productMF = 13
                if MF == None: MF = productMF
                if MF != productMF: raise Exception('MF = %s mixed with MF = %s' % (MF, productMF))
    if MF is not None:
        outputChannel = self.outputChannel
        productName = outputChannel[0].pid
        if productName != IDsPoPsModule.photon:
            raise Exception('production only supported for gamma as product: product is "%s"' % productName)

        if MT not in targetInfo['production_gammas']: targetInfo['production_gammas'][MT] = [MF]
        if MF != targetInfo['production_gammas'][MT][0]:
            raise Exception('Prior instances were MF = %d, which cannot be mixed with MF = %d' %
                            (targetInfo['production_gammas'][MT][0], MF))
        targetInfo['production_gammas'][MT].append(self)
    else:
        if MT not in targetInfo['MF8']: targetInfo['MF8'][MT] = []
        targetInfo['MF8'][MT].append(self)

productionModule.Production.toENDF6 = toENDF6
