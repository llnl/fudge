# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module is for FUDGE internal use.
The functions in this module are used to calculate product matrices for neutron elastic scattering including a non-zero
temperature of the target.

This module contains the following functions:
            
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Function                              | Description                                                                       |
    +=======================================+===================================================================================+
    | SnElasticUpScatter                    |                                                                                   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | executeCommand                        | This function is used internally to call                                          |
    |                                       | :py:func:`subprocessingModule.executeCommand`                                     |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | checkNegative_l0                      | This function counts the number of negative elements in the l=0 product matrix.   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | parseAveEnergyOutputFile              | This function parses the average product data in the file *fileName*.             |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | rescaleCrossSection                   |                                                                                   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
"""

import os, sys, subprocess

import fudge as fudgeModule
from LUPY import subprocessing as subprocessingModule
from LUPY import times as timesModule
from fudge.processing.deterministic import transferMatrices as transferMatricesModule

def SnElasticUpScatter(style, tempInfo):
    """
    This function generates input for and calls the code which includes a target temperature when
    calculating product matrices for neutron elastic scattering.

    :param style:       This is not used.
    :param tempInfo:    This is a dictionary with needed data.
    """

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']
    restart = tempInfo['restart'] and os.path.exists('%s/upscatterLegendre.out'%workDir)
    
    if not restart:
        projectileGroupBoundaries = tempInfo['groupBoundaries']
        s = '\n'.join(transferMatricesModule.GBToString('Projectile', projectileGroupBoundaries, None).split()[8:])  ### exclude some header stuff we dont want to parse in C
        dataFile = open('%s/groupStructure.dat'%workDir, 'w')
        dataFile.write(s)
        dataFile.close()

        ### run the calcUpscatter code -- ./calcUpscatterKernel 42092 3.1e-8 ==> upscatterEMuEp.out
        cmd = os.path.join(os.path.split(os.path.dirname(fudgeModule.__file__))[0],
                           'bin', 'calcUpscatterKernel')
        cmd += '  %s %g %g %g %g %d ' % (
            'groupStructure.dat', tempInfo['targetMassRatio'], tempInfo['temperature'],
            tempInfo['minEval'], tempInfo['maxEval'], tempInfo['legendreOrder'])
        logFile = executeCommand(cmd, workDir, 'legendre')
    else:
        logFile = "%s/upscatter_legendre.log" % workDir

    ### call the transferMatrix parseOutput on the outfile
    TM1, TME, paras = transferMatricesModule.parseOutputFile(
        '%s/upscatterLegendre.out'%workDir, 0, firstLine = "upscatter: version 1 \n")
    
    if TM1 is not None:
        negative_l0Counter = checkNegative_l0(TM1, "0", logFile)
    if TME is not None:
        negative_l0Counter += checkNegative_l0(TME, "E", logFile)
    if negative_l0Counter > 0:
        print('WARNING: %d negative l=0 elements found in transfer matrix' % negative_l0Counter)
    APE = parseAveEnergyOutputFile('%s/AveEnergy.out'%workDir)
    if not restart and not tempInfo.get("keepDoubleDifferential"):
        os.remove('%s/upscatterEMuEp.out'%workDir)

    return TM1, TME, APE, paras['maxIncidentEnergyGroup']
        
def executeCommand(cmd, workDir, workFile):
    """
    This functions sets up stuff and executes *cmd* using :py:func:`ubprocessingModule.executeCommand`..

    :param cmd:         The system command to execute.
    :param workDir:     The directory where input and output files are written to.
    :param workFile:    Path to file where infomation about the calculations is written.

    :returns:           The number of negative l=0 product matrix elements.
    """

    transferMatrixExecute = cmd.split()[0]
    if os.path.exists(transferMatrixExecute):
        srcPath = os.path.abspath('./')
    else:
        transferMatrixExecute = os.path.split(transferMatrixExecute)[1]
        srcPath = os.path.join(sys.prefix, 'bin', transferMatrixExecute)
        if os.path.exists(srcPath):
            transferMatrixExecute = srcPath
        else:
            srcPath = os.path.join(os.path.dirname(os.path.dirname(fudgeModule.__file__)), 'bin')
            transferMatrixExecute = os.path.join(srcPath, transferMatrixExecute)

    cmd = [transferMatrixExecute] + cmd.split()[1:]

    if not os.path.exists(workDir):
        os.mkdir(workDir)
    current_directory = os.getcwd()
    os.chdir(workDir)
    fullFileName = 'upscatter_%s' % workFile 
    infoFile = '%s.log' % fullFileName
    
    t0 = timesModule.Times()
    try:
        fCmd = open('%s.sh' % workFile, 'w')
        fCmd.write(' '.join(cmd) + '\n')
        fCmd.close()
        status, stdout, stderr = subprocessingModule.executeCommand(
            cmd, stdout = infoFile , stderr = subprocess.STDOUT)
    except:
        fErr = open(fullFileName + ".err", "w")
        fErr.close()
        raise
        
    fOut = open(infoFile, 'a')
    fOut.write(str(t0) + "\n")
    if status != 0:
        print("status = ", status)
        raise Exception('Upscatter failed for %s %s' % (cmd, fullFileName))
    fOut.close()
    
    os.chdir(current_directory)
    return os.path.join(workDir, infoFile)

### checker for resulting legendre matrices
def checkNegative_l0(TM_EEpL, weight, infoFile):
    """
    This function counts the number of negative l=0 elements in the product matrix *TM_EEpL*.

    :param TM_EEpL:     The product matrix to check.
    :param weight:      A string indicating if the weight is '1' or 'E'. Oops, need to fix the calling string of '0' should be '1'.
    :param infoFile:    Path to log file.

    :returns:           int that is the number of negative l=0 elements found.
    """

    fOut = open(infoFile, 'a')
    negative_l0Counter, largestNegative = 0, 0.
    for ig in sorted(TM_EEpL.keys()):
        TM_EpL = TM_EEpL[ig]
        for ih in sorted(TM_EpL.keys()):
            l0Cell = TM_EpL[ih][0]
            if l0Cell < 0:
                negative_l0Counter += 1
                largestNegative = min(largestNegative, l0Cell)
                fOut.write('negative l=0 for weight %s at row %3d column %3d: %s\n' %
                           (weight, ig, ih, l0Cell))
    fOut.close()
    return negative_l0Counter


def parseAveEnergyOutputFile(fileName):
    """
    For interal use. This function parses the average product data in the file *fileName* and returns the values as a list of floats.

    :param fileName:        The path of the file to parse.

    :returns:               A python list of floats.
    """
    try:
        f = open(fileName)
    except:
        raise Exception('Could not open file = %s' % fileName)
    ls = f.readlines()
    f.close()

    APE = []
    for s in ls:
        APE.append([float(x) for x in s.split()[-2:]])

    return APE    
        


def rescaleCrossSection(groupedCrossSec, UpScatterMatrix):
    """
    For each incident group in *UpScatterMatrix*, this function scales each element so that sum of the l=0 component is equal to 
    the corresponding element in *groupedCrossSec*.

    :param groupedCrossSec:     The multi-group cross sections.
    :param UpScatterMatrix:     The product matrices to scale to match *groupedCrossSec*.

    :returns:                   The scaled *UpScatterMatrix*.
    """

    for g,sigma in enumerate(groupedCrossSec):
        summedarray = [UpScatterMatrix[g][h][0] for h in range(len(UpScatterMatrix[g]))]
        sumXS = sum(summedarray) 
        if sumXS != 0.0:
            ratio = sigma/sumXS
            for l in range(len(UpScatterMatrix[g][0])):
                for h in range(len(UpScatterMatrix[g])) :
                    if (UpScatterMatrix[g][h][l] != 0.0):
                        UpScatterMatrix[g][h][l] = UpScatterMatrix[g][h][l]*ratio

    return UpScatterMatrix
