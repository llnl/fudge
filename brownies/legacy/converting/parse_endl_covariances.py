# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
Function that extracts covariance data from an ENDL *cov.xml file

For example,
python parse_endl_covariances.py /usr/gapps/data/nuclear/endl_official/endl2009.3/ascii/yi01/za094239/yo00c15i000s000_cov.xml
"""
import sys
import os
import numpy
from xml.etree import ElementTree

def covarianceDict(rootdir):  #  ACD ADDED THIS CHUNK
    """
    Creates a dictionary of covariance files in rootdir, keyed by the pair (C,S).
    """
    covar_dictionary = {}
    for folder, dirs, files in os.walk(rootdir):
        for covFile in files:
            if covFile.endswith('_cov.xml'):
                C = int(covFile.split('c')[1].split('i')[0])
                S = int(covFile.split('s')[1].split('_')[0])
                if (C,S) not in covar_dictionary:
                    covar_dictionary[(C,S)] = []
                covar_dictionary[(C,S)].append(covFile)

    return ( covar_dictionary )


def parse_endl_covariance(covFile):
    """
    Returns a list of (energy bins,  covariance matrix, covariance_type) tuples.
    The list may be empty (some ENDL cov.xml files are empty).
    """
    xdoc = ElementTree.parse(covFile)
    root = xdoc.getroot()
    if root.tag not in ('cross_section_covariance'):
        # to-be-done: support energy distribution covariances like za094239/yo01c15i005s000_cov.xml
        raise NotImplementedError("Unsupported covariance type %s" % root.tag)

    ebins, covariances, covariance_types = [],[],[]
    for child in root:
        if child.tag == 'covariance':
            covariance_types.append( child.get('type') )
            matrixNode = child.find('matrix')
            rows = int(matrixNode.get('dim1'))
            columns = int(matrixNode.get('dim2'))
            covariances.append(
                    numpy.array(list(map(float, matrixNode.text.split()))).reshape(rows,columns) )
        elif child.tag == 'histogram':
            ebins.append( numpy.array( list(map(float, child.text.split())) ) )

    return list(zip(ebins, covariances, covariance_types))


if __name__ == '__main__':
    results = parse_endl_covariance(sys.argv[1])
    print("Read %d covariances from file %s" % (len(results), sys.argv[1]))

