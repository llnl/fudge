# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import subprocess
import os
import re

# NOTE: 1) statusMessageReporting test data is ignored since it is only applicable to C code


def process_args():
    # see http://docs.python.org/lib/optparse-tutorial.html
    import argparse
    parser = argparse.ArgumentParser(description="Download the data necessary to run 'make check' for fudge/Makefile")
    parser.add_argument('fudgeRepoURL', metavar='fudgeRepoURL', type=str, help='URL to the FUDGE GIt repo, e.g. ssh://.../fudge/fudge.git@commitLabel')
    parser.add_argument('requiresFile', metavar='requiresFile', type=str, help='Path to requires.txt stored in the FUDGE egg-info folder')

    return parser.parse_args()

def archiveFromGIT(testData, repoURL):
    regexURL = re.compile(r'^([^\@]+\@[^\@]+)\@(.+)$')
    if regexURL.match(repoURL):
        repoURLNoTag, repoTag = regexURL.findall(repoURL)[0]
    else:
        repoURLNoTag, repoTag = (repoURL, 'HEAD')

    gitCommand = ['git',  'archive',  '--remote=%s' % repoURLNoTag, repoTag] + testData

    return subprocess.Popen(gitCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)

def urlFromGitModules(repoURL):
    processGitArchive = archiveFromGIT([r'.gitmodules'], repoURL)

    tarCommand = ['tar', '--to-stdout', '-xf', '-']
    processTar = subprocess.Popen(tarCommand, stdin=processGitArchive.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)

    processGitArchive.stdout.close()
    stdout, stderr = [output.decode('utf-8') for output in processTar.communicate()]
    assert stderr == '', 'Shell command "%s" failed:\n%s' % (' '.join(tarCommand), stderr)

    return re.findall(r'^\s*url\s*=\s*(ssh:\S+)\s*$', stdout, re.MULTILINE)[0]

def downloadTestData(testData, repoURL, localDestination):
    processGitArchive = archiveFromGIT(testData, repoURL)

    if localDestination:
        os.makedirs(localDestination, exist_ok=True)
        tarCommand = ['tar', '-C', localDestination, '-xf', '-']
    else:
        tarCommand = ['tar', '-xf', '-']

    processTar = subprocess.Popen(tarCommand, stdin=processGitArchive.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)

    processGitArchive.stdout.close()
    stdout, stderr = processTar.communicate()
    assert stderr.decode('utf-8') == '', 'Shell command "%s" failed:\n%s' % (' '.join(tarCommand), stderr.decode('utf-8'))


args = process_args()

# package REPO URLs
packageRepoURL = {'fudge': args.fudgeRepoURL}
targetPath = {'fudge': None}
with open(args.requiresFile) as fileObject:
    requiresLine = fileObject.readline()
    while requiresLine:
        if re.match(r'^git', requiresLine):
            gitURL, packageName = re.findall(r'^([^#]+)#egg=(\S+)', requiresLine)[0]
            packageRepoURL[packageName] = gitURL.replace(r'git+', '')
            targetPath[packageName] = packageName

        requiresLine = fileObject.readline()

# # URL for statusMessageReporting (git submodule of numericalFunctions)
# packageName = 'statusMessageReporting'
# packageRepoURL[packageName] = urlFromGitModules(packageRepoURL['numericalFunctions'])
# targetPath[packageName] = packageName


testData = {
             # fudge 
             'fudge': 
                [
                 # root level Makefile
                 'Makefile',

                 # rebuild_test_data
                 r'fudge/covariances/test/rebuild_test_data.py', r'fudge/covariances/test/*.endf',
                 r'fudge/processing/resonances/test/rebuild_test_data.py', r'fudge/processing/resonances/test/*.endf',

                 # fudge test files python scripts
                 r'fudge/core/math/test/test*.py', r'fudge/covariances/test/test*.py',
                 r'fudge/processing/resonances/test/test*.py',
                 r'fudge/productData/distributions/test/__init__.py', r'fudge/reactionData/test/test_crossSection.py'],

             # brownies
             'brownies': [r'legacy/endl/test/test_endlProject.py', r'legacy/endl/test/testdb'],

             # xData
             'xData': [r'test/test_XYs.py'],

             # pqu
             'pqu': [r'Check/*.py', r'Check/Out.checked/*.out'],

             # numericalFunctions
             'numericalFunctions': [
                 # scripts
                 r'nf_specialFunctions/Python/Test/UnitTesting/',

                 #data
                 r'nf_specialFunctions/Test/UnitTesting/exponentialIntegral/exponentialIntegralTest.dat',
                 r'nf_specialFunctions/Test/UnitTesting/gammaFunctions/gammaTest.dat',
                 r'nf_specialFunctions/Test/UnitTesting/gammaFunctions/incompleteGammaTest.dat']
            }

for packageName in testData.keys():
    downloadTestData(testData[packageName], packageRepoURL[packageName], targetPath[packageName])



