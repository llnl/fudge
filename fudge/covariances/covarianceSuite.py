# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Module with containers for covariances in several different forms
"""

import sys
import os
import pathlib

from LUPY import ancestry as ancestryModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import axes as axesModule

from fudge import enums as enumsModule
from fudge import styles as stylesModule, suites as suitesModule
from . import covarianceSection as sectionModule, summed as summedModule, mixed as mixedModule, \
    modelParameters as modelParametersModule


class CovarianceSections(suitesModule.Suite):
    """
    Most covariances are stored in 'sections', each representing either the internal covariance
    for a single quantity (i.e. cross section) or the cross-term between two quantities
    """

    moniker = 'covarianceSections'
    legacyMemberNameMapping = {'section': sectionModule.CovarianceSection.moniker}

    def __init__(self):
        suitesModule.Suite.__init__(self, (sectionModule.CovarianceSection,))


class ParameterCovariances(suitesModule.Suite):
    """
    Resolved and unresolved resonance parameters are stored inside the parameterCovariances.
    """

    moniker = 'parameterCovariances'
    legacyMemberNameMapping = {'section': modelParametersModule.AverageParameterCovariance.moniker}

    def __init__(self):
        suitesModule.Suite.__init__(self, (modelParametersModule.ParameterCovariance,
                                           modelParametersModule.AverageParameterCovariance))


class CovarianceSuite(ancestryModule.AncestryIO):
    """
    All covariances for a target/projectile combination are stored in a :py:class:`CovarianceSuite` in gnds.
    The :py:class:`CovarianceSuite` is stored in a separate file from the reactionSuite.

    Within the :py:class:`CovarianceSuite`, data is sorted into :py:class:`sections` (see the section module), each
    of which contains one section of the full covariance matrix.
    
    """

    moniker = 'covarianceSuite'
    ancestryMembers = ('externalFiles', 'styles', 'covarianceSections', 'parameterCovariances')
    childNodeOrder = {
        GNDS_formatVersionModule.version_1_10: (suitesModule.ExternalFiles.moniker, stylesModule.Styles.moniker,
                                                CovarianceSections.moniker, ParameterCovariances.moniker),
        GNDS_formatVersionModule.version_2_0_LLNL_4: (suitesModule.ExternalFiles.moniker, stylesModule.Styles.moniker,
                                                      CovarianceSections.moniker, ParameterCovariances.moniker),
        GNDS_formatVersionModule.version_2_0: (suitesModule.ExternalFiles.moniker, stylesModule.Styles.moniker,
                                               CovarianceSections.moniker, ParameterCovariances.moniker)}

    def __init__(self, projectile, target, evaluation, interaction=None,
                 formatVersion=GNDS_formatVersionModule.default, sourcePath=None):
        """
        :param projectile:  particle id
        :param target:      particle id
        :param evaluation:  evaluation id
        :param interaction: type of interaction (typically 'nuclear')
        """

        ancestryModule.AncestryIO.__init__(self)

        if formatVersion not in GNDS_formatVersionModule.allowedPlus: raise Exception(
            "Unsupported GNDS structure '%s'!" % str(formatVersion))
        self.formatVersion = formatVersion

        self.projectile = projectile  #: The projectile
        self.target = target  #: The target
        self.evaluation = evaluation

        if interaction == enumsModule.Interaction.legacyTNSL:
            interaction = enumsModule.Interaction.TNSL
        self.interaction = interaction

        self.sourcePath = sourcePath

        self.__externalFiles = suitesModule.ExternalFiles()
        self.__externalFiles.setAncestor(self)

        self.__styles = stylesModule.Styles()
        self.__styles.setAncestor(self)

        self.__covarianceSections = CovarianceSections()
        self.__covarianceSections.setAncestor(self)

        self.__parameterCovariances = ParameterCovariances()
        self.__parameterCovariances.setAncestor(self)

    @property
    def sourcePath(self):
        """Returns the sourcePath member which is the path to the covarianceSuite file for self if self is from a file."""

        return (self.__sourcePath)

    @sourcePath.setter
    def sourcePath(self, path):
        """
        This method sets *self*'s *sourcePath* to *path*.

        :param path:        The new *sourcePath*.
        """

        if isinstance(path, pathlib.Path):
            path = str(path)
        if path is not None:
            if not isinstance(path, str):
                raise ValueError('Path must be a str or None, got "%s"' % type(path))

        self.__sourcePath = path

    @property
    def styles(self):
        return self.__styles

    @property
    def externalFiles(self):
        return self.__externalFiles

    @property
    def covarianceSections(self):
        return self.__covarianceSections

    @property
    def parameterCovariances(self):
        return self.__parameterCovariances

    @property
    def domainMin(self):
        """Returns the minimum of the projectile energy for the evaluation. This needs to be fixed to handle multiple evaulation styles."""

        return (self.styles[0].projectileEnergyDomain.min)

    @property
    def domainMax(self):
        """Returns the maximum of the projectile energy for the evaluation. This needs to be fixed to handle multiple evaulation styles."""

        return (self.styles[0].projectileEnergyDomain.min)

    @property
    def domainUnit(self):
        """Returns the unit of the projectile energy for the evaluation. This needs to be fixed to handle multiple evaulation styles."""

        return (self.styles[0].projectileEnergyDomain.unit)

    @property
    def interaction(self):
        """Returns self's interaction."""

        return (self.__interaction)

    @interaction.setter
    def interaction(self, value):

        if value is None:
            value = enumsModule.Interaction.nuclear
            if self.formatVersion != GNDS_formatVersionModule.version_1_10:
                print('Need to specify interaction when calling CovarianceSuite.__init__. Setting it to "%s". Please update your code as in the future this will execute a raise.' % value)

        self.__interaction = enumsModule.Interaction.checkEnumOrString(value)

    def convertUnits(self, unitMap):
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        for style in self.styles: style.convertUnits(unitMap)
        for section in self.covarianceSections: section.convertUnits(unitMap)
        for mpsection in self.parameterCovariances: mpsection.convertUnits(unitMap)

    def check(self, **kwargs):
        """
        Check all covariance sections, returning a list of warnings. 
        
        :keyword bool checkUncLimits: Should we check the uncertainty limits? (default: True)
        :keyword float minRelUnc: Minimum allowable relative uncertainty (default: 0.0)
        :keyword float maxRelUnc: Maximum allowable relative uncertainty (default: 10.0)
        :keyword theData: A reference to the data for this covariance.
                          This is useful for converting between relative and absolute covariance (default: None)
        :keyword float negativeEigenTolerance: Ignore negative eigenvalues smaller than this (default: -1e-6) 
        :keyword float eigenvalueRatioTolerance: Warn if smallest eigenvalue < this value * biggest (default: 1e-8)
        :keyword float eigenvalueAbsoluteTolerance: Warn if smallest eigenvalue < this value (default: 1e-14)
        :keyword bool verbose: Be verbose while running checks
        :rtype: warning.Context
        """

        from fudge import warning

        # default input options
        options = {
            'checkUncLimits': True,
            'minRelUnc': 0.0,
            'maxRelUnc': 10.0,
            'theData': None,
            'negativeEigenTolerance': -1e-6,  # ignore smaller negative eigenvalues
            'eigenvalueRatioTolerance': 1e-8,  # warn if smallest eival < 1e-8 * biggest
            'eigenvalueAbsoluteTolerance': 1e-14,
            'verbose': False,
        }
        for key in kwargs:
            if key in options:
                options[key] = kwargs[key]
            else:
                raise KeyError("check() received unknown keyword argument '%s'" % key)

        # assemble some useful info, to be handed down to children's check() functions:
        info = {'covarianceSuite': self, 'style': 'eval'}
        info.update(options)

        warnings = []

        #: summedReactions are sums of other sections, but that opens the possibility of a cyclic dependency
        #: check for that, using algorithm taken from
        #:    http://neopythonic.blogspot.com/2009/01/detecting-cycles-in-directed-graph.html 
        def find_cycle(NODES, EDGES, style):
            todo = set(NODES)
            while todo:
                node = todo.pop()
                stack = [node]
                while stack:
                    top = stack[-1]
                    for node in EDGES(top, style):
                        if node in stack:
                            return stack[stack.index(node):]
                        if node in todo:
                            stack.append(node)
                            todo.remove(node)
                            break
                    else:  # this node is not in a cycle
                        stack.pop()
            return None

        def get_edges(section_, style):  #: return list of all pointers from this section
            form = section_[style]
            if isinstance(form, summedModule.SummedCovariance):
                return [v.link for v in form]
            elif isinstance(form, mixedModule.MixedForm):
                edges = []
                for part in form:
                    if isinstance(part, summedModule.SummedCovariance): edges += [v.link for v in part]
                return edges

        nodes = [sec for sec in self.covarianceSections if
                 get_edges(sec, info['style'])]  # sections that contain pointers

        if nodes:
            cycle = find_cycle(nodes, get_edges, info['style'])
            if cycle:
                warnings.append(warning.CyclicDependency(cycle))

        # check each section
        for section_ in self.covarianceSections:
            sectionWarnings = section_.check(info)
            if sectionWarnings:
                warnings.append(warning.Context("covarianceSection '%s':" % section_.label, sectionWarnings))

        for parameterCovariance in self.parameterCovariances:
            pcwarnings = parameterCovariance.check(info)
            if pcwarnings:
                warnings.append(warning.Context("parameterCovariance '%s':" % parameterCovariance.label, pcwarnings))

        return warning.Context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)

    def fix(self, **kwargs):
        """
        Apply basic fixes to a covariance

        :param bool removeNegativeEVs: Should we remove eigenspaces corresponding to negative eigenvalues? (Default: True)
        :param bool removeSmallEVs: Should we remove eigenspaces corresponding to small eigenvalues? (Default: False)
        :param bool fixUncLimits: Should we fix the uncertainties to lie within bounds imposed by minRelUnc and maxRelUnc? (Default: False)
        :param minRelUnc: Minimum allowable uncertainty for this covariance
        :type minRelUnc: float or None
        :param maxRelUnc: Maximum allowable uncertainty for this covariance
        :type maxRelUnc: float or None
        :param theData: Reference to the data that accompanies this covariance so that we may convert between absolute and relative covariance as needed.
        :type theData: instance or None 

        """
        from fudge import warning
        # default input options
        options = {
            'removeNegativeEVs': True,
            'removeSmallEVs': False,
            'fixUncLimits': False,
            'minRelUnc': None,
            'maxRelUnc': None,
            'theData': None,
            'negativeEigenTolerance': -1e-6,  # ignore smaller negative eigenvalues
            'eigenvalueRatioTolerance': 1e-8,  # warn if smallest eival < 1e-8 * biggest
            'eigenvalueAbsoluteTolerance': 1e-14,
        }
        for key in kwargs:
            if key in options:
                options[key] = kwargs[key]
            else:
                raise KeyError("fix() received unknown keyword argument '%s'" % key)
        # assemble some useful info, to be handed down to children's check() functions:
        info = {'covarianceSuite': self}
        info.update(options)
        # do the fixing
        warnings = []
        for section_ in self.covarianceSections: warnings += section_.fix(**info)
        return warning.Context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)

    def removeExtraZeros(self):
        """
        Checks all covariance matrices for rows/columns of all zero, removes them if found.
        """
        for section_ in self.covarianceSections:
            for form in section_:
                if hasattr(form, "removeExtraZeros"):
                    form.removeExtraZeros()
                if isinstance(form, mixedModule.MixedForm):
                    for covar in form:
                        if hasattr(covar, "removeExtraZeros"):
                            covar.removeExtraZeros()

    def toXML_strList(self, indent='', **kwargs):
        """Returns a list of GNDS/XML strings representing self."""

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        if formatVersion in (GNDS_formatVersionModule.version_2_0_LLNL_3, GNDS_formatVersionModule.version_2_0_LLNL_4):
            print('INFO: converting GNDS format from "%s" to "%s".' % (
            formatVersion, GNDS_formatVersionModule.version_2_0))
            formatVersion = GNDS_formatVersionModule.version_2_0
        kwargs['formatVersion'] = formatVersion
        if formatVersion not in GNDS_formatVersionModule.allowed:
            raise Exception("Unsupported GNDS structure '%s'!" % str(formatVersion))

        interaction = self.interaction
        if interaction == enumsModule.Interaction.TNSL:
            interaction = enumsModule.Interaction.getTNSL_interaction(formatVersion)
        xmlString = ['%s<%s projectile="%s" target="%s" evaluation="%s" interaction="%s" format="%s">'
                     % (
                     indent, self.moniker, self.projectile, self.target, self.evaluation, interaction, formatVersion)]
        xmlString += self.externalFiles.toXML_strList(indent2, **kwargs)
        xmlString += self.styles.toXML_strList(indent2, **kwargs)
        xmlString += self.covarianceSections.toXML_strList(indent2, **kwargs)
        xmlString += self.parameterCovariances.toXML_strList(indent2, **kwargs)
        xmlString.append('%s</%s>' % (indent, self.moniker))

        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)  # keep track of location in the tree, in case errors come up

        linkData = {'reactionSuite': kwargs.get('reactionSuite'), 'unresolvedLinks': []}

        sourcePath = kwargs.get('sourcePath')
        try:
            formatVersion = element.get('format')
            linkData['formatVersion'] = formatVersion
            CS = cls(element.get('projectile'), element.get('target'), element.get('evaluation'),
                     element.get('interaction'), formatVersion, sourcePath=sourcePath)
            for child in element:
                if child.tag == 'styles':
                    CS.styles.parseNode(child, xPath, linkData, **kwargs)
                elif child.tag == suitesModule.ExternalFiles.moniker:
                    CS.externalFiles.parseNode(child, xPath, linkData, **kwargs)
                elif child.tag == CovarianceSections.moniker:
                    CS.covarianceSections.parseNode(child, xPath, linkData, **kwargs)
                elif child.tag == ParameterCovariances.moniker:
                    CS.parameterCovariances.parseNode(child, xPath, linkData, **kwargs)
                else:
                    print("Warning: encountered unexpected element '%s' in covarianceSuite!" % child.tag)
        except Exception:
            print('Error encountered at xpath "%s" while parsing %s.' % ('/'.join(xPath), cls.moniker))
            raise

        # fix links:
        for link_ in linkData['unresolvedLinks']:
            if isinstance(link_.ancestor, axesModule.Grid) and link_.path[0:3] == '../': continue
            path = link_.path
            if link_.root is None:
                if path.startswith('/covarianceSuite'):
                    root = CS
                elif path.startswith('.'):
                    root = link_  # relative link
            elif link_.root == '$reactions':
                root = linkData['reactionSuite']
            else:
                root = None

            if root is None:
                continue  # FIXME links to other files (i.e. cross-material covariances) currently broken
            try:
                res = link_.follow(root)
                link_.link = res
            except ancestryModule.XPathNotFound:
                print("Warning: link %s could not be resolved!" % link_)

        xPath.pop()

        return CS

    def parseCleanup(self, node, **kwargs):

        reactionSuite = kwargs.get('reactionSuite')
        if reactionSuite is not None:
            errors = 0
            for externalLink in reactionSuite._externalLinks:
                try:
                    externalLink.link = externalLink.follow(self)
                except ancestryModule.XPathNotFound:
                    errors += 1
                    print("Warning: link %s could not be resolved!" % externalLink.path)

            if errors:
                raise ancestryModule.XPathNotFound("Encountered %d broken xpath links!" % errors)

        return CovarianceSuite

    @staticmethod
    def read(fileName, **kwargs):
        """
        Reads in the file name *fileName* and returns a **ReactionSuite** instance.
        """

        return CovarianceSuite.readXML_file(fileName, **kwargs)


def read(fileName, **kwargs):
    """
    Reads in the file name *fileName* and returns a **ReactionSuite** instance.
    """

    return CovarianceSuite.read(fileName, **kwargs)
