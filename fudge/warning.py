# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import functools

import PoPs.warning
from fudge import styles as stylesModule
from LUPY import enums as enumsModule

"""
Store and report warnings and errors in GNDS data, in order to discover problems in the library.
The reactionSuite.check() function returns a nested list of Warning objects:

    >>> warnings = reactionSuite.check()
    >>> print(warnings)

May include or exclude specific classes of Warning using the filter command.
filter() returns a new context instance:

    >>> warnings2 = warnings.filter( exclude=[Warning.EnergyImbalance] )

Or, for easier searching you may wish to flatten the list (to get warnings alone without context messages):

    >>> flat = warnings.flatten()
"""


@functools.total_ordering
class Level(enumsModule.Enum):
    """
    Define severity levels for warnings.
    """
    Pedantic = enumsModule.auto()
    Minor = enumsModule.auto()
    Moderate = enumsModule.auto()
    Severe = enumsModule.auto()
    Fatal = enumsModule.auto()

    def __le__(self, other):
        if not isinstance(other, Level):
            return NotImplemented

        return self._positions[self] < self._positions[other]


Level._positions = {x: i for i, x in enumerate(Level)}


class Context:
    """
    Store warnings in Context. This class contains location information (reactionSuite, reaction, etc)
    plus a nested list of warnings or other Context instances
    """

    def __init__(self, message='', warningList=None):
        self.message = message
        self.warningList = warningList or []

    def __len__(self):
        return len(self.warningList)

    def __getitem__(self, idx):
        return self.warningList[idx]

    def __str__(self):
        if len(self.warningList) == 0:
            return self.message + ": no problems encountered"
        return '\n'.join(self.toStringList())

    def __eq__(self, other):
        return self.message == other.message and self.warningList == other.warningList

    def filter(self, threshold=None, include=None, exclude=None):
        """
        Filter warning list to only include warnings at or greater than specified threshold,
        or to include (or exclude) specific classes of Warning.
        Returns new Context instance with screened list, plus dictionary indicating how many warnings were screened.
        Examples:

            >>> context = Context()
            >>> newWarnings, screened = context.filter( threshold=Level.Moderate )
            # or
            >>> newWarnings, screened = context.filter( exclude=[Warning.EnergyImbalance, Warning.Q_mismatch] )

        'include' takes precedence over 'exclude', 'threshold' can be used with either 'include' or 'exclude'
        """

        newWarningList = []
        screened = {}
        if threshold is None and include is None and exclude is None:
            return self, screened
        for warning in self.warningList:
            if isinstance(warning, (Context, PoPs.warning.Context)):
                newContext, newScreened = warning.filter(threshold, include, exclude)
                if newContext: newWarningList.append(newContext)
                for key in newScreened:
                    screened[key] = screened.get(key, 0) + newScreened[key]
            elif threshold is not None:
                if warning.level >= threshold:
                    if include is not None:
                        if warning.__class__ in include:
                            newWarningList.append(warning)
                        else:
                            screened[warning.__class__] = screened.get(warning.__class__, 0) + 1
                    elif exclude is not None:
                        if warning.__class__ not in exclude:
                            newWarningList.append(warning)
                        else:
                            screened[warning.__class__] = screened.get(warning.__class__, 0) + 1
                    else:
                        newWarningList.append(warning)
                else:
                    screened[warning.level] = screened.get(warning.level, 0) + 1
            elif include is not None:
                if warning.__class__ in include:
                    newWarningList.append(warning)
                else:
                    screened[warning.__class__] = screened.get(warning.__class__, 0) + 1
            else:  # exclude is not None:
                if warning.__class__ not in exclude:
                    newWarningList.append(warning)
                else:
                    screened[warning.__class__] = screened.get(warning.__class__, 0) + 1
        return Context(self.message, newWarningList), screened

    def flatten(self):
        """From a nested hierarchy of warnings, get back a flat list for easier searching:

            >>> w = reactionSuite.check()
            >>> warningList = w.flatten()"""

        List = []
        for val in self.warningList:
            if isinstance(val, Warning):
                List.append(val)
            else:
                List += val.flatten()
        return List

    def toStringList(self, indent='', dIndent='    '):
        """ Format warnings for printing. Returns a list of warning strings with indentation. """
        s = ['%s%s' % (indent, self.message)]
        for warning in self.warningList:
            s += warning.toStringList(indent + dIndent)
        return s


class DiffResult:

    def __init__(self, code, message, link1, link2):
        self.code = code
        self.message = message
        self.link1 = link1
        self.link2 = link2

    def __repr__(self):
        return "%s: %s\n    %s\n    %s" % (self.code, self.message, self.link1, self.link2)


class DiffResults:

    def __init__(self, protare1, protare2):

        self.protare1 = protare1
        self.protare2 = protare2
        self.diffs = []

    def __len__(self):
        """Returns the number of differences found."""

        return len(self.diffs)

    def __iter__(self):

        n1 = len(self.diffs)
        for i1 in range(n1): yield self.diffs[i1]

    def __repr__(self):
        """Returns a string representation of the differences found."""

        return '\n'.join([str(diff) for diff in self.diffs])

    def __copy__(self):
        """Make a copy of self."""

        diffResults1 = DiffResults(self.protare1, self.protare2)
        for item in self: diffResults1.append(DiffResult(item.code, item.message, item.link1, item.link2))

        return diffResults1

    def append(self, code, message, link1, link2):
        """Adds a DiffResult instance to self."""

        self.diffs.append(DiffResult(code, message, link1, link2))

    def cull(self, code):
        """Removes all DiffResult instance that contain the specified code."""

        cullList = []
        for i1, item in enumerate(self.diffs):
            if item.code == code: cullList.append(i1)
        cullList.reverse()
        for i1 in cullList: del self.diffs[i1]

    def patch(self, label, version, library=None):

        chains = self.protare2.styles.chains(True)
        if len(chains) != 1: raise Exception("Only one styles chain is supported.")

        chain = chains[0]
        for style in chain:
            if isinstance(style, stylesModule.Evaluated): break
        if not isinstance(style, stylesModule.Evaluated):
            raise Exception("At least one evaluated style must be present in file to patch.")

        if len(self.diffs) > 0:
            if label in self.protare2.styles: raise ValueError("Label already in file to patch.")
            if library is None: library = style.library
            evaluateStyle = stylesModule.Evaluated(label, style.label, None, None, library, version)
            if style.library == evaluateStyle.library:
                if style.version == evaluateStyle.version:
                    raise ValueError("Library and version must be different than one in patch file.")
            self.protare2.styles.add(evaluateStyle)

            for diff in self:
                if diff.code == "Distribution unspecified - 1":
                    pass
                elif diff.code == "Distribution unspecified - 2":
                    fromComponent = self.protare1.followXPath(diff.link1)
                    toComponent = self.protare2.followXPath(diff.link2)
                    for key in toComponent.keys(): toComponent.remove(key)
                    toComponent.add(fromComponent.evaluated)
                elif diff.code == "reactions missing":
                    if diff.link1 != "":
                        fromReaction = self.protare1.followXPath(diff.link1)
                        fromReaction.amendForPatch(style.label, evaluateStyle.label)
                        self.protare2.reactions.add(fromReaction)
                else:
                    print("""Unsupported diff code "%s".""" % diff.code)

        crossSectionReconstructedStyles = self.protare2.styles.getStylesOfClass(stylesModule.CrossSectionReconstructed)
        if len(crossSectionReconstructedStyles) > 0:
            if len(crossSectionReconstructedStyles) > 1:
                raise Exception("Multiple cross section reconstructed styles present which is not supported.")
            crossSectionReconstructedStyle = crossSectionReconstructedStyles[0]
            self.protare2.removeStyle(crossSectionReconstructedStyle.label)
            crossSectionReconstructedStyle = stylesModule.CrossSectionReconstructed(
                crossSectionReconstructedStyle.label, label)
            self.protare2.reconstructResonances(crossSectionReconstructedStyle, 1e-3)


class Warning(BaseException):
    """
    General Warning class. Contains link to problem object,
    xpath in case the object leaves memory,
    and information about the warning or error.
    """

    def __init__(self, level: Level, obj=None):
        self.level = level
        self.obj = obj
        self.xpath = ''
        if hasattr(obj, 'toXLink'):
            self.xpath = obj.toXLink()

    def __str__(self):
        return "Generic warning for %s" % self.xpath

    def __eq__(self, other):
        return self.xpath == other.xpath

    def toStringList(self, indent=''):
        return ['%s%s warning: %s' % (indent, self.level, self)]


#
# specific Warning classes:
#

class NotImplemented(Warning):
    def __init__(self, form, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.form = form

    def __str__(self):
        return "Checking not yet implemented for %s type data" % self.form

    def __eq__(self, other):
        return self.form == other.form and self.xpath == other.xpath


class UnorthodoxParticleNotImplemented(Warning):
    def __init__(self, particleId, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.particleId = particleId

    def __str__(self):
        return "Checking not yet implemented for unorthodox particle '%s'" % self.particleId

    def __eq__(self, other):
        return self.particleId == other.particleId


class ElementalTarget(Warning):
    def __init__(self, particleId, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.particleId = particleId

    def __str__(self):
        return "Elemental target '%s'!  Some checks (e.g. energy balance and Q-values) will be skipped" % (
            self.particleId)

    def __eq__(self, other):
        return self.particleId == other.particleId


# external files:

class MissingExternalFile(Warning):
    def __init__(self, path, obj=None):
        Warning.__init__(self, Level.Fatal, obj)
        self.path = path

    def __str__(self):
        return "External file '%s' is missing" % self.path


class WrongExternalFileChecksum(Warning):

    def __init__(self, path, expected, computed, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.path = path
        self.expected = expected
        self.computed = computed

    def __str__(self):
        s = "Computed checksum for external file '%s' does not match: %s vs. %s" % (
            self.path, self.computed, self.expected
        )
        return s


class EvaluationDomainMinTooHigh(Warning):

    def __init__(self, expected, obj):
        Warning.__init__(self, Level.Severe, obj)
        self.expected = expected

    def __str__(self):
        s = f"Evaluation starts at incident energy {self.obj.domainMin} {self.obj.domainUnit}, expected {self.expected}"
        return s


# resonance region:


class BadScatteringRadius(Warning):

    def __init__(self, factor, gotAP, expectedAP, L=None, E=None, obj=None):

        level = Level.Moderate
        if factor >= 10:
            level = Level.Severe
        Warning.__init__(self, level, obj)
        self.factor = factor
        self.gotAP = gotAP
        self.expectedAP = expectedAP
        self.L = L
        self.E = E

    def __str__(self):
        direction = 'high'
        if self.gotAP < self.expectedAP: direction = 'low'
        s = "Scattering radius (%s) at least a factor of %s too %s compared to expectations (%s)" % (
            self.gotAP, self.factor, direction, self.expectedAP)
        if self.L is not None: s += " for L=%i" % self.L
        if self.E is not None: s += " for E=%s" % str(self.E)
        return s


class BadSpinStatisticalWeights(Warning):

    def __init__(self, L, gJ, expectedgJ, reaction=None):

        diff = abs(abs(gJ) - abs(2.0 * L + 1))
        level = Level.Minor
        if diff > 0.1: level = Level.Moderate
        if diff > 0.5: level = Level.Severe
        Warning.__init__(self, level, reaction)
        self.L = L
        self.gJ = gJ
        self.expectedgJ = expectedgJ
        self.reaction = reaction

    def __str__(self):
        direction = 'many'
        if self.gJ < self.expectedgJ: direction = 'few'
        s = "The spin statistical weights for L=%i sums to %s, but should sum to %s.  You have too %s channels " % (
            self.L, str(self.gJ), str(self.expectedgJ), direction)
        if self.reaction is not None:
            s += 'for reaction ' + self.reaction
        return s


class InvalidAngularMomentaCombination(Warning):

    def __init__(self, L, S, J, name, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.L = L
        self.S = S
        self.J = J
        self.name = name

    def __str__(self):
        return 'Invalid angular momenta combination: cannot couple L = %s and S = %s up to J = %s for channel "%s"' % (
            self.L, self.S, self.J, self.name)


class InvalidParity(Warning):

    def __init__(self, L, S, J, Pi, name, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.L = L
        self.S = S
        self.J = J
        self.Pi = Pi
        self.name = name

    def __str__(self):
        return 'Invalid parity %s for channel "%s" with L=%s, S=%s, J=%s' % (
            self.Pi, self.name, self.L, self.S, self.J)


class UnknownMass(Warning):

    def __init__(self, particle, obj=None):
        Warning.__init__(self, Level.Minor, obj)
        self.particle = particle

    def __str__(self):
        return "Could not determine mass for particle '%s'" % self.particle


class UnknownSpinParity(Warning):

    def __init__(self, reaction, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.reaction = reaction

    def __str__(self):
        return "Could not determine product spin/parities for reaction '%s'" % self.reaction


class MissingResonanceChannel(Warning):

    def __init__(self, L, S, J, name, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.L = L
        self.S = S
        self.J = J
        self.name = name

    def __str__(self):
        return 'Missing a channel with angular momenta combination L = %s, S = %s and J = %s for "%s"' % (
            self.L, self.S, self.J, self.name)


class InvalidSpinCombination(Warning):

    def __init__(self, Ia, Ib, S, name, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.Ia = Ia
        self.Ib = Ib
        self.S = S
        self.name = name

    def __str__(self):
        return 'Invalid spin combination: cannot couple Ia = %s and Ib = %s up to S = %s for channel "%s"' % (
            self.Ia, self.Ib, self.S, self.name)


class PotentialScatteringNotConverged(Warning):

    def __init__(self, L, E, fom, fomTarget, obj=None):
        level = Level.Minor
        if fom > 0.1:
            level = Level.Moderate
        if fom > 0.5:
            level = Level.Severe
        Warning.__init__(self, level, obj)
        self.L = L
        self.E = E
        self.fom = fom
        self.fomTarget = fomTarget

    def __str__(self):
        s = "Potential scattering hasn't converged by L=%d at E=%s eV, xs[%d]/xs[0]=%s%% > %s%%" % (
            self.L, self.E, self.L, str(100.*self.fom), str(100.*self.fomTarget))
        return s


class RRmultipleRegions(Warning):

    def __init__(self, obj=None):
        Warning.__init__(self, Level.Moderate, obj)

    def __str__(self):
        return "Use of more than one resolved resonance regions is deprecated"


class URRdomainMismatch(Warning):

    def __init__(self, Lval, Jval, name, obj=None):
        Warning.__init__(self, Level.Fatal, obj)
        self.Lval = Lval
        self.Jval = Jval
        self.name = name

    def __str__(self):
        return "Unresolved L=%i / J=%.1f %s doesn't span URR energy limits" % (self.Lval, self.Jval, self.name)

    def __eq__(self, other):
        return self.Lval == other.Lval and self.Jval == other.Jval


class URRinsufficientEnergyGrid(Warning):

    def __init__(self, Lval, Jval, eLow, eHigh, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.Lval = Lval
        self.Jval = Jval
        self.eLow = eLow
        self.eHigh = eHigh

    def __str__(self):
        return "More points needed in L=%i J=%.1f unresolved widths between %s and %s" % (
            self.Lval, self.Jval, self.eLow, self.eHigh)

    def __eq__(self, other):
        return (self.Lval == other.Lval and self.Jval == other.Jval and self.eLow == other.eLow
                and self.eHigh == other.eHigh)


class URRunphysicalLevelSpacing(Warning):

    def __init__(self, Lval, Jval, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.Lval = Lval
        self.Jval = Jval

    def __str__(self):
        return "L=%d J=%.1f level spacing goes outside reasonable limits" % (self.Lval, self.Jval)

    def __eq__(self, other):
        return self.Lval == other.Lval and self.Jval == other.Jval


class URRunphysicalWidth(Warning):

    def __init__(self, Lval, Jval, resonanceReaction, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.Lval = Lval
        self.Jval = Jval
        self.resonanceReaction = resonanceReaction

    def __str__(self):
        return "L=%d J=%.1f average %s widths go outside reasonable limits" % (
            self.Lval, self.Jval, self.resonanceReaction)

    def __eq__(self, other):
        return self.Lval == other.Lval and self.Jval == other.Jval and self.resonanceReaction == other.resonanceReaction


class UnresolvedLink(Warning):

    def __init__(self, link):
        Warning.__init__(self, Level.Fatal, link)
        self.link = link

    def __str__(self):
        return "Unresolved link to %s" % str(self.link)


# reaction objects:

class WicksLimitError(Warning):

    def __init__(self, percentErr, energy_in, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.percentErr = percentErr
        self.energy_in = energy_in

    def __str__(self):
        return "Wick's limit too low by %.3f%% at %s" % (100 * self.percentErr, self.energy_in)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.percentErr == other.percentErr and self.energy_in == other.energy_in


class ZAbalanceWarning(Warning):

    def __init__(self, obj=None):
        Warning.__init__(self, Level.Fatal, obj)

    def __str__(self):
        return "ZA doesn't balance for this reaction!"


class Q_mismatch(Warning):

    def __init__(self, Qcalc, Qactual, obj=None):
        if Qactual.value != 0:
            ratio = Qcalc / Qactual
            if ratio > 1: ratio = 1/ratio
            level = Level.Pedantic
            if ratio < 0.999: level = Level.Minor
            if ratio < 0.99: level = Level.Moderate
            if ratio < 0.95: level = Level.Severe
        else:
            level = Level.Severe  # 0 vs non-zero
        Warning.__init__(self, level, obj)
        self.Qcalc = Qcalc
        self.Qactual = Qactual

    def __str__(self):
        return "Calculated and tabulated Q-values disagree: %s vs %s!" % (self.Qcalc, self.Qactual)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.Qcalc == other.Qcalc and self.Qactual == other.Qactual


class BadFissionEnergyRelease(Warning):

    def __init__(self, worstCase, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.worstCase = worstCase

    def __str__(self):
        return "Fission energy release seems unphysical. Worst case: %s" % self.worstCase

    def __eq__(self, other):
        return self.xpath == other.xpath and self.worstCase == other.worstCase


class NegativeDelayedRate(Warning):

    def __init__(self, rate, obj=None):
        Warning.__init__(self, Level.Fatal, obj)
        self.rate = rate

    def __str__(self):
        return "Negative delayed neutron decay constant: %s" % self.rate

    def __eq__(self, other):
        return self.xpath == other.xpath and self.rate == other.rate


class Threshold_mismatch(Warning):

    def __init__(self, threshold, thresholdCalc, obj=None):
        ratio = threshold / thresholdCalc
        if ratio > 1: ratio = 1/ratio
        level = Level.Pedantic
        if ratio < 0.999: level = Level.Minor
        if ratio < 0.99: level = Level.Moderate
        if ratio < 0.95: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.threshold = threshold
        self.thresholdCalc = thresholdCalc

    def __str__(self):
        return "Calculated and tabulated thresholds disagree: %s vs %s!" % (self.thresholdCalc, self.threshold)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.threshold == other.threshold
                and self.thresholdCalc == other.thresholdCalc)


class Coulomb_threshold_mismatch(Threshold_mismatch):

    def __init__(self, threshold, thresholdCalc, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.threshold = threshold
        self.thresholdCalc = thresholdCalc

    def __str__(self):
        return "Tabulated threshold %s below calculated threshold %s!" % (self.threshold, self.thresholdCalc)


class NonZero_crossSection_at_threshold(Warning):

    def __init__(self, value, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.value = value

    def __str__(self):
        return "First cross section point for threshold reaction should be 0, not %s" % self.value

    def __eq__(self, other):
        return self.xpath == other.xpath and self.value == other.value


class GapInCrossSection(Warning):

    def __init__(self, minGapEnergy, maxGapEnergy, obj=None):
        Warning.__init__(self, Level.Fatal, obj)
        self.minGapEnergy = minGapEnergy
        self.maxGapEnergy = maxGapEnergy

    def __str__(self):
        return "Gap in cross section data from %s to %s" % (self.minGapEnergy, self.maxGapEnergy)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.minGapEnergy == other.minGapEnergy
                and self.maxGapEnergy == other.maxGapEnergy)


class CrossSectionFlatInterpolation(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Moderate, obj)

    def __str__(self):
        return "Cross section uses unphysical flat interpolation"

class NegativeCrossSection(Warning):
    def __init__(self, energy_in, index, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.energy_in = energy_in
        self.index = index

    def __str__(self):
        return "Negative cross section encountered at %s (index %i)!" % (self.energy_in, self.index)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index)


class BadCrossSectionReference(Warning):

    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Cross section reference outside of current reactionSuite not allowed!"


class SummedCrossSectionMismatch(Warning):

    def __init__(self, maxPercentDiff, obj=None):
        ratio = maxPercentDiff / 100
        if ratio > 1: ratio = 1/ratio
        level = Level.Pedantic
        if ratio < 0.999: level = Level.Minor
        if ratio < 0.99: level = Level.Moderate
        if ratio < 0.95: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.maxPercentDiff = maxPercentDiff

    def __str__(self):
        return "Cross section does not match sum of linked reaction cross sections! Max diff: %.2f%%" % (
               self.maxPercentDiff)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.maxPercentDiff == other.maxPercentDiff


class SummedCrossSectionDomainMismatch(Warning):

    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Cross section domain does not match domain of linked reaction cross section sum"


class SummedCrossSectionZeroDivision(Warning):

    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Divide by 0 error when comparing summed cross section to components!"


class SummedMultiplicityMismatch(Warning):

    def __init__(self, maxPercentDiff, obj=None):
        ratio = maxPercentDiff / 100
        if ratio > 1: ratio = 1/ratio
        level = Level.Pedantic
        if ratio < 0.999: level = Level.Minor
        if ratio < 0.99: level = Level.Moderate
        if ratio < 0.95: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.maxPercentDiff = maxPercentDiff

    def __str__(self):
        return "Multiplicity does not match sum of linked product multiplicities! Max diff: %.2f%%" % (
               self.maxPercentDiff)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.maxPercentDiff == other.maxPercentDiff


class SummedMultiplicityDomainMismatch(Warning):

    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Multiplicity domain does not match domain of linked product multiplicity sum"


class SummedMultiplicityZeroDivision(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Divide by 0 error when comparing summed multiplicity to components!"


class NegativeMultiplicity(Warning):
    def __init__(self, value, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.value = value

    def __str__(self):
        return "Encountered negative multiplicity (%s)!" % self.value

    def __eq__(self, other):
        return self.xpath == other.xpath and self.value == other.value


class NonConstantMultiplicity(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Multiplicity should be constant but is energy-dependant!"

    def __eq__(self, other):
        return self.xpath == other.xpath


class Domain_mismatch(Warning):
    def __init__(self, lowBound, highBound, xscLowBound, xscHighBound, obj=None):
        level = Level.Severe
        if lowBound <= xscLowBound and highBound >= xscHighBound:
            level = Level.Moderate
        Warning.__init__(self, level, obj)
        self.lowBound, self.highBound = lowBound, highBound
        self.xscLowBound, self.xscHighBound = xscLowBound, xscHighBound

    def __str__(self):
        return "Domain doesn't match the cross section domain: (%s -> %s) vs (%s -> %s)" % (
            self.lowBound, self.highBound, self.xscLowBound, self.xscHighBound)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.lowBound == other.lowBound
                and self.highBound == other.highBound)


class MissingDistribution(Warning):
    def __init__(self, productName, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.productName = productName

    def __str__(self):
        return "Missing distribution (required for all '%s' products)!" % self.productName

    def __eq__(self, other):
        return self.xpath == other.xpath and self.productName == other.productName


""" eliminate this one?
class NoDistributions(Warning):
    # more serious than missingDistributions:
    def __str__(self):
        return "No distribution data available for *any* of this reaction's products!"
"""


class MissingRecoilDistribution(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Recoil distribution type specified, but recoil partner has unsupported distribution type!"


class WrongDistributionComponent(Warning):
    def __init__(self, component, reactionType='2-body', obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.component = component
        self.reactionType = reactionType

    def __str__(self):
        return "%s is not a valid distribution component for %s reaction" % (self.component, self.reactionType)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.component == other.component
                and self.reactionType == other.reactionType)


class Wrong2BodyFrame(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "2-body reaction not in center-of-mass frame!"


class UncorrelatedFramesMismatch(Warning):
    def __init__(self, angleFrame, energyFrame, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.angleFrame = angleFrame
        self.energyFrame = energyFrame

    def __str__(self):
        return ("For uncorrelated energy/angle distributions, frame must be identical! "
                "Currently, angle given in %s, but energy in %s" % (self.angleFrame, self.energyFrame))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.angleFrame == other.angleFrame
                and self.energyFrame == other.energyFrame)


class FlatIncidentEnergyInterpolation(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "For distributions, flat interpolation along incident energy is unphysical!"


class MissingInterpolationQualifier(Warning):
    def __init__(self, transportable=True, obj=None):
        if transportable:
            level = Level.Moderate
        else:
            level = Level.Pedantic
        Warning.__init__(self, level, obj)

    def __str__(self):
        return "Missing interpolationQualifier for outgoing energy spectrum (should be 'unitBase', 'correspondingPoints', etc.)"


class EnergyDistributionBadU(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "For energy distribution functional form, parameter 'U' results in negative outgoing energies!"


class NegativeDiscreteGammaEnergy(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Discrete gamma energy <= 0"


class PrimaryGammaEnergyTooLarge(Warning):
    """Primary gamma energy should be <= available energy (depending on which discrete level it ends up in)"""

    def __init__(self, energy, fraction, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.energy = energy
        self.fraction = fraction

    def __str__(self):
        return ("Primary gamma energy %s exceeds available energy by %s%%"
                % (self.energy, self.fraction))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy == other.energy
                and self.fraction == other.fraction)


class UnphysicalDiscreteOrPrimaryPhotonMultiplicity(Warning):
    """Max multiplicity for discrete or primary photons should be 1"""

    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return ("Multiplicity > 1 for primary or discrete photon")


class MadlandNixBadParameters(Warning):
    def __init__(self, EFL, EFH, minTm, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.EFL = EFL
        self.EFH = EFH
        self.minTm = minTm

    def __str__(self):
        return ("Madland-Nix fission spectrum contains negative parameters: EFL=%s, EFH=%s, min(Tm)=%s"
                % (self.EFL, self.EFH, self.minTm))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.EFL == other.EFL
                and self.EFH == other.EFH and self.minTm == other.minTm)


class WeightsDontSumToOne(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Weights don't sum to 1.0!"


class UnnormalizedDistribution(Warning):
    def __init__(self, energy_in, index, integral, obj=None):
        ratio = integral if integral < 1 else 1/integral
        level = Level.Pedantic
        if ratio < 0.99: level = Level.Minor
        if ratio < 0.9: level = Level.Moderate
        if ratio < 0.5: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.energy_in = energy_in
        self.index = index
        self.integral = integral

    def __str__(self):
        return "Unnormalized distribution! At energy_in = %s (index %i), integral = %s" % (
            self.energy_in, self.index, self.integral)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index and self.integral == other.integral)


class UnnormalizedDistributionAtMu(Warning):
    """
    Only appears when checking legacy ENDL data (I=1 and 3)
    """

    def __init__(self, mu, integral, obj=None):
        ratio = integral if integral < 1 else 1/integral
        level = Level.Pedantic
        if ratio < 0.99: level = Level.Minor
        if ratio < 0.9: level = Level.Moderate
        if ratio < 0.5: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.mu = mu
        self.integral = integral

    def __str__(self):
        return "Unnormalized distribution! At mu = %s, integral = %s" % (self.mu, self.integral)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.mu == other.mu
                and self.integral == other.integral)


class UnnormalizedKMDistribution(UnnormalizedDistribution):
    pass  # identical to UnnormalizedDistribution, except it needs a special fix() function


class KalbachMannDomainMismatch(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Kalbach-Mann terms do not span the same incident energy domain!"


class IncompleteDistribution(Warning):
    def __init__(self, energy_in, lowerMu, upperMu, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.energy_in = energy_in
        self.lowerMu = lowerMu
        self.upperMu = upperMu

    def __str__(self):
        return ("Incomplete angular coverage (mu = %.2f to %.2f) for distribution at incident energy %s"
                % (self.lowerMu, self.upperMu, self.energy_in))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.lowerMu == other.lowerMu and self.upperMu == other.upperMu)


class NegativeProbability(Warning):

    def __init__(self, value, energy_in, energy_out=None, mu=None, obj=None):
        level = Level.Minor
        if value < -0.05: level = Level.Moderate
        if value < -0.1: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.energy_in = energy_in
        self.energy_out = energy_out
        self.mu = mu
        self.value = value

    def __str__(self):
        msg = "Negative probabilities encountered. Incident energy: %s" % self.energy_in
        if self.mu is not None: msg += ", mu: %s" % self.mu
        if self.energy_out is not None: msg += ", outgoing energy: %s" % self.energy_out
        if self.value is not None: msg += ", worst case: %s" % self.value
        if self.obj is not None: msg += ', %s' % self.obj.moniker
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.energy_out == other.energy_out and self.mu == other.mu)


class ExtraOutgoingEnergy(Warning):
    """
    If an outgoing energy distribution ends with more than one energy with probability=0,
    proper unitbase/correspondingPoints treatment is unclear. Distribution should end with exactly one P=0 point.
    """

    def __init__(self, energy_in, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.energy_in = energy_in

    def __str__(self):
        msg = "Extra zero-probability outgoing energies found at incident energy %s" % self.energy_in
        return msg

    def __eq__(self, other):
        return self.xpath == other.xpath and self.energy_in == other.energy_in


class MissingCoulombIdenticalParticlesFlag(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Severe, obj)

    def __str__(self):
        return "Need 'identicalParticles=\"true\"' when target==projectile"


class IncorrectCoulombIdenticalParticlesFlag(Warning):
    def __init__(self, projectile, target, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.projectile = projectile
        self.target = target

    def __str__(self):
        return "doubleDifferentialCrossSection claims that target (%s) & projectile (%s) are identical!" \
            % (self.target, self.projectile)

    def __eq__(self, other):
        return self.target == other.target and self.projectile == other.projectile


class EnergyImbalance(Warning):
    def __init__(self, energy_in, index, availableEnergy, deposition_per_product, obj=None):
        total_deposited = sum([val[1] for val in deposition_per_product])
        ratio = total_deposited / 100
        # note: using stricter threshold if we're exceeding available energy
        level = Level.Pedantic
        if 0.95 < ratio < 1.01: level = Level.Minor
        elif 0.9 < ratio < 1.05: level = Level.Moderate
        elif ratio < 0.9 or ratio > 1.05: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.energy_in = energy_in
        self.index = index
        self.availableEnergy = availableEnergy
        self.deposition_per_product = deposition_per_product
        self.total_deposited = total_deposited
        if self.availableEnergy == 0: self.total_deposited = float('inf')

    def __str__(self):
        non_zero_products = [(key, val) for key, val in self.deposition_per_product if val]
        per_product = ', '.join(["%s = %.4g%%" % (key, val) for key, val in non_zero_products[:5]])
        if len(non_zero_products) > 5: per_product += ', ...'
        return ("Energy imbalance at incident energy %s (index %i). Total deposited = %.4g%% (%s)" %
                (self.energy_in, self.index, self.total_deposited, per_product))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index and self.availableEnergy == other.availableEnergy
                and self.deposition_per_product == other.deposition_per_product
                and self.total_deposited == other.total_deposited)


class FissionEnergyImbalance(EnergyImbalance):
    def __str__(self):
        per_product = ', '.join(["%s = %.4g%%" % (key, val) for key, val in self.deposition_per_product[:5]])
        if len(self.deposition_per_product) > 5: per_product += ', ...'
        return ("Fission energy imbalance at incident energy %s (index %i). Total deposited = %.4g%% (%s), "
                "leaving insufficient energy for fission products!" %
                (self.energy_in, self.index, self.total_deposited, per_product))


class ValueOutOfRange(Warning):
    def __init__(self, contextMessage, value, lowerBound, upperBound, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.contextMessage = contextMessage
        self.value = value
        self.lowerBound = lowerBound
        self.upperBound = upperBound

    def __str__(self):
        return ("%s. Value=%s, should be in range %s -> %s" % (self.contextMessage, self.value,
                                                               self.lowerBound, self.upperBound))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.contextMessage == other.contextMessage
                and self.value == other.value and self.lowerBound == other.lowerBound
                and self.upperBound == other.upperBound)


class TestSkipped(Warning):
    """ indicate if test was skipped due to missing information """

    def __init__(self, testName, reason, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.testName = testName
        self.reason = reason

    def __str__(self):
        return "Skipped test %s: %s" % (self.testName, self.reason)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.testName == other.testName and self.reason == other.reason


class ExceptionRaised(Warning):
    """If we run into an exception when running check(), try to exit gracefully and return this warning."""

    def __init__(self, Exception_String, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.Exception_String = Exception_String

    def __str__(self):
        return "Encountered Exception: %s" % self.Exception_String

    def __eq__(self, other):
        return self.xpath == other.xpath and self.Exception_String == other.Exception_String


class EnergyDepositionExceptionRaised(ExceptionRaised):
    def __str__(self):
        return "Exception raised when calculating energy deposition: %s" % self.Exception_String


class SkippedCoulombElasticEnergyDeposition(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, Level.Pedantic, obj)

    def __str__(self):
        return "Energy/momentum deposition cannot be computed for charged-particle elastic"


# <editor-fold desc="covarianceSuite warnings">
class CyclicDependency(Warning):
    def __init__(self, cycle, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.cycle = tuple(cycle)

    def __str__(self):
        if len(self.cycle) == 2:
            return "Cyclic dependency in summed covariances for sections %s and %s" % self.cycle
        else:
            return "Cyclic dependency in summed covariances for section %s" % self.cycle

    def __eq__(self, other):
        return self.xpath == other.xpath and self.cycle == other.cycle


class NonMonotonicGrid(Warning):
    def __init__(self, label, obj=None):
        Warning.__init__(self, Level.Fatal, obj)
        self.label = label  # 'row' or 'column'

    def __str__(self):
        return f"Covariance matrix {self.label} grid is not monotonic"

    def __eq__(self, other):
        return self.xpath == other.xpath and self.label == other.label


class MatrixDimensionMismatch(Warning):
    def __init__(self, matrixShape, gridShape, obj=None):
        Warning.__init__(self, Level.Fatal, obj)
        self.matrixShape = matrixShape
        self.gridShape = gridShape

    def __str__(self):
        return f"Matrix shape {self.matrixShape} doesn't match grid shape {self.gridShape}"

    def __eq__(self, other):
        return (self.xpath == other.xpath and
                self.matrixShape == other.matrixShape and
                self.gridShape == other.gridShape)


class NegativeVariance(Warning):
    def __init__(self, index, value, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.index = index
        self.value = value

    def __str__(self):
        return "Encountered negative variance %g at index %i!" % (self.value, self.index)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.index == other.index
                and self.value == other.value)


class VarianceTooSmall(Warning):
    def __init__(self, index, relativeUncertainty, absoluteUncertainty, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.index = index
        self.relativeUncertainty = relativeUncertainty
        self.absoluteUncertainty = absoluteUncertainty

    def __str__(self):
        absoluteString = ""
        if self.absoluteUncertainty:
            absoluteString = ", absolute = %g" % self.absoluteUncertainty
        return "Encountered very small variance at index %i. Relative uncertainty = %g%%%s" % (
            self.index, 100 * self.relativeUncertainty, absoluteString)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.index == other.index
                and self.relativeUncertainty == other.relativeUncertainty
                and self.absoluteUncertainty == other.absoluteUncertainty)


class VarianceTooLarge(Warning):
    def __init__(self, index, relativeUncertainty, absoluteUncertainty, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.index = index
        self.relativeUncertainty = relativeUncertainty
        self.absoluteUncertainty = absoluteUncertainty

    def __str__(self):
        absoluteString = ""
        if self.absoluteUncertainty: absoluteString = ", absolute = %g" % self.absoluteUncertainty
        return "Encountered very large variance at index %i. Relative uncertainty = %g%%%s" % (
            self.index, 100 * self.relativeUncertainty, absoluteString)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.index == other.index
                and self.relativeUncertainty == other.relativeUncertainty
                and self.absoluteUncertainty == other.absoluteUncertainty)


class CorrelationsOutOfRange(Warning):
    def __init__(self, badCount, worstCase, obj=None):
        Warning.__init__(self, Level.Moderate, obj)
        self.badCount = badCount
        self.worstCase = worstCase

    def __str__(self):
        return "%d correlations out of range [-1, 1]. Worst case = %g" % (self.badCount, self.worstCase)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.badCount == other.badCount
                and self.worstCase == other.worstCase)


class NegativeEigenvalues(Warning):
    def __init__(self, negativeCount, worstCase, obj=None):
        level = Level.Pedantic
        if worstCase < -0.001: level = Level.Minor
        if worstCase < -0.1: level = Level.Moderate
        if worstCase < -0.5: level = Level.Severe
        Warning.__init__(self, level, obj)
        self.negativeCount = negativeCount
        self.worstCase = worstCase

    def __str__(self):
        return "%i negative eigenvalues! Worst case = %e" % (self.negativeCount, self.worstCase)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.negativeCount == other.negativeCount
                and self.worstCase == other.worstCase)


class BadEigenvalueRatio(Warning):
    def __init__(self, ratio, obj=None):
        Warning.__init__(self, Level.Minor, obj)
        self.ratio = ratio

    def __str__(self):
        return "Ratio of smallest/largest eigenvalue (%e) is too small" % (self.ratio)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.ratio == other.ratio


class InvalidShortRangeVarianceData(Warning):
    def __init__(self, matrixType, obj=None):
        Warning.__init__(self, Level.Severe, obj)
        self.matrixType = matrixType

    def __str__(self):
        return "shortRangeSelfScalingVariance should contain a diagonal array, contains %s" % self.matrixType

    def __eq__(self, other):
        return self.xpath == other.xpath and self.matrixType == other.matrixType


class ParameterCovarianceMismatch(Warning):
    def __init__(self, nParams, matrixShape, obj=None):
        Warning.__init__(self, Level.Fatal, obj)
        self.nParams = nParams
        self.matrixShape = matrixShape

    def __str__(self):
        return "Number of linked parameters (%d) does not match covariance matrix dimension (%s)" % (
            self.nParams, self.matrixShape)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.nParams == other.nParams and self.matrixShape == other.matrixShape
# </editor-fold>
