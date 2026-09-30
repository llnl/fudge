# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Store and report warnings and errors in a PoPs database.
PoPs.check() returns a nested list of warning objects:

    >>> warnings = PoPs.check()
    >>> print( warnings )

May include or exclude specific classes of warning using the filter command.
filter() returns a new context instance:

    >>> warnings2 = warnings.filter( exclude=[warning.unnormalizedGammas] )

Or, for easier searching you may wish to flatten the list (to get warnings alone without context messages):

    >>> flat = warnings.flatten()
"""

# FIXME context class and base warning class are both identical to stuff in fudge.warning. Move to external utility?


class Context:
    """
    Store warnings in context. This class contains location information (reactionSuite, reaction, etc)
    plus a nested list of warnings or other context instances
    """

    def __init__(self, message="", warningList=None):
        self.message = message
        self.warningList = warningList or []

    def __len__(self):
        return len(self.warningList)

    def __getitem__(self, idx):
        return self.warningList[idx]

    def __str__(self):
        if len(self.warningList) == 0:
            return self.message + ": no problems encountered"
        return "\n".join(self.toStringList())

    def __eq__(self, other):
        return self.message == other.message and self.warningList == other.warningList

    def filter(self, threshold=None, include=None, exclude=None):
        """
        Filter warning list to only include (or exclude) specific classes of warning. For example:

        >>> newWarnings = warnings.filter(exclude=[DiscreteLevelsOutOfOrder])

        Note that if both 'include' and 'exclude' lists are provided, exclude is ignored.
        """

        newWarningList = []
        screened = {}
        if include is None and exclude is None and threshold is None:
            return self, screened
        for warning in self.warningList:
            if not isinstance(warning, Warning):
                newContext, newScreened = warning.filter(include, exclude)
                if newContext:
                    newWarningList.append(newContext)
                for key in newScreened:
                    screened[key] = screened.get(key, 0) + newScreened[key]
            elif threshold is not None:
                pass    # FIXME implement thresholds for PoPs warnings!
            elif include is not None:
                if warning.__class__ in include:
                    newWarningList.append(warning)
            else:  # exclude is not None:
                if warning.__class__ not in exclude:
                    newWarningList.append(warning)
        return Context(self.message, newWarningList), screened

    def flatten(self):
        """
        From a nested hierarchy of warnings, get back a flat list for easier searching:

        >>> w = PoPs.check()
        >>> warningList = w.flatten()

        :return: list containing all of warnings
        """

        List = []
        for val in self.warningList:
            if isinstance(val, Warning):
                List.append(val)
            else:
                List += val.flatten()
        return List

    def toStringList(self, indent="", dIndent="    "):
        """Format warnings for printing. Returns a list of warning strings with indentation."""
        s = ["%s%s" % (indent, self.message)]
        for warning in self.warningList:
            s += warning.toStringList(indent + dIndent)
        return s


class Warning:  # FIXME make abstract base class?
    """
    General warning class. Contains link to problem object,
    xpath in case the object leaves memory,
    and information about the warning or error.
    """

    def __init__(self, obj=None):
        self.obj = obj
        self.xpath = ""
        if hasattr(obj, "toXLink"):
            self.xpath = obj.toXLink()

    def __str__(self):
        return "Generic warning for %s" % self.xpath

    def __eq__(self, other):
        return self.xpath == other.xpath

    def toStringList(self, indent=""):
        return ["%sWARNING: %s" % (indent, self)]


#
# specific warning classes:
#


class NotImplemented(Warning):
    def __init__(self, form, obj=None):
        Warning.__init__(self, obj)
        self.form = form

    def __str__(self):
        return "Checking not yet implemented for %s type data" % self.form

    def __eq__(self, other):
        return self.form == other.form and self.xpath == other.xpath


class UnknownEnergy(Warning):
    def __init__(self, nucleus):
        Warning.__init__(self)
        self.nucleus = nucleus

    def __str__(self):
        return (
            "Could not determine excitation energy for nucleus '%s'" % self.nucleus.id
        )


class DiscreteLevelsOutOfOrder(Warning):
    def __init__(self, lidx, obj=None):
        Warning.__init__(self, obj)
        self.lidx = lidx

    def __str__(self):
        return "Discrete level %s is out of order" % self.lidx

    def __eq__(self, other):
        return self.lidx == other.lidx


class UnnormalizedDecayProbabilities(Warning):
    def __init__(self, branchingSum, obj=None):
        Warning.__init__(self, obj)

        from fudge import warning as fudgeWarningModule

        self.branchingSum = branchingSum
        difference = abs(branchingSum - 1)
        if difference < 0.005:
            self.level = fudgeWarningModule.Level.Minor
        elif difference < 0.05:
            self.level = fudgeWarningModule.Level.Moderate
        else:
            self.level = fudgeWarningModule.Level.Severe

    def __str__(self):
        return "Sum of decay probabilities = %s, should be 1.0!" % (self.branchingSum)

    def __eq__(self, other):
        return self.xpath == other.xpath and self.branchingSum == other.branchingSum


class AliasToNonExistentParticle(Warning):
    def __init__(self, id, pid, obj=None):
        Warning.__init__(self, obj)
        self.id = id
        self.pid = pid

    def __str__(self):
        return "Alias '%s' points to non-existent particle '%s'" % (self.id, self.pid)

    def __eq__(self, other):
        return self.id == other.id and self.pid == other.pid
