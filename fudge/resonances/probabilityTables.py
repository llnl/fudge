# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Container for URR probability tables.

GNDS-2.0 doesn't have an official probability table container,
so these currently reside in the <applicationData> section.

Example:
<institution label="LLNL::URR_probability_tables">
  <probabilityTables>
    <probabilityTable label="">  (one per temperature)
      <incidentEnergy value="" unit="" normalized="true">
        <table>
          <columnHeaders>...</columnHeaders>
          <data>...</data></table></incidentEnergy>
        ...
"""

from fudge import suites as suitesModule
from LUPY import ancestry as ancestryModule
from xData import table as tableModule

LLNLProbabilityTablesToken = "LLNL::URR_probability_tables"

class IncidentEnergy(ancestryModule.Ancestry):

    moniker = 'incidentEnergy'
    keyName = 'value'

    def __init__(self, value, unit, table_ = None, normalized=True):

        self.value = value
        self.unit = unit
        self.normalized = normalized   # whether table entries have been divided by mean cross section
        self.table = table_

    @property
    def table(self):
        return self.__table

    @table.setter
    def table(self, value):
        if not isinstance(value, tableModule.Table):
            raise TypeError(f"Expected Table instance, got {type(value)}")
        self.__table = value
        self.__table.setAncestor(self)

    def convertUnits(self, unitMap):
        from pqu import PQU as PQUModule
        if self.unit in unitMap:
            newUnit = unitMap[self.unit]
            factor = PQUModule.PQU(1, self.unit).getValueAs(newUnit)
            self.value *= factor
            self.unit = newUnit
        self.table.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + '  '
        attrs = ""
        if not self.normalized:
            attrs = ' normalized="false"'
        xmlString = ['%s<%s value="%r" unit="%s"%s>' %
                        (indent, self.moniker, self.value, self.unit, attrs)]
        xmlString += self.table.toXML_strList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        value = element.get('value')
        if value is None:
            value = element.get('incidentEnergy')                  # Support legacy for now, but not much longer.
        xPath.append('%s[@value="%r"]' % (cls.moniker, value))

        value = float(value)
        unit = element.get("unit")
        if unit is None:
            unit = element.get("unit")                              # Support legacy for now, but not much longer.
        normalized = bool(element.get("normalized", "true").title())
        table = tableModule.Table.parseNodeUsingClass(
                    element.find(tableModule.Table.moniker), xPath, linkData, **kwargs)

        incidentEnergy = cls(value, unit, table, normalized)

        xPath.pop()
        return incidentEnergy


class ProbabilityTable(ancestryModule.Ancestry):

    moniker = "probabilityTable"
    keyName = "label"

    def __init__(self, label, _incidentEnergies=None):

        self.label = label
        self.__incidentEnergies = _incidentEnergies or []

    def __len__(self):
        return len(self.__incidentEnergies)

    def __getitem__(self, index):
        return self.__incidentEnergies[index]

    def add(self, _incidentEnergy):
        assert isinstance(_incidentEnergy, IncidentEnergy), f"Expected IncidentEnergy instance, got {type(_incidentEnergy)} instead"
        _incidentEnergy.setAncestor(self)
        self.__incidentEnergies.append(_incidentEnergy)

    def convertUnits(self, unitMap):
        """
        Convert all incidentEnergy units in this probability table.
        """

        for incidentEnergy in self:
            incidentEnergy.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + '  '
        xmlString = ['%s<%s label="%s">' % (indent, self.moniker, self.label)]
        for incidentEnergy in self:
            xmlString += incidentEnergy.toXML_strList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append('%s[@label="%s"]' % (cls.moniker, element.get('label')))

        tmp = cls(element.get('label'))
        for child in element:
            tmp.add(IncidentEnergy.parseNodeUsingClass(child, xPath, linkData, **kwargs))
        xPath.pop()
        return tmp


class ProbabilityTables(suitesModule.ExclusiveSuite):

    moniker = "probabilityTables"

    def __init__( self ) :

        suitesModule.ExclusiveSuite.__init__( self, [ ProbabilityTable ] )
