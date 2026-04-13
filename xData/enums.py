# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains enum classes needed by xData.

This module contains the following classes:

    +---------------------------+-------------------------------------------------------------------------------+
    | Class                     | Description                                                                   |
    +===========================+===============================================================================+
    | StorageOrder              | This enum class represents the way data can be stored in a GNDS table.        |
    +---------------------------+-------------------------------------------------------------------------------+
    | ValueType                 | This enum class represents the allowed value types for a GNDS values node.    |
    +---------------------------+-------------------------------------------------------------------------------+
    | Interpolation             | This enum class represents the allowed interpolation types for functions.     |
    +---------------------------+-------------------------------------------------------------------------------+
    | InterpolationQualifier    | This enum specifies the interpolation qualifier values for functions.         |
    +---------------------------+-------------------------------------------------------------------------------+
    | Frame                     | This enum specifies the allowed frames for particle's data.                   |
    +---------------------------+-------------------------------------------------------------------------------+
    | FixDomain                 | This enum specifies if the lower, upper or both domain limits are to be       |
    |                           | adjusted to make data for a reaction consistence.                             |
    +---------------------------+-------------------------------------------------------------------------------+
    | Extrapolation             | This enum specifies supported extrapolations for some functions.              |
    +---------------------------+-------------------------------------------------------------------------------+
    | GridStyle                 | This enum specifies allowed styles for the Grid class.                        |
    +---------------------------+-------------------------------------------------------------------------------+
"""

from LUPY import enums as enumsModule

class StorageOrder(enumsModule.Enum):
    """
    This call is an enum class representing the way data can be stored in a table or an array. The two values
    for StorageOrder is 'row-major' or 'column-major' and are represented by the names rowMajor and columnMajor.
    """

    rowMajor = 'row-major'
    columnMajor = 'column-major'

class ValueType(enumsModule.Enum):
    """
    This enum class represents the allowed value types for a GNDS values node.
    """

    integer32 = 'Integer32'
    float32 = 'Float32'
    float64 = 'Float64'

class Interpolation(enumsModule.Enum):
    """
    This enum specifies the interpolation values for functions.
    """

    linlin = 'lin-lin'
    linlog = 'lin-log'
    loglin = 'log-lin'
    loglog = 'log-log'
    flat = 'flat'
    chargedParticle = 'charged-particle'

class InterpolationQualifier(enumsModule.Enum):
    """
    This enum specifies the interpolation qualifier values for functions.
    """
    # We also need direct, correspondingEnergies and correspondingEnergies-unscaled.

    none = ''
    unitBase = 'unitbase'
    unitBaseUnscaled = 'unitbase-unscaled'
    correspondingPoints = 'correspondingPoints'
    correspondingPointsUnscaled = 'correspondingPoints-unscaled'
    cumulativePoints = 'cumulativePoints'
    cumulativePointsUnscaled = 'cumulativePointsUnscaled'

class Frame(enumsModule.Enum):
    """
    This enum specifies the allowed frames for particle's data (e.g., a paricles distribution)..
    """

    none = enumsModule.auto()
    lab = enumsModule.auto()
    centerOfMass = enumsModule.auto()
    product = enumsModule.auto()

class FixDomain(enumsModule.Enum):
    """
    This enum specifies if the lower, upper or both domain limits are to be adjusted to make data for a reaction consistence.
    """

    upper = enumsModule.auto()
    lower = enumsModule.auto()
    both = enumsModule.auto()

class Extrapolation(enumsModule.Enum):
    """
    This enum specifies supported extrapolations for some functions.
    """

    none = enumsModule.auto()
    flat = enumsModule.auto()

class GridStyle(enumsModule.Enum):
    """
    This enum specifies allowed styles for the Grid class.
    """

    none = enumsModule.auto()
    points = enumsModule.auto()
    boundaries= enumsModule.auto()
    parameters= enumsModule.auto()

class BorderAddition(enumsModule.Enum):
    '''
    This enum specifies how to treat the y-value at the border between two 1d curves is to be added.
    Description of options:

    lower - use the y-value from the lower curve.
    upper - use the y-value from the upper curve.
    average - use the average of the lower and upper curve's y-values.
    sum - use the sum of the lower and upper curve's y-values.
    '''

    average = enumsModule.auto()
    lower = enumsModule.auto()
    upper = enumsModule.auto()
    sum = enumsModule.auto()
