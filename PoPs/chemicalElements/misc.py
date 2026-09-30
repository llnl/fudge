# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Module containing dictionaries and helper functions for handling particle ids and symbols.

Defines the following dictionaries:
    symbolFromZ, nameFromZ          # look up chemical symbol or name by Z (nameFromZ returns default name)
    ZfromSymbol, nameFromSymbol     # look up Z or name by chemical symbol (nameFromSymbol returns default name)
    ZfromName, symbolFromName       # look up Z or symbol by chemical name (accepts all alternative spellings)
    namesFromZ, namesFromSymbol     # look up chemical name by Z or by Symbol (allows alternate names)
"""

import math

from pqu import PQU as PQUModule
from .. import misc as miscModule
from .. import IDs as IDsModule

MaxA = 400
continuumID = 1000000
sumID = continuumID + 1

chemicalElementZSymbolNames = (
    (   1, "H",  ("Hydrogen",) ),           (   2, "He", ("Helium",) ),             (   3, "Li", ("Lithium",) ),
    (   4, "Be", ("Beryllium",) ),          (   5, "B",  ("Boron",) ),              (   6, "C",  ("Carbon",) ),
    (   7, "N",  ("Nitrogen",) ),           (   8, "O",  ("Oxygen",) ),             (   9, "F",  ("Fluorine",) ),
    (  10, "Ne", ("Neon",) ),               (  11, "Na", ("Sodium",) ),             (  12, "Mg", ("Magnesium",) ),
    (  13, "Al", ("Aluminium", "Aluminum") ), (  14, "Si", ("Silicon",) ),          (  15, "P",  ("Phosphorus",) ),
    (  16, "S",  ("Sulfur", "Sulphur") ),   (  17, "Cl", ("Chlorine",) ),           (  18, "Ar", ("Argon",) ),
    (  19, "K",  ("Potassium",) ),          (  20, "Ca", ("Calcium",) ),            (  21, "Sc", ("Scandium",) ),
    (  22, "Ti", ("Titanium",) ),           (  23, "V",  ("Vanadium",) ),           (  24, "Cr", ("Chromium",) ),
    (  25, "Mn", ("Manganese",) ),          (  26, "Fe", ("Iron",) ),               (  27, "Co", ("Cobalt",) ),
    (  28, "Ni", ("Nickel",) ),             (  29, "Cu", ("Copper",) ),             (  30, "Zn", ("Zinc",) ),
    (  31, "Ga", ("Gallium",) ),            (  32, "Ge", ("Germanium",) ),          (  33, "As", ("Arsenic",) ),
    (  34, "Se", ("Selenium",) ),           (  35, "Br", ("Bromine",) ),            (  36, "Kr", ("Krypton",) ),
    (  37, "Rb", ("Rubidium",) ),           (  38, "Sr", ("Strontium",) ),          (  39, "Y",  ("Yttrium",) ),
    (  40, "Zr", ("Zirconium",) ),          (  41, "Nb", ("Niobium",) ),            (  42, "Mo", ("Molybdenum",) ),
    (  43, "Tc", ("Technetium",) ),         (  44, "Ru", ("Ruthenium",) ),          (  45, "Rh", ("Rhodium",) ),
    (  46, "Pd", ("Palladium",) ),          (  47, "Ag", ("Silver",) ),             (  48, "Cd", ("Cadmium",) ),
    (  49, "In", ("Indium",) ),             (  50, "Sn", ("Tin",) ),                (  51, "Sb", ("Antimony",) ),
    (  52, "Te", ("Tellurium",) ),          (  53, "I",  ("Iodine",) ),             (  54, "Xe", ("Xenon",) ),
    (  55, "Cs", ("Caesium", "Cesium") ),   (  56, "Ba", ("Barium",) ),             (  57, "La", ("Lanthanum",) ),
    (  58, "Ce", ("Cerium",) ),             (  59, "Pr", ("Praseodymium",) ),       (  60, "Nd", ("Neodymium",) ),
    (  61, "Pm", ("Promethium",) ),         (  62, "Sm", ("Samarium",) ),           (  63, "Eu", ("Europium",) ),
    (  64, "Gd", ("Gadolinium",) ),         (  65, "Tb", ("Terbium",) ),            (  66, "Dy", ("Dysprosium",) ),
    (  67, "Ho", ("Holmium",) ),            (  68, "Er", ("Erbium",) ),             (  69, "Tm", ("Thulium",) ),
    (  70, "Yb", ("Ytterbium",) ),          (  71, "Lu", ("Lutetium",) ),           (  72, "Hf", ("Hafnium",) ),
    (  73, "Ta", ("Tantalum",) ),           (  74, "W",  ("Tungsten",) ),           (  75, "Re", ("Rhenium",) ),
    (  76, "Os", ("Osmium",) ),             (  77, "Ir", ("Iridium",) ),            (  78, "Pt", ("Platinum",) ),
    (  79, "Au", ("Gold",) ),               (  80, "Hg", ("Mercury",) ),            (  81, "Tl", ("Thallium",) ),
    (  82, "Pb", ("Lead",) ),               (  83, "Bi", ("Bismuth",) ),            (  84, "Po", ("Polonium",) ),
    (  85, "At", ("Astatine",) ),           (  86, "Rn", ("Radon",) ),              (  87, "Fr", ("Francium",) ),
    (  88, "Ra", ("Radium",) ),             (  89, "Ac", ("Actinium",) ),           (  90, "Th", ("Thorium",) ),
    (  91, "Pa", ("Protactinium",) ),       (  92, "U",  ("Uranium",) ),            (  93, "Np", ("Neptunium",) ),
    (  94, "Pu", ("Plutonium",) ),          (  95, "Am", ("Americium",) ),          (  96, "Cm", ("Curium",) ),
    (  97, "Bk", ("Berkelium",) ),          (  98, "Cf", ("Californium",) ),        (  99, "Es", ("Einsteinium",) ),
    ( 100, "Fm", ("Fermium",) ),            ( 101, "Md", ("Mendelevium",) ),        ( 102, "No", ("Nobelium",) ),
    ( 103, "Lr", ("Lawrencium",) ),         ( 104, "Rf", ("Rutherfordium",) ),      ( 105, "Db", ("Dubnium",) ),
    ( 106, "Sg", ("Seaborgium",) ),         ( 107, "Bh", ("Bohrium",) ),            ( 108, "Hs", ("Hassium",) ),
    ( 109, "Mt", ("Meitnerium",) ),         ( 110, "Ds", ("Darmstadtium",) ),       ( 111, "Rg", ("Roentgenium",) ),
    ( 112, "Cn", ("Copernicium",) ),        ( 113, "Nh", ("Nihonium",) ),           ( 114, "Fl", ("Flerovium",) ),
    ( 115, "Mc", ("Moscovium",) ),          ( 116, "Lv", ("Livermorium",) ),        ( 117, "Ts", ("Tennessine",) ),
    ( 118, "Og", ("Oganesson",) ) )

symbolFromZ = {}
for Z, symbol, names in chemicalElementZSymbolNames: symbolFromZ[Z] = symbol

nameFromZ = {}
for Z, symbol, names in chemicalElementZSymbolNames: nameFromZ[Z] = names[0]

namesFromZ = {}
for Z, symbol, names in chemicalElementZSymbolNames: namesFromZ[Z] = names

ZFromSymbol = {}
for Z, symbol, names in chemicalElementZSymbolNames: ZFromSymbol[symbol] = Z

nameFromSymbol = {}
for Z, symbol, names in chemicalElementZSymbolNames: nameFromSymbol[symbol] = names[0]

namesFromSymbol = {}
for Z, symbol, names in chemicalElementZSymbolNames: namesFromSymbol[symbol] = names

ZFromName = {}
for Z, symbol, names in chemicalElementZSymbolNames:
    for name in names:
        ZFromName[name] = Z

symbolFromName = {}
for Z, symbol, names in chemicalElementZSymbolNames :
    for name in names:
        symbolFromName[name] = symbol


def checkZ(Z):

    if not isinstance(Z, int):
        raise TypeError('Z not an int object.')
    if 0 < Z <= chemicalElementZSymbolNames[-1][0]:
        return Z
    raise ValueError('Z = "%s" out of range.' % Z)


def checkSymbol(symbol):

    if not isinstance(symbol, str):
        raise TypeError('symbol not a str object.')
    if symbol not in ZFromSymbol:
        raise ValueError('Invalid symbol: %s.' % miscModule.toLimitedString(symbol))


def checkA(A):

    if not isinstance(A, int):
        raise TypeError('A not an int')
    if not 0 <= A < MaxA:
        raise ValueError('integer A out of range: A = %s' % A)
    return A


def checkIndex(index):

    if not isinstance(index, int):
        raise TypeError('index attribute not int object.')
    if index < 0:
        raise ValueError('Negative level index = %s' % index)
    return index


def chemicalElementALevelIDsAndAnti(id, qualifierAllowed=False):
    """
    Parse a particle id to extract the following information::

        baseId: particle id with any qualifiers removed
        chemicalElementSymbol: symbol if id is a chemical element, isotope, nuclide or nucleus. None otherwise
        A: total nucleon number (as integer) if id is for an isotope, nuclide or nucleus. None otherwise
        levelId: integer level index if id is a nuclide or nucleus. None otherwise
        isNucleus: True if id is a nucleus, False if id is a chemical element, isotope or nuclide, None otherwise
        anti: '_anti' if id contains that string, empty string otherwise
        qualifier: string qualifier if one is found (qualifiers appear inside {curly brackets}. Empty string otherwise

    Unless qualifierAllowed = True, ValueError will be raised if a qualifier is detected in id

    :param id: string particle id or symbol
    :param qualifierAllowed: boolean, whether qualifiers are allowed in the id
    :return: tuple(baseId, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier)
    """

    baseID, anti, qualifier = miscModule.baseAntiQualifierFromID(id, qualifierAllowed=qualifierAllowed)

    chemicalElementSymbol = None
    A = None
    levelID = None
    isNucleus = False
    if baseID.count('_e') < 2:        # Can have 0 or 1 '_e' (e.g., 'O16', 'O16_e3') if chemicalElement derived.
        isotopeID, sep, levelID = baseID.partition('_e')
        if sep == '_e':
            try:
                levelID = int(levelID)
                checkIndex(levelID)
            except:
                levelID = None
                isotopeID = ''
        else:
            levelID = None

        if len(isotopeID) > 0:
            for i1, character in enumerate(isotopeID):
                if character.isdigit( ): break
            isotopeIDUpper = isotopeID[0].upper( ) + isotopeID[1:]
            if character.isdigit():
                chemicalElementSymbol = isotopeIDUpper[:i1]
                try:
                    A = int(isotopeIDUpper[i1:])
                    checkSymbol(chemicalElementSymbol)
                    checkA(A)
                    isNucleus = isotopeIDUpper != isotopeID
                    if levelID is None: levelID = 0
                except:
                    chemicalElementSymbol = None
            else:
                try:
                    checkSymbol(isotopeID)
                    chemicalElementSymbol = isotopeID
                except:
                    chemicalElementSymbol = None

    if chemicalElementSymbol is None: A, levelID, isNucleus = None, None, None
    return baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier


def ZAInfo(particle):
    """
    Compute Z, A, ZA and level index for a particle instance.

    :param particle: PoPs particle instance. If the particle is not a baryon, isotope, nuclide or nucleus,
        Z, A, ZA and level index will all be 0.
    :return: tuple( Z, A, ZA, levelIndex ).
    """

    from .. import alias as aliasModule
    from .. import database as databaseModule
    from ..families import nuclide as nuclideModule
    from ..families import nucleus as nucleusModule
    from ..families import baryon as baryonModule
    from ..families import unorthodox as unorthodoxModule
    from . import isotope as isotopeModule

    level = 0
    Z = 0
    A = 0
    if isinstance(particle, (isotopeModule.Isotope, nuclideModule.Particle, nucleusModule.Particle)):
        Z = particle.Z
        A = particle.A
        if isinstance(particle, (nuclideModule.Particle, nucleusModule.Particle)):
            level = particle.index
    elif isinstance(particle, aliasModule.BaseAlias):
        return ZAInfo(particle.findClassInAncestry(databaseModule.Database).final(particle.pid))
    elif isinstance(particle, baryonModule.Particle):
        if particle.id == IDsModule.neutron:
            A = 1
        if particle.id == IDsModule.proton:
            Z = 1
            A = 1
    elif isinstance(particle, unorthodoxModule.Particle):
        if particle.id == IDsModule.FissionProductENDL99120:
            Z = 99
            A = 120
        elif particle.id == IDsModule.FissionProductENDL99125:
            Z = 99
            A = 125
        elif len(particle.charge) > 0:
            Z = particle.charge[0].value

    try:
        return Z, A, 1000 * Z + A, level            # FIXME raise if it didn't match anything?
    except:
        raise


def ZAInfo_fromString(particleName):
    """
    Parses an isotope name (e.g., 'Mn55_e3') to get its Z, A, ZA and nuclear level index.
    If the level is negative, the particle is a metastable name (e.g., 'Am242_m1').
    The name can be an isotope which returns ( Z, A, 1000 * Z + A, level ) (e.g., 'O16' returns (8, 16, 8016, 0)), 
    a natural isotope which returns ( Z, 0, 1000 * Z, level) (e.g., 'O0' returns (8, 0, 8000, 0)) 
    or just a chemical element's symbol which returns ( Z, 0, 0, 0 ) (e.g., 'O' returns (8, 0, 0, 0)).
    Also handles PoPs (GNDS) nucleus name which is the nucleus' isotope name but with the first character lower-cased
    (e.g., 'o16' for the nucleus of Oyxgen-16) as well as 'n', 'p', 'd', 't', 'h' and 'a'.
    Handles nuclear excited levels (e.g., 'Am242_e2' returns (95, 242, 95242, 2)) and nuclear metastables 
    (e.g., 'Am242_m1' returns (95, 242, 95242, -1)).

    :param particleName: string
    :return: tuple( Z, A, ZA, levelIndex )
    """

    import re

    if particleName == IDsModule.FissionProductENDL99120:
        return 99, 120, 99120, 0
    elif particleName == IDsModule.FissionProductENDL99125:
        return 99, 125, 99125, 0

    if len(particleName) > 1:
        particleName = particleName.capitalize( )

    try:
        return {IDsModule.neutron: (0, 1, 1, 0),                IDsModule.proton: (1, 1, 1001, 0),
                IDsModule.familiarDeuteron: (1, 2, 1002, 0),    IDsModule.familiarTriton: (1, 3, 1003, 0),
                IDsModule.familiarHelion: (2, 3, 2003, 0),      IDsModule.familiarAlpha: (2, 4, 2004, 0)   }[particleName]
    except:
        pass

    pattern = re.compile( "([A-Za-z]+)([0-9]+)(_[em])?([0-9]+)?" )
    result = pattern.match( particleName )
    if result is None:
        try:
            Z = ZFromSymbol[particleName]                   # Assumes that particleName is just the chemical element part (e.g., 'Xe').
        except:
            Z = 0                                           # Not a valid isotope or chemical element name so return all 0's.
        return Z, 0, 0, 0

    symbol, A, levelIndicator, level = result.groups()

    reconstructed = symbol + A
    if levelIndicator is not None:
        if level is not None:
            reconstructed += levelIndicator + level
    if particleName != reconstructed:
        return 0, 0, 0, 0

    A = int(A)
    if level is None:
        level = 0
    else:
        level = int(level)
        if levelIndicator == '_m':
            level = -level

    Z = ZFromSymbol[symbol]

    return Z, A, 1000 * Z + A, level


def ZA(particle):
    """
    Uses ZAInfo to compute the particle ZA. See ZAInfo for more detail
    """

    return ZAInfo(particle)[2]


def idFromZAndA(Z, A, allowNeutron=True):
    """
    Compute the PoPs id for a particle given Z and A.

    :param Z:               Atomic number of the particle.
    :param A:               Atomic mass number of the particle.
    :param allowNeutron:    If True, Z = 0 and A = 1 will return the neutron id, otherwise a raise is executed.

    :return:                String particle id
    """

    from .. import IDs as IDsModule

    if Z == 0 and A == 1 and allowNeutron:
        return IDsModule.neutron

    return isotopeSymbolFromChemicalElementIDAndA(symbolFromZ[Z], A)


def idFromZA(ZA, allowNeutron=True):
    """
    Computes the PoPs id for a particle given its ZA.

    :param ZA:              1000 * Z + A where Z and N are the particle's atomic number and mass number, respectively.
    :param allowNeutron:    If True, Z = 0 and A = 1 will return the neutron id, otherwise a raise is executed.

    :return:                Particle's PoPs id.
    """

    return idFromZAndA(ZA // 1000, ZA % 1000, allowNeutron)


def nucleusIDFromZAndA(Z, A, allowNeutron=True):
    """
    Equivalent to calling nucleusIDFromNuclideID(idFromZAndA(Z, A, allowNeutron)).

    :param Z:               Atomic number of the particle.
    :param A:               Atomic mass number of the particle.
    :param allowNeutron:    If True, Z = 0 and A = 1 will return the neutron id, otherwise a raise is executed.

    :return:                String particle id
    """

    return nucleusIDFromNuclideID(idFromZAndA(Z, A, allowNeutron))


def hasNucleus(particle, nucleusReturnsTrue = False):
    """
    Checks whether particle contains a nucleus.

    :param particle: PoPs particle instance
    :param nucleusReturnsTrue: boolean. If False, hasNucleus returns False when particle is a nucleus instance
    """

    from . import isotope as isotopeModule
    from ..families import nuclide as nuclideModule
    from ..families import nucleus as nucleusModule

    if isinstance(particle, (isotopeModule.Isotope, nuclideModule.Particle)):
        return True
    return nucleusReturnsTrue and isinstance(particle, nucleusModule.Particle)


def isotopeSymbolFromChemicalElementIDAndA(elementID, A):

    checkSymbol(elementID)
    checkA(A)
    return "%s%s" % (elementID, A)


def nuclideIDFromIsotopeSymbolAndIndex(isotopeSymbol, index):

    checkIndex(index)
    if index == 0:
        return isotopeSymbol
    return "%s_e%s" % (isotopeSymbol, index)


def nuclideIDFromZ_A_AndIndex(Z, A, index):
    """
    This function returns the nuclide id from a nuclide's Z, A and nuclear level index.

    :param Z:       The atomic number of the nuclide.
    :param A:       The atomic mass number of the nuclide.
    :param index:   The nuclear level index of the nuclide.
    """

    return nuclideIDFromIsotopeSymbolAndIndex(idFromZAndA(Z, A), index)


def nucleusIDFromNuclideID(nuclideID):
    """
    Returns the nucleus PoPs id from a nuclide PoPs id by converting 1st letter to lower case.
    No check is performed on the validity of *nuclideID* as a nuclide PoPs id.

    :param nuclideID:       PoPs id for a nuclide.

    :return:                The nucleus PoPs id for the specified nuclide PoPs id.
    """

    return nuclideID[0].lower() + nuclideID[1:]


def nuclideIDFromNucleusID(nucleusID):
    """ The nuclide id is computed from nucleus id by converting 1st letter to upper case """

    return nucleusID[0].upper( ) + nucleusID[1:]


def nuclearBindingEnergyPerNucleonSemiEmpirical(Z, A, unit):
    """
    Returns the nuclear binding energy per nucleon from a semi empirical formula (see https://en.wikipedia.org/wiki/Nuclear_binding_energy).
    """

    N = A - Z
    pairingTerm = 0
    if A % 2 == 0:
        pairingTerm = 33.0
        if Z % 2 == 1:
            pairingTerm *= -1

    energy_MeV = 14.0 - 13.0 * math.pow(A, -1/3) - 0.585 * Z**2 * math.pow(A, -4/3) - 19.3 * (N - Z)**2 / A**2 + pairingTerm * math.pow(A, -7/4)

    return PQUModule.PQU(energy_MeV, 'MeV').getValueAs(unit)
