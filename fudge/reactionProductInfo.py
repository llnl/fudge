
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import copy

from LUPY import enums as enumsModulon
from PoPs import IDs as IDsPoPsModule
from PoPs import specialNuclearParticleID as specialNuclearParticleIDModule
from PoPs import IDs as IDsPoPsModule
from PoPs.chemicalElements import misc as miscModule
from fudge import enums as enumsModule
from fudge.productData import multiplicity as multiplicityModule

__doc__ = """
This module contains the ReactionProductInfo class which is used labels for a reaction.
"""

class ReactionProductInfo:
    """
    This class analyzes a Reaction class instance and stores data to calculate the reactions tag, short tag, or
    ask if the reaction is a elastic, capture, fission or sumOfRemainingOutputChannels reaction.

    A reaction's tag is the familiar formula for the reaction of the format T(P,Ps)R where T is the protare's
    target, P is the protare's projectile, Ps is the list of products excluding the heaviest product and R is 
    the heaviest product. With the exception of capture, photons are not shown it the tag. Each product in Ps 
    is either the product's PoPs id if its multplicity is 1 or its mmultplicity followed by its PoPs id. 
    The list of products in the tag can be either the initial products or the final products. The short tag is
    the string bewteen the ',' and ')' of the tag For fission, the tag and short tag are 'fission'. 
    The following table list some example:

        +-------+-------------------------------------------------------+-------------------------------+-----------+
        | item  | Reaction label                                        | tag                           | short tag |
        +=======+=======================================================|===============================|===========+
        |    1  | n + Li6 --> n + d + He4                               | Li6(n,nd)He4                  | nd        |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    2  | n + Li6 --> 2n + p + He4                              | Li6(n,2np)He4                 | 2np       |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    3  | n + Be9 --> 2n + 2He4                                 | Be9(n,2na)He4                 | 2na       |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    4  | n + Be9 --> C14 + photon [inclusive]                  | C13(n,g)C14                   | g         |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    5  | n + U235 --> m(E)*n + photon [total fission]          | fission                       | fission   |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    6  | n + U235 --> n + (U235_e3 -> U235 + photon)           | U235(n,n)U235_e3              | n         |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    7  | n + U235 --> n + (U235_e3 -> U235 + photon)           | U235(n,n)U235                 | n         |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    8  | n + B10 --> n + (B10_e12 -> d + 2He4)                 | B10(n,n)B10_e12               | n         |
        +-------+-------------------------------------------------------+-------------------------------+-----------+
        |    9  | n + B10 --> n + (B10_e12 -> d + 2He4)                 | B10(n,nda)He4                 | nda       |
        +-------+-------------------------------------------------------+-------------------------------+-----------+

    In the table, all tags are for the final products (i.e., initialProducts=False) except for items 6 and 8.
    """

    def __init__(self, reaction):
        """
        :param reaction:        A Reaction class instance.
        """

        def getMultiplicity(product):
            """
            For internal use.
            This function returns the multiplicity of *product*. If the product's multiplicity is not
            an integer or it is energy dependent, a 0 is returned.

            :param product:     An instance of the Product class.
            """

            multiplicity = product.multiplicity[0]

            if isinstance(multiplicity, multiplicityModule.Constant1d):
                value = multiplicity.value
                if value.is_integer():
                    return int(value)
                return 0
            elif isinstance(multiplicity, (multiplicityModule.XYs1d, multiplicityModule.Regions1d,
                    multiplicityModule.Branching1d, multiplicityModule.Polynomial1d)):
                return 0

            raise TypeError('Need to implement ' + multiplicity.moniker)

        def addProduct(products, product):
            """
            For internal use.
            This function adds/updates *product* to the dictionary *products*. If *product* is not in *products*,
            it GNDS PoPs id is added as a key and the key's value is it multiplicity as returned by the function
            **getMultiplicity**.  If *product* is a key in *products*, its multiplicity is added to the current 
            multiplicity for that product.  If the multiplicity or a prior multiplicity is energy dependent or 
            not an integer, the multiplicity is set to 0.
            """

            pid = product.pid
            multiplicityValue = getMultiplicity(product)
            if pid in products:
                if multiplicityValue == 0:
                    products[pid] = 0
                else:
                    currentMultiplicityValue = products[pid]
                    if currentMultiplicityValue != 0:
                        products[pid] += multiplicityValue
            else:
                products[pid] = multiplicityValue

        def processOutputChannel(outputChannel, initialProducts, finalProducts):
            """
            For internal use.
            This function loops over the product in *outputChannel* and adds them to initialProducts if it is
            not **None**. If a product does not have an output channel, it is added to *finalProducts*; otherwise,
            this function is called with the product's output channel as the first argument.

            :param outputChannel:       The output channel whose products are processed.
            :param initialProducts:     None of the dictionary of initial products.
            :param finalProducts:       The dictionary of final products.
            """

            if outputChannel is None:
                return

            for product in outputChannel.products:
                if initialProducts is not None:
                    addProduct(initialProducts, product)
                if product.outputChannel is None:
                    addProduct(finalProducts, product)
                else:
                    processOutputChannel(product.outputChannel, None, finalProducts)

        def sortProducts(__products, pops):
            """
            This function inserts the products in **__products** into a new dictionary ordered by mass 
            except for the neutron.  If a neutron is present, it will always be the first product inserted 
            into the returned dictionary. Note, since Python 3.7, the standard Python dictionary maintains 
            insertion order when iterating over the entries.

            :param __products:      Dictionary whose products are sorted.
            :param pops:            A PoPs instance for getting a particle's mass.
            """

            products = {}
            pids = list(__products.keys())
            for pid in ['n', 'H1', 'H2', 'H3', 'He3', 'He4', 'photon']:
                if pid in __products:
                    products[pid] = __products.pop(pid)
                    pids.pop(pids.index(pid))

            masses = {}
            for pid in pids:
                try:
                    mass = pops.estimateMass(pid, 'amu')
                except:
                    mass = 99999
                masses[pid] = mass

            while len(__products) > 0:
                pid1 = pids[0]
                mass1 = masses[pid1]
                for pid2 in pids:
                    mass2 = masses[pid2]
                    if mass2 < mass1:
                        pid1 = pid2
                        mass1 = mass2
                products[pid1] = __products.pop(pid1)
                pids.pop(pids.index(pid1))

            return products

        protare = reaction.rootAncestor

        self.projectileId = protare.projectile
        self.targetId = protare.target
        self.interaction = protare.interaction

        initialProducts = {}
        finalProducts = {}

        processOutputChannel(reaction.outputChannel, initialProducts, finalProducts)

        self.__isElastic = False
        self.__isCapture = False
        self.__isFission = reaction.isFission()
        self.__isSumOfRemainingOutputChannels = reaction.label == str(enumsModule.Genre.sumOfRemainingOutputChannels)
        if initialProducts == finalProducts:
            if specialNuclearParticleIDModule.sameSpecialNuclearParticle(self.projectileId, self.targetId):
                if len(initialProducts) == 1:               # Special check for projectile and target being the same particle.
                    pid = list(initialProducts.keys())[0]
                    self.__isElastic = specialNuclearParticleIDModule.sameSpecialNuclearParticle(self.projectileId, pid) \
                            and initialProducts[pid] == 2
            if len(initialProducts) == 2 and not self.__isElastic:
                if self.projectileId in initialProducts and self.targetId in initialProducts:
                    projectileMultiplity = initialProducts[self.projectileId]
                    targetMultiplity = initialProducts[self.targetId]
                    self.__isElastic = projectileMultiplity == targetMultiplity == 1
                if not self.__isElastic:
                    Zp, Ap, ZA, levelIndex = miscModule.ZAInfo_fromString(self.projectileId)
                    Zt, At, ZA, levelIndex = miscModule.ZAInfo_fromString(self.targetId)
                    residual = miscModule.idFromZAndA(Zp + Zt, Ap + At)
                    self.__isCapture = IDsPoPsModule.photon in initialProducts and residual in initialProducts

        pops = protare.PoPs
        self.initialProducts = sortProducts(initialProducts, pops)
        self.finalProducts = sortProducts(finalProducts, pops)

    @property
    def isElastic(self):

        return self.__isElastic

    @property
    def isCapture(self):

        return self.__isCapture

    @property
    def isFission(self):

        return self.__isFission

    @property
    def isSumOfRemainingOutputChannels(self):

        return self.__isSumOfRemainingOutputChannels

    def tag(self, initialProducts, mode=specialNuclearParticleIDModule.Mode.familiar):
        """
        This method returns the tag of the reaction.

        :param initialProducts:     A boolean. If **True**, the initial products tag is returned; otherwise, the final products tag is returned.
        :param mode:                Determines which id is used for the the light particles (e.g., for 'H3', is is 'H3', 'h3' or 't').
        """

        if self.__isFission:
            return 'fission'
        elif self.__isSumOfRemainingOutputChannels:
            return 'sumOfRemainingOutputChannels'

        __tag = '%s(%s,' % (self.targetId, specialNuclearParticleIDModule.specialNuclearParticleID(self.projectileId, mode))

        products = self.initialProducts if initialProducts else self.finalProducts
        products2 = copy.copy(products)

        if IDsPoPsModule.photon in products2:
            products2.pop(IDsPoPsModule.photon)

        pids = [pid for pid in products2]
        if len(pids) == 1:
            others, rid = [], pids[0]
        else:
            others, rid = pids[:-1], pids[-1]

        if self.__isCapture:
            return __tag + '%s)%s' % (specialNuclearParticleIDModule.specialNuclearParticleID(IDsPoPsModule.photon, mode), 
                rid)

        residual = {rid: products2.pop(rid)}

        if residual[rid] > 1:
            others.append(rid)
            products2[rid] = residual[rid] - 1

        for pid in others:
            spid = specialNuclearParticleIDModule.specialNuclearParticleID(pid, mode)
            if products2[pid] == 1:
                __tag += '%s' % spid
            else:
                __tag += '%d%s' % (products2[pid], spid)

        __tag += ')%s' % rid

        return __tag

    def shortTag(self, initialProducts, mode=specialNuclearParticleIDModule.Mode.familiar):
        """
        This method returns the short tag for of the reaction which is the list of products between the ',' and ')'
        of the tag.

        :param initialProducts:     A boolean, if **True**, the initial products tag is returned; otherwise, the final products tag is returned.
        :param mode:                Determines whether id is used for the the light particles (e.g., for 'H3', is is 'H3', 'h3' or 't').
        """

        tag1 = self.tag(initialProducts, mode)

        return tag1.split(',')[-1].split(')')[0]
