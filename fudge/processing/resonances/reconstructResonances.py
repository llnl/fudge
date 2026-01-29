#! /usr/bin/env python3
#encoding: utf-8

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
See ENDF-102 (endf documentation) appendix D for equations

Basic usage:
------------

    if protare is a reactionSuite instance,::

        >>protare.reconstructResonances(tolerance=0.01)

    This reconstructs all resolved/unresolved resonance sections. In each section,
    the results are accurate under linear interpolation to tolerance of 1% or better.
    Sections are summed together and then added to the appropriate background cross
    sections to form new pointwise, lin-lin cross sections that are stored in the
    appropriate reaction inside protare.

Alternate uses:
---------------

    You can skip the final step (adding background cross sections) and just get the
    resonance parameter contribution to the cross section by doing::

        >>xsecs = reconstructResonances.reconstructResonances(protare, 0.01)

    Here,::

        * xsecs = dictionary containing cross sections: {'total':XYs, 'elastic':XYs, ...}

    each cross section is an XYs class instance, containing data and also axes with units

    Another option would be to only reconstruct a single section::

        >> resCls = reconstructResonances.RMcrossSection(protare, energyUnit='eV')   # for Reich_Moore
        >> energy_grid = s.generateEnergyGrid()
        >> crossSections = s.getCrossSection( energy_grid )
        # the input to getCrossSection is the energy (or list of energies) in self.energyUnit
        # crossSections are returned as a dictionary {'total':,'elastic':,'capture':,'fission':,}

        # improve grid to desired tolerance for linear interpolation:
        >> new_energy_grid, new_crossSections, messages = s.refineInterpolation(energy_grid, crossSections, tolerance=0.01)

"""
import sys
import numpy, collections, math
import abc

from pqu import PQU as PQUModule

from fudge import reactionSuite as reactionSuiteModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData.distributions import angular as angularModule

from fudge.resonances import resolved as resolvedModule, common as commonResonancesModule

import fudge.processing.resonances.getCoulombWavefunctions as getCoulombWavefunctions
from fudge.core.math.fudgemath import RoundToSigFigs

from xData import XYs1d as XYs1dModule

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.families import nuclide as nuclidePoPsModule

if sys.version_info.major > 2:
    import multiprocessing
    try:
        multiprocessing.set_start_method("fork")
    except RuntimeError:    # multiprocessing context may have already been set
        pass

debug = False   # recommend setting to True before debugging (disables multiprocessing)

VERBOSE = False

FISSIONCHANNEL, NEUTRONCHANNEL, CPCHANNEL, GAMMACHANNEL, COMPETITIVECHANNEL = range(5)

class ChannelDesignator:
    __slots__=[ 'l', 'J', 'reaction', 'index', 's', 'gfact', 'particleA', 'particleB', 'Xi', 'isElastic', 'channelClass', 'useRelativistic', 'eliminated', 'externalRMatrix' ]

    def __init__(self, l, J, reaction, index, s, gfact=None, particleA=None, particleB=None, Xi=0.0, isElastic=None, channelClass=None, useRelativistic=False, eliminated=False, externalRMatrix=None ):
        """
        Simple widget to specify a channel

        RRBaseClass.setResonanceParametersByChannel() creates these mostly and getAllowedTotalSpins() is used to determine the spins in there.
        So, if you're wondering why s is really 2*s, unlike l and J, check in those routines.

        :param l: orbital angular momentum of channel, an int
        :param J: total angular momentum of a channel, can be an int or 1/2 int (use the _eq_ below to avoid floating point comparisons)
        :param reaction: string, name of the reaction associated with this channel
        :param index: int, an optional index of the channel (useful if using this in an unordered dictionary)
        :param s: total spin of a channel, can be an int or 1/2 int (use the _eq_ below to avoid floating point comparisons)
        :param gfact: statistical factor
        :param particleA: the projectile, either a string ('n') or None if you don't need it
        :param particleB: the target, either a string ('Pu239') or None if you don't need it
        :param Xi: threshold energy in lab frame
        :param isElastic: flag to denote that this channel is, in fact, also incoming channel so this is an elastic channel
        :param channelClass: nature of channel: gamma, fission, neutron or charged particle
        :param useRelativistic: flag for using relativistic kinematics
        :param eliminated: flag for whether to eliminate this channel using Reich-Moore approximation (useful for gamma channels)
        :param externalRMatrix: optional external function (SAMMY or Froehner parameterization)
        :return:
        """
        self.l=l
        self.J=J
        self.reaction=reaction
        self.index=index
        self.s=s
        self.gfact=gfact
        self.particleA=particleA
        self.particleB=particleB
        self.Xi=Xi
        self.isElastic=isElastic
        self.channelClass=channelClass
        self.useRelativistic=useRelativistic
        self.eliminated=eliminated
        self.externalRMatrix=externalRMatrix

    def __eq__(self, other):
        if self.gfact is None or other.gfact is None:
            gfact_identical = self.gfact is other.gfact
        else:
            gfact_identical = math.isclose(self.gfact, other.gfact, rel_tol=1e-6, abs_tol=1e-5)
        return self.l==other.l and \
               int(2.0*self.J)==int(2.0*other.J) and \
               int(2.0*self.s)==int(2.0*other.s) and \
               gfact_identical and \
               self.particleA==other.particleA and \
               self.particleB==other.particleB and \
               self.isElastic==other.isElastic and \
               self.channelClass==other.channelClass and \
               self.useRelativistic==other.useRelativistic and \
               self.eliminated==other.eliminated and \
               self.externalRMatrix==other.externalRMatrix

    def __hash__(self):
        theHash=0
        for s in ChannelDesignator.__slots__:
            if s == 'Xi': continue # is computable and is real, so lousy for hashing
            if s == 'gfact': continue # is computable and is real, so lousy for hashing
            if s == 'reaction': continue # LRF=7 calls 'elastic' 'n+Target', so the name is kinda meaningless
            if s == 'eliminated': continue # different approximations treat this field differently
            if s == 'index': continue # this isn't used except for printing out
            theHash += hash(getattr(self, s))
        return theHash

    def __repr__(self):
        return "ChannelDesignator(l=%i, J=%s, reaction='%s', index=%i, s=%s, gfact=%s, particleA='%s', particleB='%s', Xi=%s, isElastic=%s, channelClass=%s, useRelativistic=%s, eliminated=%s)" \
               % (self.l, self.J, self.reaction, self.index, self.s, self.gfact, self.particleA, self.particleB, str(self.Xi), self.isElastic, self.channelClass, self.useRelativistic, self.eliminated)

    def has_same_J(self,other):
        return abs(self.J-other.J)<0.001

    def is_open(self,Ein): return Ein >= self.Xi


def spins_equal(s1,s2): return int(2.0*s1)==int(2.0*s2)


def getResonanceReconstructionClass( formalism ):
    if formalism.moniker == resolvedModule.BreitWigner.moniker:
        gndsClass = resolvedModule.BreitWigner
        proc = {
            gndsClass.Approximation.SingleLevel: SLBWcrossSection,
            gndsClass.Approximation.MultiLevel: MLBWcrossSection
        }[ formalism.approximation ]
        return proc
    elif formalism.moniker == resolvedModule.RMatrix.moniker:
        return RMatrixLimitedcrossSection
    else:
        raise TypeError("Don't recognize resonance type %s" % formalism.moniker)


def getAllowedTotalSpins( L, S, useFactor2Trick=True ):
    """
    Returns a list of allowed J values from summing angular momenta L and S, where
    :math:`\\vec{J}=\\vec{L}+\\vec{S}`
    which implies :math:`|L-S| \\leq J \\leq L+S`

    The useFactor2Trick flag tells the routine whether we are summing real angular momenta or momenta * 2.
    If the useFactor2Trick flag is true, then momenta are really momenta*2, meaning they can be pure integers,
    even if the real momenta refer to 1/2-integer values (e.g. spin).  The default is to useFactor2Trick because
    most C/C++/Fortran codes that compute angular momentum-stuff use the trick so they can use integer math.
    Also, it makes the use of the Python range() function possible.
    """
    if useFactor2Trick: return list( range( abs(L-S), L+S+2, 2 ) )
    else: return [ j/2.0 for j in getAllowedTotalSpins( int(L*2), int(S*2) ) ]


def reconstructResonances(reactionSuite, tolerance, verbose=False, significantDigits=None, disableUnresolvedWidthInterpolation=False):
    """
    Reconstruct all resonance cross sections (resolved and unresolved) in reactionSuite,
    and add results together for full (resonance region) pointwise cross section.

    Optional arguments:

        - 'tolerance': fractional tolerance, used to refine interpolation grid.
          That is, if tolerance = 0.001, points are added until lin-lin interpolation is good to 0.1% everywhere.

        - 'verbose': print status messages during reconstruction

        - 'significantDigits': Controls how many digits can be used to represent the energy grid. For example,
          if significantDigits=4 the resulting energy grid can contain 1.034, 1.035, 1.036 but not 1.0345.
          Using significantDigits=8 should allow data to be written back to ENDF-6 without loss of precision.

        - 'disableUnresolvedWidthInterpolation': set 'True' to reproduce old ENDF recommendation to interpolate
          in cross sections rather than widths. That recommendation was reversed before ENDF-VIII release

    Note that this does not add the background contribution if any (ENDF-MF3 portion).
    Use the method reactionSuite.reconstructResonances which adds the background contribution.
    """

    energyUnit = reactionSuite.domainUnit
    egrids, xsecs = [], []

    def resolvedReconstruct( resolvedForm ):
        """
        Helper function to reconstruct one region, used for single & multiple regions as well as URR ensembles.
        """

        resCls = getResonanceReconstructionClass(resolvedForm)
        reconstructClass = resCls( resolvedForm, enableAngDists=False, verbose=verbose, energyUnit=energyUnit )
        egrid = reconstructClass.generateEnergyGrid()
        if significantDigits is not None:
            egrid = sorted( set( RoundToSigFigs(egrid,significantDigits).tolist() ) )
        xsecs_now = reconstructClass.getCrossSection( egrid )
        egrid, xsecs_now, messages = reconstructClass.refineInterpolation(numpy.array(egrid), xsecs_now, tolerance,
                significantDigits=significantDigits)
        if verbose:
            for message in messages: print (message)
        return egrid, xsecs_now

    # Resolved resonance reconstruction
    if reactionSuite.resonances.resolved:
        # For multiple regions, we need to do each region separately, then add them to the unified xs table & egrid
        if isinstance(reactionSuite.resonances.resolved.evaluated, commonResonancesModule.EnergyIntervals):
            if verbose:
                print( "      WARNING! Multiple resolved/unresolved energy intervals are deprecated\n"
                        +"        and should be consolidated into single section")
            for interval in reactionSuite.resonances.resolved.evaluated:
                egrid_now, xsecs_now = resolvedReconstruct( interval.evaluated )
                # merge with 'full' energy grid:
                egrids.append( egrid_now )
                xsecs.append( xsecs_now )
        # Single region, everything goes on unified grid
        else:
            resolvedForm = reactionSuite.resonances.resolved.evaluated
            egrid_now, xsecs_now = resolvedReconstruct( resolvedForm )
            # merge with 'full' energy grid:
            egrids.append( egrid_now )
            xsecs.append( xsecs_now )

    # Unresolved region has more complex logic, we may do one of a couple of things:
    #   a) do nothing
    #   b) reconstruct the average cross section assuming SLBW resonances (traditional ENDF)
    #   c) reconstruct the average xs and the PDF for the xs
    if reactionSuite.resonances.unresolved:
        if reactionSuite.resonances.unresolved.evaluated.useForSelfShieldingOnly: # don't reconstruct
            if verbose:
                print ("Skipping unresolved: for self shielding only")
        else:
            unresolvedForm = reactionSuite.resonances.unresolved.evaluated
            reconstructClass = URRcrossSection( unresolvedForm, verbose )
            egrid, interpolateWidths = reconstructClass.generateEnergyGrid()
            xsecs_now = reconstructClass.getCrossSection( egrid )

            if not disableUnresolvedWidthInterpolation:
                # ENDF manual now states that we should interpolate cross sections rather than parameters for URR.
                egrid, xsecs_now, messages = reconstructClass.refineInterpolation(numpy.array(egrid), xsecs_now, tolerance)
                if verbose:
                    for message in messages: print (message)

            # meld with resolved region:
            egrids.append( egrid )
            xsecs.append( xsecs_now )

    # 'xsecs' list now holds one or more regions. Convert to pointwise or Regions1d cross section instance:
    xsecs_final = {}
    crossSectionAxes = crossSectionModule.defaultAxes( energyUnit )
    if len(xsecs)==1:   # only one region: treat as XYs1d
        for key in xsecs[0]:
            xsecs_final[key] = crossSectionModule.XYs1d( axes = crossSectionAxes, data=(egrids[0], xsecs[0][key]), dataForm="XsAndYs" )
    else:               # multiple regions: treat as Regions1d
        for key in xsecs[0]:
            pwxs = crossSectionModule.Regions1d( axes = crossSectionAxes )
            for idx in range(len(xsecs)):
                xys = crossSectionModule.XYs1d( axes = crossSectionAxes, data=(egrids[idx], xsecs[idx][key]), dataForm="XsAndYs" )
                pwxs.append( xys )
            xsecs_final[key] = pwxs

    return xsecs_final

def reconstructAngularDistributions(reactionSuite, tolerance=None, verbose=False):
    """
    Reconstruct all pure two-body angular distributions from the resonance region.
    For MLBW and RM, that means 'elastic' channel only, but R-Matrix evaluations can support additional channels.
    SLBW cannot be used for angular distributions.

    Returns a Python dict.  They key is the reaction and the value is a reconstructed
    angularModule.XYs2d instance containing Legendre expansions for each incident energy.
    """

    energyUnit = reactionSuite.domainUnit
    angdists = {}
    angdistRegions = []
    egrids = []

    # Helper function to compute the per-region angular distributions
    def resolvedReconstruct( resolvedForm ):

        # Get the correct class
        resCls = getResonanceReconstructionClass(resolvedForm)
        reconstructClass = resCls( resolvedForm, enableAngDists=True, verbose=verbose, energyUnit=energyUnit )

        # Deal with potential thresholds by splitting up the energy grid at the thresholds
        fullEgrid = reconstructClass.generateEnergyGrid()
        if hasattr(reconstructClass,'_thresholds'):
            thresholdIndices = [0]
            for Xi in reconstructClass._thresholds:
                if Xi in fullEgrid:
                    thresholdIndices.append(list(fullEgrid).index(Xi))
            thresholdIndices.append(-1)
            subgrids = [fullEgrid[thresholdIndices[i]:thresholdIndices[i+1]] for i in range(len(thresholdIndices[:-1]))]
        else: subgrids=[fullEgrid]

        # Now do the calculation & merge the different regions
        angularResults={}
        for egrid in subgrids:
            print("Working egrid: [%s, ..., %s]"%(str(egrid[0]),str(egrid[-1])))
            thisAngularResults = reconstructClass.getAngularDistribution( egrid, keepL0Term=True, renormalize=True )
            if tolerance:
                raise NotImplementedError("Refining interpolation grid for angular distributions")
            for reaction, results in thisAngularResults.items():
                if reaction not in angularResults: angularResults[reaction] = []
                angularResults[reaction] += list( zip( egrid, numpy.array( results ).squeeze().T ) )
        return angularResults

    if reactionSuite.resonances.resolved:

        # Compute the per-region angular distributions & make a list of regional distributions
        if isinstance(reactionSuite.resonances.resolved.evaluated, commonResonancesModule.EnergyIntervals):
            raise NotImplementedError("Angular reconstruction for multiple resolved regions, untested")
            if verbose:
                print( "      WARNING! Multiple resolved/unresolved energy regions are deprecated\n"
                        +"        and should be consolidated into single section")
            for RRidx in range(len(reactionSuite.resonances.resolved.regions)):
                formalism = reactionSuite.resonances.resolved.regions[RRidx].evaluated
                angdistRegions.append( resolvedReconstruct( formalism ) )
        else:
            formalism = reactionSuite.resonances.resolved.evaluated
            angdistRegions.append( resolvedReconstruct( formalism ) )

        # Check that egrids don't overlap
        if len( angdistRegions ) ==1: pass
        else:
            for iRegion, this_region in enumerate(angdistRegions[:-1]):
                next_region = angdistRegions[iRegion+1]
                if this_region[-1][0] > next_region[0][0]:
                    raise ValueError( "this_region[-1][0] > next_region[0][0] for region %i" % iRegion )

        # 'angdistRegions' list holds one or more regions.  Convert to one Legendre coefficient table
        # right now only do elastic scattering (OK for MLBW & RM, but (unimplemented) RML can do more)
        angdists_final={}
        angularAxes = angularModule.defaultAxes(energyUnit)
        if len(angdistRegions) == 1:    # only one region, treat as XYs2d
            for key in angdistRegions[0]:
                form = angularModule.XYs2d(axes=angularAxes)
                for i1, E_Bs in enumerate(angdistRegions[0][key]):
                    form.append(angularModule.Legendre(coefficients=E_Bs[1], outerDomainValue=E_Bs[0], axes=angularAxes))
                angdists_final[key] = form
        else:
            raise NotImplementedError("Still TBD")

    return angdists_final


"""
@blockwise: function decorator for improving performance in resolved region.
Each 'getCrossSection' and 'getAngularDistribution' method is wrapped by this function.
If we have lots of incident energies, this splits up a calculation using the multiprocessing module.
May need to tweak the 'NEmax' and 'numTasks' variables for best performance.
"""
def blockwise(function):
    NEmax = 50000
    def wrapped(self,E,**kwargs):
        if numpy.isscalar(E):
            NE = 1
            E = numpy.array([[E]])
        else:
            NE = len(E)
            E = numpy.array(E).reshape(NE,1)

        if debug: # disable multiprocessing
            return function(self,E,**kwargs)
        elif NE < NEmax: # faster to run directly
            return function(self,E,**kwargs)
        else:
            from multiprocessing import Process, Queue
            queue = Queue()

            def enqueue_result(elist, Slice, queue):
                # perform the calculation, put result in the queue
                try:
                    result = function(self, elist, **kwargs)
                    queue.put( (Slice, result) )
                except Exception as err:
                    queue.put( (Slice, err) )

            jobs = []
            from fudge.defaults import numTasks # number of processes to spawn
            errorsEncountered = False
            result = None
            # how many energies does each process calculate?
            chunk = int(numpy.ceil((NE+0.) / numTasks))
            for pid in range(numTasks): # start the calculations
                Slice = slice(pid*chunk, (pid+1)*chunk)
                p = Process(target=enqueue_result, args=(E[Slice], Slice, queue))
                jobs.append(p)
                p.start()

            for pid in range(numTasks): # collect results as they finish
                Slice, vals = queue.get()
                if isinstance(vals, Exception):
                    print("Error in child process #%d: %s" % (pid,vals))
                    errorsEncountered = True
                    continue
                if function.__name__ == 'getCrossSection':
                    if result is None:
                        result = dict.fromkeys( vals )
                        for key in result: result[key] = numpy.zeros(NE)
                    for key in result:
                        result[key][Slice] = vals[key]
                elif function.__name__ == 'getAngularDistribution':
                    if result is None:
                        result = dict.fromkeys( vals )
                        for key in result: result[key] = [numpy.zeros((NE,1)) for L in range(len(vals[key]))]
                    for key in result:
                        for L in range(len(vals[key])):
                            result[key][L][Slice] = vals[key][L]
                else:
                    raise NotImplementedError("Blockwise computation for function %s" % function.__name__)

            # allow child processes to exit:
            for job in jobs: job.join()

            if errorsEncountered:
                raise Exception("Aborting due to errors in child processes")
            return result
    return wrapped


class ResonanceReconstructionBaseClass(abc.ABC):
    """Base class common to resolved and unresolved resonance reconstruction."""

    def __init__( self, resonanceForm, **kw ) :

        reactionSuite = resonanceForm.findClassInAncestry(reactionSuiteModule.ReactionSuite)
        try:
            neutron = resonanceForm.PoPs['n']
            self.neutronMass_amu = neutron.getMass('amu')
        except KeyError:
            self.neutronMass_amu = 1.00866491574

        self.energyUnit = reactionSuite.domainUnit

        self.projectile = reactionSuite.projectile
        projectile = resonanceForm.PoPs[reactionSuite.projectile]
        self.projectileMass_amu = projectile.getMass('amu')
        if isinstance(projectile, nuclidePoPsModule.Particle):
            self.projectileSpin = projectile.nucleus.spin[0].value
        else:
            self.projectileSpin = projectile.spin[0].value

        targetID = reactionSuite.target
        if( resonanceForm.PoPs.hasAlias( targetID ) ) : targetID = resonanceForm.PoPs[targetID].pid
        self.target = targetID
        target = resonanceForm.PoPs[targetID]
        if( isinstance( target, nuclidePoPsModule.Particle ) ) :
            nuclide = target.ancestor[0]
            self.targetMass_amu = nuclide.getMass('amu')
        else :
            self.targetMass_amu = target.getMass('amu')
        self.targetSpin = target.nucleus.spin[0].value

        self.targetToNeutronMassRatio = self.targetMass_amu / self.neutronMass_amu

        self.kConstant1 = math.sqrt(2 * PQUModule.PQU(self.projectileMass_amu, 'amu * c**2').getValueAs(self.energyUnit))
        self.kConstant1 /= PQUModule.PQU(1,'hbar * c').getValueAs('%s * 1e-14 * m' % self.energyUnit)   # 1e-14 m is sqrt(b).
        self.kConstant2 = self.kConstant1 * self.targetMass_amu / ( self.projectileMass_amu + self.targetMass_amu )

        self.energyFactor = math.sqrt(PQUModule.PQU(1, self.energyUnit).getValueAs('eV'))

        self.reactionSuite = reactionSuite

    @abc.abstractmethod
    def getCrossSection(self, E): pass

    def k(self, energy):
        """
        Calculates the wave number for the projectile in the center of mass given the projectile energy in the lab frame.
        The returned value is in b**-1/2.
        """

        return self.kConstant2 * numpy.sqrt(energy)

    """ refer to SAMMY manual page 9 for penetration/shift/phase factor equations: """
    def penetrationFactor(self, L, rho):
        if L==0: return rho
        elif L==1: return rho**3 / (1+rho**2)
        elif L==2: return rho**5 / (9+3*rho**2+rho**4)
        elif L==3: return rho**7 / (225+45*rho**2+6*rho**4+rho**6)
        elif L==4: return rho**9 / (11025+1575*rho**2+135*rho**4+10*rho**6+rho**8)
        elif L>4:
            # find it recursively
            P_Lminus1 = self.penetrationFactor(L-1,rho)
            S_Lminus1 = self.shiftFactor(L-1,rho)
            return rho**2 * P_Lminus1 / ((L-S_Lminus1)**2 + P_Lminus1**2)

    def shiftFactor(self, L, rho):
        """ calculate shift factor used in SLBW and MLBW formalisms """
        if L==0: return 0.0*rho   # rho may be an array
        elif L==1: return -1.0 / (1+rho**2)
        elif L==2: return -(18.0+3*rho**2) / (9+3*rho**2+rho**4)
        elif L==3: return -(675.0+90*rho**2+6*rho**4) / (225+45*rho**2+6*rho**4+rho**6)
        elif L==4: return -( (44100.0+4725*rho**2+270*rho**4+10*rho**6) /
                (11025+1575*rho**2+135*rho**4+10*rho**6+rho**8) )
        elif L>4:
            P_Lminus1 = self.penetrationFactor(L-1,rho)
            S_Lminus1 = self.shiftFactor(L-1,rho)
            return (rho**2 * (L-S_Lminus1) / ((L-S_Lminus1)**2+P_Lminus1**2) - L)

    def phi(self, L, rho):
        # calculate hard-sphere phase-shift
        if   L==0: return rho
        elif L==1: return rho-numpy.arctan(rho)
        elif L==2: return rho-numpy.arctan(3*rho/(3-rho**2))
        elif L==3: return rho-numpy.arctan(rho*(15-rho**2)/(15-6*rho**2))
        elif L==4: return rho-numpy.arctan(rho*(105-10*rho**2)/(105-45*rho**2+rho**4))
        elif L>4:
            P_Lminus1 = self.penetrationFactor(L-1,rho)
            S_Lminus1 = self.shiftFactor(L-1,rho)
            phi_Lminus1 = self.phi(L-1,rho)
            return phi_Lminus1 - numpy.arctan(P_Lminus1 / (L-S_Lminus1))

    # for use after the cross section has been calculated on the initial grid:
    def refineInterpolation(self, egrid, xsecs, tolerance=0.01, significantDigits=None):
        """ generateEnergyGrid may not give a fine enough grid to linearly interpolate to desired tolerance.
        My solution to that: for all consecutive points (x0,y0), (x1,y1) and (x2,y2) do a linear interpolation between
        (x0,y0) and (x2,y2). If the interpolation doesn't agree with (x1,y1) within tolerance,
        subdivide up the region by adding two more calculated points.  Iterate until interpolation agrees within tolerance.

        This means that in the end we will have more points than required for given tolerance.
        The results can be thinned (thinning implemented in xData.XYs1d)
        """
        def checkInterpolation(x,y):
            """
            Does a linear interpolation of each point from its two nearest neighbors,
            checks where the interpolation grid is insufficient and returns a list of indices (in array x)
            where additional points are needed.

            @type x: numpy.multiarray.ndarray
            @type y: numpy.multiarray.ndarray
            :return:
            """
            m = (y[2:]-y[:-2])/(x[2:]-x[:-2])
            delta_x = x[1:-1] - x[:-2]
            b = y[:-2]

            # add first and last points back in, for easier comparison with
            # original y values:
            interpolated = numpy.zeros_like(x)
            interpolated[1:-1] = m*delta_x+b
            interpolated[0] = y[0]; interpolated[-1] = y[-1]

            # find where original and interpolated grids differ by more than tolerance
            # silence div/0 warnings for this step, since xsc = 0 case is explicitly handled below
            with numpy.errstate( divide='ignore', invalid='ignore' ):
                delt = interpolated / y
                mask = (delt>1+tolerance) + (delt<1-tolerance) # boolean array

            badindices = numpy.arange(len(mask))[ mask ]    # points where finer mesh is needed

            # switch to absolute convergence condition for very small cross sections (i.e. near thresholds):
            smallXSec = 1e-50
            zeros = (y[badindices-1]<smallXSec) + (y[badindices]<smallXSec) + (y[badindices+1]<smallXSec)
            if any(zeros):
                ignore = []
                for idx in badindices[ zeros ]:
                    if abs(y[idx]-y[idx-1])<1e-3 and abs(y[idx+1]-y[idx])<1e-3:
                        mask[ idx ] = False
                        ignore.append( idx )
                badindices = list(badindices)
                for idx in ignore[::-1]:
                    badindices.remove(idx)

            return badindices

        messages = []
        reactionDone = dict.fromkeys(xsecs, False)
        n_iter = 0
        addedPoints = 0
        while True:
            newIdx = set()
            for key in xsecs:
                if not any(xsecs[key]):
                    continue
                if not reactionDone[key]:
                    badindices = checkInterpolation(egrid,xsecs[key])
                    if len(badindices)==0: reactionDone[key] = True
                    newIdx.update( badindices )
            newIdx = sorted(newIdx)

            mask = numpy.zeros( len(egrid), dtype=bool )
            mask[newIdx] = True
            midpoints = (egrid[:-1]+egrid[1:])/2

            energies_needed = sorted( set( list(midpoints[mask[:-1]]) + list(midpoints[mask[1:]]) ) )
            if significantDigits is not None:
                rounded = RoundToSigFigs( energies_needed, significantDigits )
                rounded = set( rounded.tolist() )
                rounded.difference_update( egrid.tolist() ) # remove any rounded values that were already computed
                energies_needed = sorted( rounded )

            if len(energies_needed)==0:    # success!
                break
            if n_iter > 20:
                messages.append("Iteration limit exceeded when refining interpolation grid!")
                break
            n_iter += 1
            if debug:
                print( "Iteration #%d: adding %d points to interpolation grid" % (n_iter, len(energies_needed)) )
            addedPoints += len(energies_needed)
            newY = self.getCrossSection( energies_needed )
            # merge new x/y values with original list:
            fulllist = numpy.append( egrid, energies_needed )
            order = numpy.argsort( fulllist, kind='merged' )
            egrid = fulllist[order]
            for key in xsecs:
                xsecs[key] = numpy.append( xsecs[key], newY[key] )[ order ]

        messages.append("%i points were added (for total of %i) to achieve tolerance of %s%%" %
            (addedPoints, len(egrid), tolerance*100))
        return egrid, xsecs, messages


#### base class for resolved resonance reconstruction ####
class RRBaseClass(ResonanceReconstructionBaseClass, abc.ABC):

    def __init__(self, resolvedForm, lowerBound, upperBound, **kw):

        ResonanceReconstructionBaseClass.__init__(self, resolvedForm, **kw)
        """ store resonance parameters in convenient structure for quick cross section calculations: """

        self.RR = resolvedForm

        self.lowerBound = lowerBound
        self.upperBound = upperBound

        self.removeNegatives = True # set to False to allow returning negative cross section values

    def sortLandJ(self):
        """
        SLBW, MLBW and Reich_Moore formalisms have similar structure
        it's convenient to sort their resonances by L and J
        This method should NOT be used for R-Matrix Limited
        """
        # make a unified table of resonance parameter data, convert to self.energyUnit if necessary
        nRes = len( self.RR.resonanceParameters.table )
        params = ('energy','L','J','channelSpin','totalWidth','neutronWidth','captureWidth',
                'fissionWidth','fissionWidthB')
        EU = self.energyUnit
        units = (EU,'','','',EU,EU,EU,EU,EU)
        data = [self.RR.resonanceParameters.table.getColumn( quant,unit ) for quant,unit in zip(params,units) ]
        missingColumns = []
        for i in range(len(data)):
            if data[i] is None:
                if params[i] in ('energy','L','J','neutronWidth','captureWidth'):
                    missingColumns.append(params[i])
                data[i] = [0]*nRes
        if missingColumns:
            raise KeyError(f"Required columns {', '.join(missingColumns)} missing from BreitWigner parameters table!")
        table = numpy.array( data )

        # sort resonances by L and J, store parameters in numpy arrays
        self.Ls = []
        Llist = table[ params.index('L') ]
        for L in range(int(max( Llist ))+1):
            Lres = table[ :, table[ params.index('L') ]==L ]
            Jlist = Lres[ params.index('J') ]
            Js = []
            for J in sorted(set(Jlist)):
                LJres = Lres[ :, Lres[ params.index('J') ]==J ]
                spinList = LJres[ params.index('channelSpin') ]
                spins = []
                for spin in sorted(set(spinList)):
                    spinRes = LJres[ :, LJres[ params.index('channelSpin') ]==spin ]
                    energies = spinRes[ params.index('energy') ]
                    neutronWidth = ( spinRes[ params.index('neutronWidth') ] /
                        self.penetrationFactor(L,self.rho(abs(energies))) )
                    # FIXME for easier-to-read code, should I wait until getCrossSection
                    # to adjust the neutronWidth for penetrability at resonances?

                    spindict = {
                            'channelSpin': spin,
                            'energy': energies,
                            'neutronWidth': neutronWidth,
                            'captureWidth': spinRes[params.index('captureWidth')],
                            'fissionWidth': spinRes[params.index('fissionWidth')],
                            'fissionWidthB': spinRes[params.index('fissionWidthB')],
                            'shiftFactor': 0.5*self.shiftFactor(L,self.rho(numpy.abs(energies))),
                            }
                    spins.append( spindict )
                Jdict = {
                        'J': J,
                        'gfact': (2.0*abs(J)+1)/((2*self.projectileSpin+1)*(2*self.targetSpin+1)),
                        'channelSpins': spins,
                        }
                Js.append( Jdict )

            Ldict = { 'L': L, 'Js': Js }
            self.Ls.append( Ldict )

        # for each L, gfactors should sum to 2*L+1.
        self.missingGfactor = {}
        for L in self.Ls:
            gsum = sum( [J['gfact'] * len(J['channelSpins']) for J in L['Js']] )
            self.missingGfactor[L['L']] = 2*L['L']+1 - gsum

        # for energy grid generation:
        self._energies = table[ params.index('energy') ]
        totalWidths = table[ params.index('totalWidth') ]
        if not any(totalWidths): # total widths aren't specified in Reich_Moore case
            totalWidths = ( table[params.index('neutronWidth')] + table[params.index('captureWidth')]
                    + table[params.index('fissionWidth')] + table[params.index('fissionWidthB')] )
        self._widths = totalWidths

    def getAPByChannel(self,c,trueOrEffective='true'):
        """get the channel radius, rho. If L is specified try to get L-dependent value"""
        if self.RR.calculateChannelRadius:
            a = 0.123 * self.targetMass_amu**( 1. / 3. ) + 0.08  # eq. D.14 in ENDF manual, a in b^-1/2
        else:
            if self.RR.getScatteringRadius().isEnergyDependent():
                a = self.RR.getScatteringRadius().getValueAs('10**-12*cm')   # a in b^-1/2
            else:
                a = self.RR.getScatteringRadius().getValueAs( '10**-12*cm', L=c.l )   # a in b^-1/2
        return a

    def rho(self, E):
        """
        Get the channel radius rho. May be computed from a simple expression or from the tabulated scattering radius.

        :param E: incident energy grid, necessary if the tabulated scattering radius is energy-dependent
        :param forceUseTabulated: override 'calculateChannelRadius' flag (should be used for computing phase shift)
        :return:
        """
        if self.RR.calculateChannelRadius:
            a = 0.123 * self.targetMass_amu**( 1. / 3. ) + 0.08  # eq. D.14 in ENDF manual, a in b^-1/2
        else:
            if self.RR.getScatteringRadius().isEnergyDependent():
                a = self.RR.getScatteringRadius().getValueAs('10**-12*cm', energy_grid=numpy.ravel(E))   # a in b^-1/2
                a = numpy.array(a).reshape(E.shape)
            else:
                a = self.RR.getScatteringRadius().getValueAs( '10**-12*cm' )   # a in b^-1/2
        return self.k(E) * a # dimensionless

    def rhohat(self, E):
        """
        Get the effective channel radius rhohat.

        :param E: incident energy grid, necessary if the tabulated scattering radius is energy-dependent
        :return:
        """
        if self.RR.getHardSphereRadius().isEnergyDependent():
            a = self.RR.getHardSphereRadius().getValueAs('10**-12*cm', energy_grid=numpy.ravel(E))   # a in b^-1/2
            a = numpy.array(a).reshape(E.shape)
        else:
            a = self.RR.getHardSphereRadius().getValueAs( '10**-12*cm' )   # a in b^-1/2
        return self.k(E) * a # dimensionless

    def generateEnergyGrid(self, lowBound=None, highBound=None, stride=1):
        """ Create an initial energy grid by merging a rough mesh for the entire region (~10 points / decade)
        with a denser grid around each resonance. For the denser grid, multiply the total resonance width by
        the 'resonancePos' array defined below. """
        energies, widths = self._energies, self._widths
        thresholds = []
        if hasattr(self, '_thresholds'):
            thresholds = self._thresholds
        if lowBound is None:
            lowBound = self.lowerBound
        if highBound is None:
            highBound = self.upperBound
        # ignore negative resonances
        for lidx in range(len(energies)):
            if energies[lidx] > 0: break
        energies = energies[lidx:]
        widths = widths[lidx:]
        # generate grid for a single peak, should be good to 1% using linear interpolation using default stride
        resonancePos = numpy.array([
            0.000e+00, 1.000e-03, 2.000e-03, 3.000e-03, 4.000e-03, 5.000e-03, 6.000e-03, 7.000e-03, 8.000e-03, 9.000e-03, 1.000e-02, 2.000e-02,
            3.000e-02, 4.000e-02, 5.000e-02, 6.000e-02, 7.000e-02, 8.000e-02, 9.000e-02, 1.000e-01, 1.100e-01, 1.200e-01, 1.300e-01, 1.400e-01,
            1.500e-01, 1.600e-01, 1.700e-01, 1.800e-01, 1.900e-01, 2.000e-01, 2.100e-01, 2.200e-01, 2.300e-01, 2.400e-01, 2.500e-01, 2.600e-01,
            2.800e-01, 3.000e-01, 3.200e-01, 3.400e-01, 3.600e-01, 3.800e-01, 4.000e-01, 4.200e-01, 4.400e-01, 4.600e-01, 4.800e-01, 5.000e-01,
            5.500e-01, 6.000e-01, 6.500e-01, 7.000e-01, 7.500e-01, 8.000e-01, 8.500e-01, 9.000e-01, 9.500e-01, 1.000e+00, 1.050e+00, 1.100e+00,
            1.150e+00, 1.200e+00, 1.250e+00, 1.300e+00, 1.350e+00, 1.400e+00, 1.450e+00, 1.500e+00, 1.550e+00, 1.600e+00, 1.650e+00, 1.700e+00,
            1.750e+00, 1.800e+00, 1.850e+00, 1.900e+00, 1.950e+00, 2.000e+00, 2.050e+00, 2.100e+00, 2.150e+00, 2.200e+00, 2.250e+00, 2.300e+00,
            2.350e+00, 2.400e+00, 2.450e+00, 2.500e+00, 2.600e+00, 2.700e+00, 2.800e+00, 2.900e+00, 3.000e+00, 3.100e+00, 3.200e+00, 3.300e+00,
            3.400e+00, 3.600e+00, 3.800e+00, 4.000e+00, 4.200e+00, 4.400e+00, 4.600e+00, 4.800e+00, 5.000e+00, 5.200e+00, 5.400e+00, 5.600e+00,
            5.800e+00, 6.000e+00, 6.200e+00, 6.400e+00, 6.500e+00, 6.800e+00, 7.000e+00, 7.500e+00, 8.000e+00, 8.500e+00, 9.000e+00, 9.500e+00,
            1.000e+01, 1.050e+01, 1.100e+01, 1.150e+01, 1.200e+01, 1.250e+01, 1.300e+01, 1.350e+01, 1.400e+01, 1.450e+01, 1.500e+01, 1.550e+01,
            1.600e+01, 1.700e+01, 1.800e+01, 1.900e+01, 2.000e+01, 2.100e+01, 2.200e+01, 2.300e+01, 2.400e+01, 2.500e+01, 2.600e+01, 2.700e+01,
            2.800e+01, 2.900e+01, 3.000e+01, 3.100e+01, 3.200e+01, 3.300e+01, 3.400e+01, 3.600e+01, 3.800e+01, 4.000e+01, 4.200e+01, 4.400e+01,
            4.600e+01, 4.800e+01, 5.000e+01, 5.300e+01, 5.600e+01, 5.900e+01, 6.200e+01, 6.600e+01, 7.000e+01, 7.400e+01, 7.800e+01, 8.200e+01,
            8.600e+01, 9.000e+01, 9.400e+01, 9.800e+01, 1.020e+02, 1.060e+02, 1.098e+02, 1.140e+02, 1.180e+02, 1.232e+02, 1.260e+02, 1.300e+02,
            1.382e+02, 1.550e+02, 1.600e+02, 1.739e+02, 1.800e+02, 1.951e+02, 2.000e+02, 2.100e+02, 2.189e+02, 2.300e+02, 2.456e+02, 2.500e+02,
            2.600e+02, 2.756e+02, 3.092e+02, 3.200e+02, 3.469e+02, 3.600e+02, 3.892e+02, 4.000e+02, 4.200e+02, 4.367e+02, 4.600e+02, 4.800e+02,
            5.000e+02, 6.000e+02, 7.000e+02, 8.000e+02, 9.000e+02, 1.000e+03][::stride])

        grid = []
        # get the midpoints (on log10 scale) between each resonance:
        # emid = [lowBound] + list(10**( ( numpy.log10(energies[1:])+numpy.log10(energies[:-1]) ) / 2.0)) + [highBound]
        # or get midpoints on linear scale:
        emid = [lowBound] + [(e1+e2)/2.0 for e1, e2 in zip(energies[1:], energies[:-1])] + [highBound]
        for e, w, lowedge, highedge in zip(energies, widths, emid[:-1], emid[1:]):
            points = e-w*resonancePos
            grid += [lowedge] + list(points[points>lowedge])
            points = e+w*resonancePos[1:]
            grid += list(points[points < highedge])
        # also add rough grid, to cover any big gaps between resonances, should give at least 10 points per decade:
        npoints = int(numpy.ceil(numpy.log10(highBound)-numpy.log10(lowBound)) * 10)
        grid += list(numpy.logspace(numpy.log10(lowBound), numpy.log10(highBound), npoints))[1:-1]
        grid += [lowBound, highBound, 0.0253]   # region boundaries + thermal
        # if threshold reactions present, add dense grid at and above threshold
        for threshold in thresholds:
            grid += [threshold]
            grid += list(threshold + resonancePos * 1e-2)
        grid = sorted(set(grid))
        # toss any points outside of energy bounds:
        grid = grid[grid.index(lowBound) : grid.index(highBound)+1]
        return grid

    def setResonanceParametersByChannel( self, multipleSScheme='ENDF', useReichMooreApproximation=False, Ein=None, warnOnly=False ):
        """
        Reorganize member data into channels (relies heavily on groundwork in sortLandJ).
        """
        self.nResonances = len( self.RR.resonanceParameters.table )
        self.lMax = 0   # updated below
        self.channelConstantsBc = []
        allowedSs = getAllowedTotalSpins( self.projectileSpin, self.targetSpin, useFactor2Trick=False )
        warnings=[]

        # Make a unified table of resonance parameter data, convert to self.energyUnit if necessary
        lList = []
        params = ('energy','L','J','channelSpin','totalWidth','neutronWidth','captureWidth', 'fissionWidth','fissionWidthB')
        EU = self.energyUnit
        units = (EU,'','','',EU,EU,EU,EU,EU)
        data = [ self.RR.resonanceParameters.table.getColumn( quant,unit ) for quant,unit in zip(params,units) ]
        for i in range(len(data)):
            if data[i] is None: data[i] = [0]*self.nResonances

        def addOrUpdateDict( theDict, theKey, theValue ):
            if not theKey in theDict: theDict[ theKey ] = {}
            theDict[ theKey ][ theValue[0] ] = theValue[1]

        particlePairs={
            'capture':('gamma', self.RR.rootAncestor.compound_nucleus),
            'elastic':(self.RR.rootAncestor.projectile, self.RR.rootAncestor.target)}
        pA102 = particlePairs['capture'][0]
        pB102 = particlePairs['capture'][1]
        pA2 = particlePairs['elastic'][0]
        pB2 = particlePairs['elastic'][1]

        # Now make the dictionary of channels
        channelDict = collections.OrderedDict()
        for iR in range( self.nResonances ):
            ER = data[0][iR]
            L =  data[1][iR]
            J =  data[2][iR]
            S =  None # not given in MLBW, RM or RML, I just list it here for completeness
            GT = data[4][iR]
            GN = data[5][iR]
            GG = data[6][iR]
            GF = data[7][iR]
            GX = data[8][iR]
            JisAllowed = False
            self.lMax = max( self.lMax, L )
            if L not in lList: lList.append(L)
            allowedSByThisJl = []
            for SS in allowedSs:
                if abs(J) in getAllowedTotalSpins( L, SS, useFactor2Trick=False ):
                    JisAllowed = True
                    allowedSByThisJl.append( SS )
            if not JisAllowed:
                if warnOnly:
                    from fudge import warning
                    warnings.append(warning.InvalidAngularMomentaCombination(L, S, abs(J), 'all channels'))
                else:
                    raise ValueError( "Invalid total angular momentum and/or orbital momentum: cannot couple up to J = "
                                      +str(J)+" with I = "+str(self.targetSpin)+", i = " + str(self.projectileSpin) + " and L = "+str(L) )
            gfact = (2.0*abs(J)+1.)/((2.*self.projectileSpin+1)*(2.*self.targetSpin+1.))

            # Way recommended by NJOY:
            #   * One of valid S's get width, other(s) get zero width
            #   * Zero width channels left in sum to get potential scattering correct
            #   * Gets good (cs from U) to (cs from Fudge) comparison; good (cs from U) to (cs from BB)
            if multipleSScheme == 'NJOY':
                if allowedSByThisJl:
                    addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'elastic', 0, allowedSByThisJl[0], gfact, pA2, pB2, 0.0, True, NEUTRONCHANNEL, False, False ), ( iR, GN ) )
                    addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'capture', 0, allowedSByThisJl[0], gfact, pA102, pB102, 0.0, False, GAMMACHANNEL, False, useReichMooreApproximation ), ( iR, GG ) )
                    if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'fission', 0, allowedSByThisJl[0], gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False ), ( iR, GF ) )
                    if GX != 0.0: addOrUpdateDict(channelDict, ChannelDesignator(L, abs(J), 'competitive', 0, allowedSByThisJl[0], gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False), (iR, GX))
                    for SS in allowedSByThisJl[1:]:
                        addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'elastic', 0, SS, gfact, pA2, pB2, 0.0, True, NEUTRONCHANNEL, False, False ), ( iR, 0.0 ) )
                        addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'capture', 0, SS, gfact, pA102, pB102, 0.0, False, GAMMACHANNEL, False, useReichMooreApproximation ), ( iR, 0.0 ) )
                        if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'fission', 0, SS, gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False ), ( iR, 0.0 ) )
                        if GX != 0.0: addOrUpdateDict(channelDict, ChannelDesignator(L, abs(J), 'competitive', 0, SS, gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False), (iR, 0.0))

            # Way mandated by ENDF manual:
            #   * Width divided between channels with valid S's based on sign of J:
            #       ** if J>0, S is max possible value
            #       ** if J<0, S is min possible value
            #   * All channels left in sum to get potential scattering correct
            #   * Gets poor (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB)
            elif multipleSScheme == 'ENDF':
                if allowedSByThisJl:
                    allowedSByThisJl.sort()
                    if len(allowedSByThisJl)==1: SS=allowedSByThisJl[0]
                    else: SS=allowedSByThisJl[int(abs(J)==J)]
                    addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'elastic', 0, SS, gfact, pA2, pB2, 0.0, True, NEUTRONCHANNEL, False, False ), ( iR, GN ) )
                    addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'capture', 0, SS, gfact, pA102, pB102, 0.0, False, GAMMACHANNEL, False, useReichMooreApproximation ), ( iR, GG ) )
                    if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'fission', 0, SS, gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False ), ( iR, GF ) )
                    if GX != 0.0: addOrUpdateDict(channelDict, ChannelDesignator(L, abs(J), 'competitive', 0, SS, gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False), (iR, GX))

            # Ignore problem:
            #   * Ignore spin of channels and the fact may be multiple valid spins
            #   * Gets best (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB)
            else :
                addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'elastic', 0, S, gfact, pA2, pB2, 0.0, True, NEUTRONCHANNEL, False, False ), ( iR, GN ) )
                addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'capture', 0, S, gfact, pA102, pB102, 0.0, False, GAMMACHANNEL, False, useReichMooreApproximation ), ( iR, GG ) )
                if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, abs(J), 'fission', 0, S, gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False ), ( iR, GF ) )
                if GX != 0.0: addOrUpdateDict(channelDict, ChannelDesignator(L, abs(J), 'competitive', 0, S, gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False), (iR, GX))

        # Take care of extra l's needed to get potential scattering correct
        # Make a single elastic channel with a zero width resonance at the energy of the first real resonance
        # This won't add to the resonance scattering part, but it will make a potential scattering part.
        # We add in channels for all spins and total angular momenta to make sure get the angular momentum algebra correct.
        for l in range( 0, self.lMax+1 ):
            for s in allowedSs:
                for j in getAllowedTotalSpins( l, s, useFactor2Trick=False ):
                    newChannel = ChannelDesignator( l, j, 'elastic', 0, s, (2.0*abs(j)+1)/((2*self.projectileSpin+1)*(2*self.targetSpin+1)), pA2, pB2, 0.0, True, NEUTRONCHANNEL, False, False )
                    if not newChannel in channelDict: addOrUpdateDict( channelDict, newChannel, ( 0, 0 ) )

        # Set up the channel->resonance mappings
        self.allChannels = channelDict
        self.channels = collections.OrderedDict()  # just the kept ones
        self.eliminatedChannels = collections.OrderedDict()
        for c in self.allChannels:
            if c.eliminated:    self.eliminatedChannels[ c ] = self.allChannels[ c ]
            else:               self.channels[ c ] = self.allChannels[ c ]
        self.nChannels = len( self.channels )
        self.identityMatrix = numpy.identity( self.nChannels, dtype=complex )

        if warnOnly: return warnings

    def resetResonanceParametersByChannel(self, multipleSScheme='ENDF', useReichMooreApproximation=False, Ein=None ): pass

    def penetrationFactorByChannel(self,c,Ein):
        if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]:
            if type(Ein)==numpy.ndarray:
                return numpy.ones(len(Ein))
            else:
                return 1.0
        return self.penetrationFactor(c.l, self.rho(Ein))

    def phiByChannel(self,c,Ein):
        if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]: return 0.0
        return self.phi( c.l, self.rho(Ein) )

    def shiftFactorByChannel(self,c,Ein):
        if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]: return 0.0
        return self.shiftFactor(c.l, self.rho(Ein))

    def getChannelConstantsBc( self ):
        raise NotImplementedError( "Override in derived classes if used" )

    def getRMatrix( self, Ein ): raise NotImplementedError( "Override in derived classes if used" )

    def getKMatrix( self, Ein ): raise NotImplementedError( "Override in derived classes if used" )

    def getAMatrix( self, Ein ): raise NotImplementedError( "Override in derived classes if used" )

    def getL0Matrix( self, Ein ):
        """
        Get the L0 matrix of Froehner:

            ..math::
                {\\bf L^0}_{cc'} = \\delta_{cc'} (L_c-B_c)

        where

            ..math::
                L_c = S_c + i P_c

        """
        if self.channelConstantsBc == []: self.channelConstantsBc = self.getChannelConstantsBc()
        L0 = numpy.zeros( (len(Ein), self.nChannels, self.nChannels), dtype=complex )
        for ic,c in enumerate(self.channels):
            rho = self.rho( Ein-c.Xi )
            shift = self.shiftFactor( c.l, rho ).flatten()
            penet = self.penetrationFactor( c.l, rho ).flatten()

            # Note, you have to add the numpy.array's this way to make a complex number
            # because the numpy array trickery doesn't work with the complex() or
            # numpy.complex() constructors
            L0[:,ic,ic] = shift + 1j* penet - self.channelConstantsBc[ ic ]
        return L0

    def getXMatrix( self, Ein ):
        """
        Get the X matrix for use in computing W.  X is:

            ..math::
                {\\bf X}_{cc'} = P^{-1/2}_c ( ( {\\bf I} - {\\bf R}{\\bf L^0} )^{-1}{\\bf R} )_{cc'} P_{c'}^{-1/2}\\delta_{JJ'}
        """
        L0 = self.getL0Matrix( Ein )
        penFact = numpy.ones( ( len(Ein), self.nChannels ), dtype=float )
        X = numpy.zeros( ( len(Ein), self.nChannels, self.nChannels ), dtype=complex )

        # Precompute phase factor, sqrt(Penetrability)
        for ic, c in enumerate( self.channels ):
            if c.channelClass in [NEUTRONCHANNEL, CPCHANNEL]: penFact[:,ic] = numpy.sqrt( L0[:,ic,ic].imag )

        # Compute the (unscaled) X matrix
        # Hafta be real careful about combining numpy's broadcasting rules with numpy's linear algebra
        R = self.getRMatrix( Ein )
        for iE in range(len(Ein)):
            thisR = numpy.array(R[iE])
            thisL0 = numpy.array(L0[iE])
            X[iE] = numpy.linalg.inv( numpy.array(self.identityMatrix) - thisR * thisL0 ) * thisR

        # Apply the sqrt(Penetrability) rescaling
        for ic1, c1 in enumerate( self.channels ):
            for ic2, c2 in enumerate( self.channels ):
                if c1.J != c2.J: continue
                X[:, ic1, ic2] *= penFact[:, ic1] * penFact[:, ic2]

        return X

    def getWMatrix( self, Ein ):
        """
        Get the W matrix for use in computing U.  W is:

            ..math::
                {\bf W} = {\bf I} + 2i{\bf X}
        """
        W = numpy.zeros( ( len(Ein), self.nChannels, self.nChannels ), dtype=complex )
        W[ :, numpy.identity( self.nChannels, dtype=bool ) ] = 1
        W += 2.0j * self.getXMatrix( Ein )
        return W

    def getEiPhis( self, Ein, useTabulatedScatteringRadius=True, enableExtraCoulombPhase=False ):
        """
        The phase factor for the collision matrix:

            ..math::
                \\Omega_c = e^{-\\varphi_c}
        """

        # Precompute phase factor
        eiphis = []
        for ic, c in enumerate( self.channels ):
            if c.reaction == 'elastic':
                # For calculating phi, default is to use tabulated scattering radius:
                if useTabulatedScatteringRadius:
                    if self.RR.getScatteringRadius().isEnergyDependent():
                        rhohat = self.RR.getScatteringRadius().getValueAs('10*fm', Ein-c.Xi) * self.k(Ein-c.Xi)
                    else:
                        rhohat = self.RR.getScatteringRadius().getValueAs('10*fm') * self.k(Ein-c.Xi)
                else: rhohat = self.rho(Ein-c.Xi)
                phi = self.phi( c.l, rhohat )
                eiphis.append( numpy.exp( -1j * phi ) )
            else:
                eiphis.append( numpy.ones( Ein.shape, dtype=complex ) )
        return eiphis

    def getScatteringMatrixU( self, Ein, useTabulatedScatteringRadius=True, enableExtraCoulombPhase=False ):
        """
        Compute the scattering matrix using
        """
        # Make sure all channels are open at all requested energies
        if not numpy.all([c.is_open(Ein) for c in self.channels]): raise ValueError("One or more channels are not open on all energies in requested energy grid")

        # Initialize lists and matrices
        eiphis =  self.getEiPhis( Ein,
                                  useTabulatedScatteringRadius=useTabulatedScatteringRadius,
                                  enableExtraCoulombPhase=enableExtraCoulombPhase )
        W = self.getWMatrix( Ein )

        # This adds in the phase factors to get us U
        U = numpy.zeros( ( len(Ein), self.nChannels, self.nChannels ), dtype=complex )
        for ic1 in range( self.nChannels ):
            for ic2 in range( self.nChannels ):
                U[:, ic1, ic2 ] = W[:, ic1, ic2 ] * (eiphis[ic1][:] * eiphis[ic2][:]).flatten()
        return U

    def getScatteringMatrixT( self, Ein, useTabulatedScatteringRadius=True, enableExtraCoulombPhase=False ):
        # Make sure all channels are open at all requested energies
        if not numpy.all([c.is_open(Ein) for c in self.channels]): raise ValueError("One or more channels are not open on all energies in requested energy grid")

        T = numpy.zeros( ( len(Ein), self.nChannels, self.nChannels ), dtype=complex )
        U = self.getScatteringMatrixU( Ein,
                                       useTabulatedScatteringRadius=useTabulatedScatteringRadius,
                                       enableExtraCoulombPhase=enableExtraCoulombPhase )
        for ic1 in range( self.nChannels ):
            for ic2 in range( self.nChannels ):
                T[:, ic1, ic2 ] = (ic1==ic2) - U[:, ic1, ic2 ]
        return T

    def getLMax(self,maxLmax=10):
        """
        LMax is determined by the behavior of the Blatt-Biedenharn Zbar coefficients.  Inside each one, there is
        a Racah coefficient and a Clebsh-Gordon coefficient.  The CG coefficient looks like this::

            ( l1 l2  L )
            (  0  0  0 )

        So, this means two things.  First, the CG coeff (and hence Zbar) will be zero if l1+l2+L=odd.
        Second, the maximum value of L will be l1max+l2max.  Hence, Lmax=2*lmax.
        """
        return min( max( [ 2*xL['L'] for xL in self.Ls ] ), maxLmax )

    def getParticleParities(self,rxn):
        """
        Wrapper around getParticleSpinParities to extract just the parities

        :param rxn: the reaction string for this resonanceReaction (or equivalent).  We'll process this string to
                    out what the light particle and heavy "residual" is and then look up their JPi.
        :return: a tuple: (JPi_light, JPi_heavy).  JPi itself is a tuple (J,Pi).
                 J is float (either integer or 1/2 integer) and Pi is -1 or +1.
        """
        r=self.getParticleSpinParities(rxn)
        return r[0][1], r[1][1]

    def getParticleSpins(self,rxn):
        """
        Wrapper around getParticleSpinParities to extract just the spins

        :param rxn: the reaction string for this resonanceReaction (or equivalent).  We'll process this string to
                    out what the light particle and heavy "residual" is and then look up their JPi.
        :return: a tuple: (JPi_light, JPi_heavy).  JPi itself is a tuple (J,Pi).
                 J is float (either integer or 1/2 integer) and Pi is -1 or +1.
        """
        r=self.getParticleSpinParities(rxn)
        return r[0][0], r[1][0]

    def getParticleSpinParities(self, rxn):
        """
        6 cases:
            rxn=='capture'
            rxn=='elastic' (equivalent to 'projectile + target')
            rxn=='fission' or 'fissionA' or 'fissionB'
            rxn=='competitive'
            rxn=='something + something' where one of somethings is photon (same as 'capture')
            rxn=='something + something' where neither something is a photon

            A word about capture channels:
            Because of the use of the Reich-Moore approximation, photon channels may be quasi-channels consisting
            of many photon channels lumped together.  This makes determination of the residual nucleus spin & parity
            difficult since it really corresponds to a bunch of states.  Here, we just compute the lowest spin
            consistent with the target & photon and leave it up to the code to override these values as needed.

            A word about fission and competitive channels:
            Here we don't really know or need the spins & parities of the particles.  Both are many channels
            lumped together, in a way that makes capture look simple.  We'll just return the elastic parameters
            under the assumption that they're not needed.

        :param rxn: the reaction string for this resonanceReaction (or equivalent).  We'll process this string to
                    out what the light particle and heavy "residual" is and then look up their JPi.
        :return: a tuple: (JPi_a, JPi_b).  JPi itself is a tuple (J,Pi).
                 J is float (either integer or 1/2 integer) and Pi is -1 or +1.
                 It doesn't matter whether particle a or b is heaviest
        """

        def get_spin_parity(_p):
            s, par = None, None
            if hasattr(_p, 'nucleus'):
                try:
                    s = _p.nucleus.spin.float('hbar')
                except KeyError as err:
                    raise KeyError(str(err) + ' for particle ' + _p.id)
                except IndexError as err:
                    # Sometimes (mainly with MLBW and SLBW), there's no spin given for a nucleus.
                    # In that case, an IndexError gets raised.  We have to punt.
                    pass
                try:
                    par = _p.nucleus.parity[0].value
                except IndexError as err:
                    # Sometimes (mainly with MLBW and SLBW), there's no parity given for a nucleus.
                    # In that case, an IndexError gets raised.  We have to punt.
                    pass
            else:
                try:
                    s=_p.spin.float('hbar')
                    par=_p.parity[0].value
                except IndexError as err:
                    #print(_p)
                    raise err
            return s, par

        # Treat elastic, fission and competitive channels as elastic
        # for elastic is appropriate, for others it doesn't matter
        if rxn=='elastic' or 'fission' in rxn or "competitive" in rxn:
            spinA, parityA = get_spin_parity(self.RR.PoPs[self.reactionSuite.projectile])
            target = self.RR.PoPs.final(self.target)
            spinB, parityB = get_spin_parity(target)
            return (spinA, parityA), (spinB, parityB)

        # Capture...
        if rxn=='capture' or 'photon' in rxn:
            spinA = min(self.projectileSpin+0.5, abs(self.projectileSpin-0.5))
            parityA = 1
            compound = self.RR.PoPs.final(self.reactionSuite.compound_nucleus)
            spinB, parityB = get_spin_parity(compound)
            return (spinA, parityA), (spinB, parityB)

        # Everything else...
        if ' + ' in rxn: # gotta search
            spinA, parityA = None, None
            spinB, parityB = None, None
            pa,pb=rxn.split(' + ')
            for p in self.RR.PoPs:
                if pa==p.id: spinA,parityA=get_spin_parity(p)
                if pb==p.id: spinB,parityB=get_spin_parity(p)

            return (spinA, parityA), (spinB, parityB)

        else:
            raise ValueError('Cannot determine spins for reactants in reaction "%s"'%rxn)

    @blockwise
    def getAngularDistribution(self, E, keepL0Term=True, renormalize=True, outChannelNames=None,
                               useTabulatedScatteringRadius=True, enableExtraCoulombPhase=False ):
        """

        :param numpy.array(type=float) E:
        :param bool keepL0Term:
        :param bool renormalize:
        :param list(str) outChannelNames:
        :param bool useTabulatedScatteringRadius:
        :param bool enableExtraCoulombPhase:
        :return:
        :rtype: dict
        """
        from numericalFunctions import angularMomentumCoupling as nf_amc

        # Make sure the incident energies are in the correct form for vectorization
        if type(E) in [float, numpy.float64]: E = numpy.array(E).reshape(1,1)
        else: E = numpy.array(E).reshape(len(E),1)

        # Make sure all channels are open at all requested energies
        self.resetResonanceParametersByChannel( Ein=E )
        if not numpy.all([c.is_open(E) for c in self.channels]): raise ValueError("One or more channels are not open on all energies in requested energy grid")

        # Get the names of the input channel and all output channels to consider
        channelNames = set([c.reaction for c in self.channels])
        inChannelName = '%s + %s' % (self.projectile, self.target)
        if inChannelName not in channelNames:   # FIXME change MLBW and R-M to label channels by particle pairs like LRF=7
            if 'elastic' in channelNames:
                inChannelName = 'elastic'
        if outChannelNames is None:
            outChannelNames = set([chan.reaction for chan in self.channels if not chan.eliminated and chan.reaction not in ['capture','fission','fissionA','fissionB','competitive']])

        # Set up variables, first the easy ones
        k = self.k(E)                                           # Relative momentum, comes out in b^-1/2
        B = {key:[] for key in outChannelNames}                 # List of Blatt-Biedenharn coefficients, used in Legendre series expansion of the angular distribution
        results = {key:[] for key in outChannelNames}
        relTol = 1e-8                                           # Relative tolerance used to decide when to cut off Legendre series
        Lmax = self.getLMax( )                                  # Lmax also used to decide when to cut off Legendre series

        # Now, the T matrix, the most important variable of all...
        T = self.getScatteringMatrixT( E,                       # The scattering matrix
                                       useTabulatedScatteringRadius=useTabulatedScatteringRadius,
                                       enableExtraCoulombPhase=enableExtraCoulombPhase )

        # The in channel stuff...
        projSpin, targSpin = self.getParticleSpins(inChannelName)
        prefactor = 1.0/k/k/(2*projSpin+1)/(2*targSpin+1)       # Prefactor for cross section, gets units, main energy dependence, and spin dep. prefactor
        allowedSpins = getAllowedTotalSpins( targSpin,          # Possible values of the total channel spin, based on the projectile spin (projSpin) and target spin (targSpin)
                                             projSpin,
                                             False )
        incidentChannels = [(i,chan) for i,chan in enumerate(self.channels)
                            if not chan.eliminated and chan.reaction == inChannelName]

        # Loop over out channels
        for reactionp in outChannelNames:
            ejectSpin, resSpin = self.getParticleSpins(reactionp)
            allowedSpinsp = getAllowedTotalSpins( resSpin,          # Possible values of the total channel spin, based on the ejectile spin (ejectSpin) and residual nucleus spin (resSpin)
                                                  ejectSpin,
                                                  False )
            outgoingChannels = [(i,chan) for i,chan in enumerate(self.channels)
                                if not chan.eliminated and chan.reaction == reactionp]

            # Main algorithm
            # Some words on the looping over L here.  There are 3 reasons to keep looping:
            #   1. We have less than 3 Legendre moments (L=0,1,2).  At threshold we have L=0 only, but otherwise had better have more.
            #   2. If we hit LMax or all the biggest L moments are too small, then stop, otherwise keep going
            #   3. We end on an odd L index.  That means Lmax must be even because LMax = last L-1
            # By the way, we must have an even LMax.  See the getLMax() documentation for reasoning.
            L = 0
            while len(B[reactionp])<=3 or ( numpy.any( numpy.abs(B[reactionp][-1])>=relTol*B[reactionp][0] ) and L<=Lmax ) or L%2==0:
                B[reactionp].append(numpy.zeros_like(E))
                for S in allowedSpins:
                    for Sp in allowedSpinsp:
                        spinFactor =  pow( -1, S - Sp ) / 4.0
                        for i1, c1 in incidentChannels:
                            if c1.s is not None and not spins_equal(c1.s, S): continue # c1.s != S: continue
                            for i2, c2 in incidentChannels:
                                if c2.s is not None and not spins_equal(c2.s, S): continue # c2.s != S: continue
                                Z = nf_amc.zbar_coefficient( int(2.*c1.l), int(2.*c1.J), int(2.*c2.l), int(2.*c2.J), int(2.*S), int(2.*L) )
                                if Z == 0.0: continue
                                for i1p, c1p in outgoingChannels:
                                    if c1.J != c1p.J: continue
                                    if c1p.s is not None and not spins_equal(c1p.s, Sp): continue # c1p.s != Sp: continue
                                    for i2p, c2p in outgoingChannels:
                                        if c2.J != c2p.J: continue
                                        if c2p.s is not None and not spins_equal(c2p.s, Sp): continue # c2p.s != Sp: continue
                                        Zp = nf_amc.zbar_coefficient( int(2*c1p.l), int(2*c1p.J), int(2*c2p.l), int(2*c2p.J), int(2*Sp), int(2*L) )
                                        if Zp == 0.0: continue
                                        thisB = spinFactor * Z * Zp * ( T[:,i1,i1p].conj() * T[:,i2,i2p] ).real
                                        B[reactionp][-1] += thisB[:,numpy.newaxis]
                L += 1

            for l,coefs in enumerate( B[reactionp] ):
                results[reactionp].append( prefactor / ( 2.0*l+1.0 ) * coefs )
            if renormalize:
                L0term=results[reactionp][0]
                if any(L0term<1e-16): L0term[L0term<1e-16]=1e-16 #takes care of questionable behavior near threshold
                results[reactionp] = [ x/L0term for x in results[reactionp] ]

        return results

    def getAverageQuantities(self, nBins=1, binScheme='linspacing', computeUncertainty=False, verbose=False):
        """
        Computes average widths and level spacings from the set of resonance parameters in self, on a per-channel basis

        The averages are computed in equal lethargy bins starting at the lowest resonance energy in a sequence up to the
        upperBound of the resonance region.  I tried to keep on average 10 resonances/logrithmic bin so I can get a
        reasonable average.

        :param nBins: the number of resonances per logrithmic bin to aim for, on average
        :param binScheme: the scheme to use to determine the binning
        :param computeUncertainty: toggle the calculation of the uncertainty of the quantities
        :return: a dictionary of results, sorted by channel.  If an entry is None, this indicates it could not be computed
        """
        from numpy import mean,var
        results={}
        if not binScheme in ['linspacing','logspacing']:
            raise ValueError("Allowed binning schemes are 'linspacing' and 'logspacing'")

        for c in self.allChannels:
            # initialize the acumulators
            numResonancesUsed=0
            results[c] = {'energyGrid': [], 'widths': [], 'spacings': []}

            # determine if we have resonances to average
            if len(self.allChannels[c]) == 0: continue
            try:                ERmin = min( [ xx for xx in self.allChannels[c].values() if xx>0 ] )
            except ValueError:  ERmin = 1.0

            # get the resonances for this channel
            res=[ (self._energies[iR], self.allChannels[c][iR]) for iR in self.allChannels[c] ]
            res.sort(key=lambda x:x[0])

            # compute the number of bins to use in practice, must be at least 1 bin
            nBins=max(min(nBins,len(res)),1)

            # compute the bin boundaries
            lowerBound = min(max(self.lowerBound, 25.85e-3), ERmin)
            upperBound = self.upperBound
            if binScheme=='logspacing':
                bins = list(numpy.logspace(math.log10(lowerBound), math.log10(upperBound), nBins+1))
            else:
                bins = list(numpy.linspace(lowerBound, upperBound, nBins+1))
            results[c]['energyGrid'] = bins

            if verbose: print(c.reaction, c.l, c.s, c.J, 'bins:',bins)

            # compute averages in a bin
            for i1 in range(len(bins)-1):
                # note filters out negative energy resonances and high energy resonances whose tails make a background
                EGam = [(self._energies[i], v) for i, v in self.allChannels[c].items() if
                        self._energies[i] >= bins[i1] and self._energies[i] <= bins[i1+1]]
                EGam.sort(key=lambda x: x[0])

                # Compute mean level spacing
                spacingsList=[abs(EGam[i + 1][0] - EGam[i][0]) for i in range(len(EGam) - 1)]
                widthList=[x[1] for x in EGam]

                # initialization
                actualResonancesPerBin =  len(widthList)
                numResonancesUsed      += actualResonancesPerBin
                D     = None
                Gamma = None

                if verbose: print(spacingsList, widthList)

                # enough resonances to make an average D, Gamma
                if len(spacingsList)>1:
                    if computeUncertainty:
                        D     = PQUModule.PQU(value=mean(spacingsList), uncertainty=math.sqrt(var(spacingsList)), unit='eV')
                        Gamma = PQUModule.PQU(value=mean(widthList), uncertainty=math.sqrt(var(widthList)), unit='eV')
                    else:
                        D     = PQUModule.PQU(value=mean(spacingsList), unit='eV')
                        Gamma = PQUModule.PQU(value=mean(widthList), unit='eV')

                # only one width, we can "average" it, but the spacing is meaningless
                elif len(widthList)==1:
                    Gamma = PQUModule.PQU(value=widthList[0], unit='eV')

                # save results
                results[c]['spacings'].append(D)
                results[c]['widths'].append(Gamma)

        return results

    def getTransmissionCoefficients(self, resonancesPerBin=10, computeUncertainty=False, method='weakCoupling'):
        """
        Compute the effective transmission coefficient in the resolved resonance regime.
        The method argument controls the schemes to do this.

        In each case, define the variable tau as::

            ..math::
                x_c = \\pi\\overline{\\Gamma}_c/D

        These are the options:

            * `method=='weakCoupling'`:  This is the traditional weak coupling limit::

                ..math::
                    T_c=2x_c

            * `method=='1stOrder'`: This is the same as the weak coupling limit.

            * `method=='2ndOrder'`: This is the 2nd order correction to the weak coupling limit::

                ..math::
                    T_c=2x_c(1-x_c)

            * `method=='sumRule'`: Use Moldauer's sum rule to extract the transmission coefficients directly from
              the RRR tables P.A. Moldauer Phys. Rev. Lett. 19, 1047-1048 (1967).  The equation is::

                ..math::
                    T_c=2x_c\\left[\\sqrt{1+x_c^2}-x_c\\right]

            * `method=='opticalModel'`: Use Moldauer's 'optical model' transmission coefficient::

                ..math::
                    T_c=1-\\exp{(-2x_c)}

            * `method=='SPRT'`: Use the

                ..math::
                    T_c=\\frac{2x_c}{ (1 + x_c/2)^2 + (P_c R^\\inf_c)^2 }

        :param computeUncertainty: toggle the calculation of the uncertainty of the transmission coefficients
        :param resonancesPerBin: we'll try to get at least this number of resonances in each energy bin
        :param method: method to compute the transmission coefficient
        :return: a dictionary of results, sorted by channel.  Tc entries are all PQU's or None if it cannot be computed
        """
        allowedMethods=["weakCoupling", '1stOrder', '2ndOrder', 'sumRule', 'opticalModel', 'SPRT']
        if method not in allowedMethods:
            raise ValueError('Transmission coefficient calculation method must be one of '+', '.join(allowedMethods))
        nBins=min([len(self.allChannels[c]) for c in self.allChannels if len(self.allChannels[c]) > 0])//resonancesPerBin
        aves=self.getAverageQuantities(nBins=nBins,computeUncertainty=computeUncertainty)
        results={}
        for c in aves:
            results[c]={'Tc':[],'energyGrid':aves[c]['energyGrid']}
            for i in range(len(aves[c]['energyGrid'])-1):
                if not aves[c]['spacings'][i] or not aves[c]['widths'][i]:
                    results[c]['Tc'].append(None)
                    continue
                try:
                    twoPiGOverD=2.0*math.pi*aves[c]['widths'][i]/aves[c]['spacings'][i]
                except ZeroDivisionError:
                    continue
                if method in ["weakCoupling",'1stOrder']:
                    results[c]['Tc'].append(twoPiGOverD)
                elif method == '2ndOrder':
                    results[c]['Tc'].append(twoPiGOverD*(1.0-twoPiGOverD/2.0))
                elif method == 'sumRule':
                    if twoPiGOverD < 1e-16: results[c]['Tc'].append(PQUModule.PQU(value=0.0,unit=''))
                    else: results[c]['Tc'].append(0.5*twoPiGOverD*twoPiGOverD*( math.sqrt(1.0+4.0/twoPiGOverD/twoPiGOverD) - 1.0))
                elif method=='opticalModel':
                    results[c]['Tc'].append(PQUModule.PQU(1.0-math.exp(-twoPiGOverD), ""))
                else: #method=="SPRT"
                    xx=1.0+twoPiGOverD*twoPiGOverD/4.0
                    results[c]['Tc'].append(twoPiGOverD/xx/xx)
        return results

    def getPoleStrength(self, computeUncertainty=False):
        """
        Computes the (usu. neutron) pole strength from Eq. (209) from JEFF Report 18::

            ..math::
                s_c(E) = \\frac{\\left<\\gamma_c^2\\right>}{D} = \\frac{\\overline{\\Gamma_c}(ER)}{2 D P_c(ER)}

        In the weak coupling limit, we could also use the transmission coefficient :math:`s_c(E)=T_c(ER)/4\\pi P_c`

        :param computeUncertainty:
        :return:  a dictionary of results, sorted by channel.
                  Pole strength entries are all PQU's or None if it cannot be computed
        """
        import numbers
        results = {}
        aves = self.getAverageQuantities(computeUncertainty=computeUncertainty)
        for c in aves:
            results[c]={'sc':[], 'energyGrid':aves[c]['energyGrid']}
            for i in range(len(aves[c]['energyGrid'])-1):
                if (not aves[c]['spacings'][i] or not aves[c]['widths'][i]) or len(self.allChannels[c]) == 0:
                    results[c]['sc'].append(None)
                    continue
                Emin=aves[c]['energyGrid'][i]
                Emax=aves[c]['energyGrid'][i+1]
                redWidths=self.getReducedWidths(c, Emin=Emin, Emax=Emax)
                if isinstance(redWidths, numbers.Number):  # Why do I need this test?
                    redWidths = [redWidths]
                if len(redWidths) == 0:
                    continue
                gg=numpy.mean([g*g for g in redWidths])
                ggUnc=numpy.sqrt(numpy.var([g*g for g in redWidths]))
                results[c]['sc'].append(PQUModule.PQU(gg, unit='eV', uncertainty=ggUnc)/aves[c]['spacings'][i])
        return results

    def getScatteringLength(self, L=None):
        """
        Said Mughabghab (Atlas of Neutron Resonances, 5th edition (2006) p. 19, eq. (1.54))
        defines R' as the s-wave scattering length, calculable from the potential
        scattering cross section sigPot = 4 Pi (R')^2.

        We compute the potential scattering cross section at E=1e-5 eV via

        ..math::
            \\sigma_{pot} = \\sum_L 4\\pi g_J\\frac{\\sin(\\phi_L(E))}{k^2}

        As Mughabghab only gives the S-wave scattering length, we provide the option of choosing an L as an argument
        to the function.  If L is supplied by the user, only that L is used in the calculation of R', otherwise all
        L's are used.

        R' is given in in b**1/2 and should be close to AP.  It won't be exactly the same because of the
        contribution of distant levels.  Froehner (JEFF Report 18) Eq. (214) and the Atlas Eq. (1.51) give

        ..math::
            R' = AP*(1-R^\\infty)

        where AP is the scattering length given in the file for the R-matrix parameterization and R^\\infty is the
        distant level parameter, related to the pole strength (call getPoleStrength) as a Hilbert transform, or
        by direct summation such as Atlas Eq. (1.52):

        ..math::
            R^\\infty=\\sum_n \\frac{\\gamma_n^2}{E_n-E}

        :param L: either an integer (the angular momentum to use) or None (all angular momentum)
        :return: R', the scattering radius, in b**1/2
        """
        E=numpy.array([ PQUModule.PQU(1e-5, 'eV').getValueAs(self.energyUnit) ])
        sigPotL = []
        for c in self.allChannels:
            if not c.isElastic: continue
            if L is None or c.l == L:
                sigPotL.append( c.gfact * 4.0 * numpy.pi * ( numpy.sin( self.phiByChannel(c,E) ) )**2 / self.k(E)**2 )
        sigPot = sum(sigPotL)
        return numpy.sqrt(sigPot/4.0/numpy.pi)[0]

    def getReducedWidths(self, channel, Emin=0.0, Emax=None):
        """
        Convert the widths corresponding to `channel` to reduced widths in the usual way (JEFF Report 18, Eq. (174))::

            ..math::
                \\gamma_c = \\sqrt{\\Gamma_c(E_R)/2P_c(E_R)}

        Here :math:`\\gamma_c` is the reduced width, :math:`\\Gamma_c` is the regular width, tabulated in the ENDF file,
        :math:`E_R` is the resonance energy (ignoring shifts), and :math:`P_c(E)` is the penetrability factor.

        :param Emin: The minimum resonance energy to convert to reduced width (in eV), defaults to 0.0.
                     If set to None, all resonances will be converted.
        :param Emax: The maximum resonance energy to convert to reduced width (in eV), defaults to None.
                     If set to None, all resonances will be converted.
        :param channel: The ChannelDesignator of the widths to be converted.
        :return: the list of reduced widths corresponding to the resonances in the range [Emin, Emax]
        """
        reducedWidthList = []
        ERs=[]
        for iR in self.allChannels[channel]:
            ER=self._energies[iR]
            if Emin is not None and ER < Emin: continue
            if Emax is not None and ER > Emax: continue
            ERs.append(ER)
            reducedWidthList.append(numpy.sqrt(
                numpy.abs(self.allChannels[channel][iR])/2.0/self.penetrationFactorByChannel(channel, ER)))
        return reducedWidthList

    def getPorterThomasFitToWidths(self, Emin=0.0, Emax=None, verbose=False):
        """
        Perform a channel-by-channel fit of the histogram of widths to a Porter-Thomas distribution.

        :param verbose:
        :return:
        """
        try:
            import scipy.optimize
        except ImportError:
            print("WARNING: scipy.optimize not imported, Porter-Thomas analysis not done")
            return {}
        try: import numpy
        except ImportError:
            print("WARNING: numpy's histogram() and mean() functions not imported, Porter-Thomas analysis not done")
            return {}
        results={}

        for c in self.allChannels:
            if verbose: print()
            if verbose: print(c)
            if len(self.allChannels[c]) < 2:
                if verbose: print( '    not enough data' )
                results[c]={'dof':1.0,'ddof':0.0}
            else:
                reducedWidthList=self.getReducedWidths(c, Emin=Emin, Emax=Emax)
                aveWidth = numpy.mean(reducedWidthList)
                if not isinstance(reducedWidthList, numpy.ndarray):
                    norm=1
                else:
                    norm=len(reducedWidthList)
                if verbose:
                    print( '    widths', reducedWidthList )
                if norm == 1:
                    if verbose: print('    not enough widths to fit')
                    results[c] = {'dof': 1.0, 'ddof': 0.0}
                elif all([w==reducedWidthList[0] for w in reducedWidthList]):
                    if verbose: print( '    widths identical' )
                    results[c]={'dof':1.0,'ddof':0.0}
                elif aveWidth == 0.0 or numpy.isnan(aveWidth):
                    if verbose: print( '    widths are zero or nan' )
                    results[c]={'dof':1.0,'ddof':0.0}
                else:
                    if verbose: print( '    scaled widths', [w/aveWidth for w in reducedWidthList] )
                    hist, bin_edges=numpy.histogram([w/aveWidth for w in reducedWidthList])
                    #if verbose:  print '    hist',hist, any(hist)
                    #if verbose:  print '    bins', bin_edges, any(bin_edges)
                    if verbose: print( '    fitting:' )

                    if any(hist):

                        def PorterThomasDistribution(y, nu):
                            """
                            Porter-Thomas distribution
                            Froehner's eq. (277).
                            Really just a chi^2 distribution with nu degrees of freedom

                            :param y:
                            :param nu:
                            :return:
                            """
                            if verbose: print('        y', y)
                            if verbose: print('        norm', norm)
                            if verbose: print('        nu',nu)
                            gam = math.gamma(nu/2.0)
                            if nu < 2.0:
                                ycut=pow(1e4*gam,-1.0/(1.0-nu/2.0))
                                y[y < ycut]=ycut
                            return norm*numpy.exp(-y)*pow(y,nu/2.0-1.0)/gam

                        popt,pcov = scipy.optimize.curve_fit( PorterThomasDistribution, [(bin_edges[i+1]+bin_edges[i+2])/2 for i in range(len(bin_edges)-2)], hist[1:], bounds=[0.0,numpy.inf] ) # dump 1st bin, missing resonances there
                        if verbose: print('    popt',popt)
                        if verbose: print('    pcov',pcov)
                        results[c]={'dof':popt[0],'ddof':math.sqrt(pcov[0,0]),'bin_edges':bin_edges,'hist':hist}
        return results


#### Single-level Breit-Wigner ###
class SLBWcrossSection(RRBaseClass):
    """
    Resonance reconstructor for Singe-level Breit Wigner.

    Note, the resonances in the SLBW format each correspond to different "levels" and so do not
    interfere.  This is unlike all of the other resolved resonance formats.  So, while one resonance energy
    and one width go with one channel and all the channels associated with one resonance energy
    go together to make one reaction, the reactions for each resonance are added incoherently (no
    interference effects at all).
    """

    def __init__(self, BreitWignerForm, enableAngDists=False, verbose=True, **kw):
        parent = BreitWignerForm.ancestor
        RRBaseClass.__init__(self, BreitWignerForm, parent.domainMin, parent.domainMax, **kw)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f %s, reconstructing using Single-level Breit-Wigner" %
                    (self.lowerBound,self.upperBound,self.energyUnit))
        if enableAngDists: self.setResonanceParametersByChannel()

    def setResonanceParametersByChannel( self, multipleSScheme='ENDF', useReichMooreApproximation=False, Ein=None, warnOnly=False ):
        """
        Reorganize member data into channels (relies heavily on groundwork in sortLandJ).

        Note, unlike the getResonanceParametersByChannel() function in MLBW, RM or RML,
        the fact that different resonances are entirely different reactions means that the
        channelDicts have to have an additional layer of sorting that corresponds to the SLBW "level".
        """
        # make a unified table of resonance parameter data, convert to self.energyUnit if necessary
        self.nResonances = len( self.RR.resonanceParameters.table )
        params = ('energy','L','J','channelSpin','totalWidth','neutronWidth','captureWidth',
                'fissionWidth','fissionWidthB')
        EU = self.energyUnit
        units = (EU,'','','',EU,EU,EU,EU,EU)
        data = [self.RR.resonanceParameters.table.getColumn( quant,unit ) for quant,unit in zip(params,units) ]
        warnings = []
        for i in range(len(data)):
            if data[i] is None: data[i] = [0]*self.nResonances

        allowedSs = getAllowedTotalSpins( self.projectileSpin, self.targetSpin, useFactor2Trick=False )

        # Now make the dictionary of channels, grouped by independent level
        channelDicts = []
        for iR in range( self.nResonances ):
            ER = data[0][iR]
            L =  data[1][iR]
            J =  data[2][iR]
            S =  None # not given in SLBW, I just list it here for completeness
            GT = data[4][iR]
            GN = data[5][iR]
            GG = data[6][iR]
            GF = data[7][iR]
            GX = data[8][iR]
            JisAllowed = False
            allowedSByThisJl = []
            for SS in allowedSs:
                if abs(J) in getAllowedTotalSpins( L, SS, useFactor2Trick=False ):
                    JisAllowed = True
                    allowedSByThisJl.append( SS )
            if not JisAllowed:
                if warnOnly:
                    allowedSByThisJl.append( SS )
                    from fudge import warning
                    warnings.append(warning.InvalidAngularMomentaCombination(L, SS, abs(J), 'SLBW'))
                else: raise ValueError( "Invalid total angular momentum and/or orbital momentum: cannot couple up to J = "+str(J)+" with I = "+str(self.targetSpin)+", i = 0.5 and L = "+str(L) )
            gfact = (2.0*abs(J)+1)/((2*self.projectileSpin+1)*(2*self.targetSpin+1))
            channelDict=collections.OrderedDict()
            numAllowedSByThisJl = len( allowedSByThisJl )

            # Way recommended by NJOY:
            #   * One of valid S's get width, other(s) get zero width
            #   * Zero width channels left in sum to get potential scattering correct
            #   * Gets good (cs from U) to (cs from Fudge) comparison; good (cs from U) to (cs from BB)
            if multipleSScheme == 'NJOY':
                channelDict[ ChannelDesignator( L, abs(J), 'elastic', iR, allowedSByThisJl[0], gfact, None, None, 0.0, True, NEUTRONCHANNEL, False, False )] = [ ( ER, GN ) ]
                channelDict[ ChannelDesignator( L, abs(J), 'capture', iR, allowedSByThisJl[0], gfact, None, None, 0.0, False, GAMMACHANNEL, False, False )] = [ ( ER, GG ) ]
                if GF != 0.0: channelDict[ ChannelDesignator( L, abs(J), 'fission', iR, allowedSByThisJl[0], gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False )]     = [ ( ER, GF ) ]
                if GX != 0.0: channelDict[ ChannelDesignator(L, abs(J), 'competitive', iR, allowedSByThisJl[0], gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False)] = [(ER, GX)]
                for SS in allowedSByThisJl[1:]:
                    channelDict[ ChannelDesignator( L, abs(J), 'elastic', iR, SS, gfact, None, None, 0.0, True, NEUTRONCHANNEL, False, False )] = [ ( ER, 0.0 ) ]
                    channelDict[ ChannelDesignator( L, abs(J), 'capture', iR, SS, gfact, None, None, 0.0, False, GAMMACHANNEL, False, False )] = [ ( ER, 0.0 ) ]
                    if GF != 0.0: channelDict[ ChannelDesignator( L, abs(J), 'fission', iR, SS, gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False )]     = [ ( ER, 0.0 ) ]
                    if GX != 0.0: channelDict[ ChannelDesignator(L, abs(J), 'competitive', iR, SS, gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False)] = [(ER, 0.0)]

            # Way mandated by ENDF manual:
            #   * Width divided between channels with valid S's based on sign of J:
            #       ** if J>0, S is max possible value
            #       ** if J<0, S is min possible value
            #   * All channels left in sum to get potential scattering correct
            #   * Gets poor (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB)
            elif multipleSScheme == 'ENDF':
                if allowedSByThisJl:
                    allowedSByThisJl.sort()
                    if len(allowedSByThisJl)==1: SS=allowedSByThisJl[0]
                    else: SS=allowedSByThisJl[int(abs(J)==J)]
                    channelDict[ ChannelDesignator( L, abs(J), 'elastic', iR, SS, gfact, None, None, 0.0, True, NEUTRONCHANNEL, False, False )] = [ ( ER, GN/numAllowedSByThisJl ) ]
                    channelDict[ ChannelDesignator( L, abs(J), 'capture', iR, SS, gfact, None, None, 0.0, False, GAMMACHANNEL, False, False )] = [ ( ER, GG/numAllowedSByThisJl ) ]
                    if GF != 0.0: channelDict[ ChannelDesignator( L, abs(J), 'fission', iR, SS, gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False )]     = [ ( ER, GF/numAllowedSByThisJl ) ]
                    if GX != 0.0: channelDict[ ChannelDesignator(L, abs(J), 'competitive', iR, SS, gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False)] = [(ER, GX / numAllowedSByThisJl)]

            # Ignore problem:
            #   * Ignore spin of channels and the fact may be multiple valid spins
            #   * Gets best (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB)
            else :
                channelDict[ ChannelDesignator( L, abs(J), 'elastic', iR, S, gfact, None, None, 0.0, True, NEUTRONCHANNEL, False, False )] = [ ( ER, GN ) ]
                channelDict[ ChannelDesignator( L, abs(J), 'capture', iR, S, gfact, None, None, 0.0, False, GAMMACHANNEL, False, False )] = [ ( ER, GG ) ]
                if GF != 0.0: channelDict[ ChannelDesignator( L, abs(J), 'fission', iR, S, gfact, None, None, 0.0, False, FISSIONCHANNEL, False, False )]     = [ ( ER, GF ) ]
                if GX != 0.0: channelDict[ ChannelDesignator(L, abs(J), 'competitive', iR, S, gfact, None, None, 0.0, False, COMPETITIVECHANNEL, False, False)] = [(ER, GX)]
            channelDicts.append( channelDict )

        self.channels = channelDicts
        self.allChannels = []
        self.nChannels=sum(map(len,self.channels))
        if warnOnly: return warnings

    def getScatteringMatrixU( self, Ein, useTabulatedScatteringRadius=True, enableExtraCoulombPhase=False ):
        """
        Compute the scattering matrix

        Note, unlike the getScatteringMatrixU() function in other resonance classes,
        the fact that different resonances are entirely different reactions means that the
        channelDicts have to have an additional layer of sorting that corresponds to the SLBW "level".
        """
        U = collections.OrderedDict(
            [ ( iL, numpy.diag(len(cs)*[complex(1.0)]) ) for iL,cs in enumerate(self.channels) ] )
        rho = self.rho(Ein)

        # For calculating phi, default is to use tabulated scattering radius:
        if useTabulatedScatteringRadius:
            if self.RR.getScatteringRadius().isEnergyDependent():
                rhohat = self.RR.getScatteringRadius().getValueAs('10*fm', Ein) * self.k(Ein)
            else:
                rhohat = self.RR.getScatteringRadius().getValueAs('10*fm') * self.k(Ein)
        else: rhohat = self.rho(Ein)

        for iL,channels in enumerate(self.channels):
            ER = self._energies[iL]
            ERp = ER
            GamTot = 0.0
            eiphis = [ ]
            sqrtgams = [ ] # all are real
            for c in channels:
                if c.reaction == 'elastic':
                    Gam = ( channels[c][0][1] *
                        self.penetrationFactor( c.l, self.rho( Ein ) ) /
                        self.penetrationFactor( c.l, self.rho( abs( ER ) ) ) )
                    ERp += (
                            self.shiftFactor( c.l, self.rho( abs( ER ) ) ) -
                            self.shiftFactor( c.l, self.rho( Ein ) )
                        ) * channels[c][0][1] / (2.0 * self.penetrationFactor( c.l, self.rho( abs(ER) ) ) )
                    eiphis.append( numpy.exp( complex( 0.0, -self.phi( c.l, rhohat ) ) ) )
                else:
                    Gam = channels[c][0][1]
                    eiphis.append( complex( 1.0, 0.0 ) )
                sqrtgams.append( numpy.sqrt( Gam ) )
                GamTot += Gam
            for i1, c1 in enumerate(channels):
                denominator = complex( ERp - Ein, -GamTot/2.0 ) # c1 & c2 have same l
                for i2, c2 in enumerate(channels):
                    if c2.index != c1.index: continue
                    U[iL][i1][i2] += complex( 0.0, sqrtgams[i1]*sqrtgams[i2] ) / denominator
                    U[iL][i1][i2] *= eiphis[i1] * eiphis[i2]
        return U

    def getAngularDistribution(self, E, **kwargs ):

        raise NotImplementedError("Angular distributions cannot be safely reconstructed using Single-Level Breit-Wigner approximation")

    @blockwise
    def getCrossSection(self, E):
        captureSum = 0
        elasticSum = 0
        fissionSum = 0
        rho = self.rho(E)
        rhohat = self.rhohat(E)    # for calculating phase shift phi
        for L in self.Ls:
            l = L['L']
            phi = self.phi(l,rhohat)
            P = self.penetrationFactor(l,rho)
            S = 0.5*self.shiftFactor(l,rho)
            for J in L['Js']:
                gfactor = J['gfact']
                for spin in J['channelSpins']:
                    try:                    # use C extension if available
                        from . import getBreitWignerSums
                        capt, elas, fiss = getBreitWignerSums.getBreitWignerSums(E, spin['energy'],
                                spin['captureWidth'], spin['neutronWidth'], spin['fissionWidth'], spin['shiftFactor'],
                                P, S, phi, 'SingleLevel')
                        captureSum += gfactor * capt
                        elasticSum += gfactor * elas
                        fissionSum += gfactor * fiss
                    except ImportError:     # numpy version (memory-intensive)
                        dE = (E-(spin['energy']+spin['neutronWidth']*(spin['shiftFactor']-S)))
                        totalWidth = P*spin['neutronWidth'] + spin['captureWidth'] + spin['fissionWidth']
                        denominator = dE**2 + totalWidth**2 / 4
                        captureSum += gfactor * numpy.sum( ( P * spin['neutronWidth'] * spin['captureWidth'] ) / denominator , axis=1)
                        fissionSum += gfactor * numpy.sum( ( P * spin['neutronWidth'] * spin['fissionWidth'] ) / denominator , axis=1)
                        elasticSum += gfactor * numpy.sum( ( P*spin['neutronWidth'] * ( P*spin['neutronWidth']
                            - 2 * numpy.sin(phi)**2 * totalWidth + 4 * numpy.sin(phi)*numpy.cos(phi) * (E-spin['energy']) ) )
                            / denominator , axis=1)
            # numpy.sum(..., axis=1) returns row vector, so also convert first term to row vector:
            elasticSum += 4*(2*l+1)*numpy.sin(phi)[:,0]**2

        # get common factor 'beta' as a row vector:
        beta = numpy.pi / self.k(E)[:,0]**2
        capture = beta * captureSum
        elastic = beta * elasticSum
        fission = beta * fissionSum
        if self.removeNegatives:
            for reaction in (capture,elastic,fission):
                reaction[ reaction<=0 ] = 0
        total = elastic + capture + fission
        nonelastic = capture + fission
        result = {'total':total, 'elastic':elastic, 'capture':capture, 'fission':fission, 'nonelastic':nonelastic}
        return result

    def getScatteringLength(self, **kwargs):
        """
        Compute R' in b^2/2, should be close to AP.
        The potential scattering cross section sigPot = 4 Pi (R')^2, so we compute the potential scattering cross section at E=1e-5 eV.
        :param kwargs:
        :return:
        """

        E=numpy.array([ PQUModule.PQU(1e-5, 'eV').getValueAs(self.energyUnit) ])
        sigPotL = []
        for l in self.channels:
            for c in l:
                if not c.isElastic: continue
                sigPotL.append( ( ( 2.0 * c.l + 1 ) ) * 4.0 * numpy.pi * ( numpy.sin( self.phiByChannel(c,E) ) )**2 / self.k(E)**2 )
        sigPot = sum(sigPotL)
        return numpy.sqrt(sigPot/4.0/numpy.pi)[0]

#### Multi-level Breit-Wigner ###
class MLBWcrossSection(RRBaseClass):
    """
    Resonance reconstructor for Multi-level Breit-Wigner parameters.
    Only the elastic channel differs from SLBW.
    """

    def __init__(self, BreitWignerForm, enableAngDists=False, verbose=True, **kw):
        parent = BreitWignerForm.ancestor
        RRBaseClass.__init__(self, BreitWignerForm, parent.domainMin, parent.domainMax, **kw)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f %s, reconstructing using Multi-level Breit-Wigner" %
                    (self.lowerBound,self.upperBound,self.energyUnit))
        if enableAngDists: self.setResonanceParametersByChannel()

    def getChannelConstantsBc( self ):
        """
        For ENDF's MLBW, should be :math:`B_c = S_\\ell(|E_\\lambda|)`
        where :math:`\\ell` is the channel angular momentum and :math:`\\lambda` is the resonances index for the channel
        """
        raise NotImplementedError( "write me" )

    def getScatteringMatrixU( self, Ein, useTabulatedScatteringRadius=True, enableExtraCoulombPhase=False ):
        """
        Compute the scattering matrix.  We could have used the generic U function in the base class,
        but Froehner has "simplifications" that we took advantage of here (that and I don't know what the
        R matrix is exactly in the case of MLBW).

        :param Ein: numpy.array(type=float)
        :param boolean useTabulatedScatteringRadius:
        :param boolean enableExtraCoulombPhase:
        :return:
        """
        # Initialize lists and matrices
        sqrtGam = numpy.zeros( ( len(Ein), self.nChannels, self.nResonances ) )
        U = numpy.zeros( ( len(Ein), self.nChannels, self.nChannels ), dtype=complex )
        U[ :, numpy.identity( self.nChannels, dtype=bool ) ] = 1
        Gtots = self.nResonances * [ 0.0 ]
        ERp = list( self._energies )
        eiphis = self.getEiPhis(  Ein, useTabulatedScatteringRadius = useTabulatedScatteringRadius )

        # Precompute sqrt(Gamma) and sum up to get the total Gamma
        # Total Gamma may be different from what is given in the ENDF data file
        # Also, work out energy shifts
        for ic, c in enumerate( self.channels ):
            if c.reaction == 'elastic':
                for iR, G in self.channels[c].items():
                    ER = self._energies[iR]
                    Gam = ( G *
                        self.penetrationFactor( c.l, self.rho( Ein ) ) /
                        self.penetrationFactor( c.l, self.rho( abs( ER ) ) ) )
                    ERp[iR] += (
                            self.shiftFactor( c.l, self.rho( abs( ER ) ) ) -
                            self.shiftFactor( c.l, self.rho( Ein ) )
                        ) * G / 2.0 / self.penetrationFactor( c.l, self.rho( abs(ER) ) )
                    sqrtGam[ :, ic, iR ] = numpy.sqrt(Gam).T
                    Gtots[ iR ] += Gam
            else:
                for iR, G in self.channels[c].items():
                    ER = self._energies[iR]
                    sqrtGam[ :, ic, iR ] = numpy.sqrt(G).T
                    Gtots[ iR ] += G

        # If the width in the ENDF file is bigger than the sum of parts, there is a competitive width
        # otherwise, the ENDF value is wrong, use the one we get from the sum of parts
        for iR in range( self.nResonances ):
            mask = Gtots[ iR ] < self._widths[ iR ]
            Gtots[ iR ][ mask ] = self._widths[ iR ]

        # Compute the denominator
        DEN = ERp - Ein - 1j*numpy.array(Gtots)/2.0

        # Compute U itself.
        for i1 in range( self.nChannels ):
            for i2 in range( i1, self.nChannels ):
                U[:,i1,i2] += 1j * numpy.sum( sqrtGam[:,i1,:] * sqrtGam[:,i2,:] / DEN[:,:,0].T, axis=1 )
                U[:,i1,i2] *= (eiphis[i1] * eiphis[i2]).flatten()
                U[:,i2,i1] = U[:,i1,i2]     # symmetrize
        return U

    @blockwise
    def getCrossSection(self, E):
        captureSum = 0
        elasticSum = 0
        fissionSum = 0
        rho = self.rho(E)
        rhohat = self.rhohat(E)    # for phase shift phi
        for L in self.Ls:
            l = L['L']
            phi = self.phi(l,rhohat)
            P = self.penetrationFactor(l,rho)
            S = 0.5*self.shiftFactor(l,rho)
            for J in L['Js']:
                gfactor = J['gfact']
                for spin in J['channelSpins']:
                    try:                    # use C extension if available
                        from . import getBreitWignerSums
                        capt, elas, fiss = getBreitWignerSums.getBreitWignerSums(E, spin['energy'],
                                spin['captureWidth'], spin['neutronWidth'], spin['fissionWidth'], spin['shiftFactor'],
                                P, S, phi, 'MultiLevel')
                        captureSum += gfactor * capt
                        elasticSum += gfactor * elas
                        fissionSum += gfactor * fiss
                    except ImportError:     # numpy version (memory-intensive)
                        dE = (E-(spin['energy']+spin['neutronWidth']*(spin['shiftFactor']-S)))
                        totalWidth = P*spin['neutronWidth'] + spin['captureWidth'] + spin['fissionWidth']
                        denominator = dE**2 + totalWidth**2 / 4
                        commonFactor = P * spin['neutronWidth'] / denominator
                        captureSum += gfactor * numpy.sum( commonFactor * spin['captureWidth'] , axis=1)
                        fissionSum += gfactor * numpy.sum( commonFactor * spin['fissionWidth'] , axis=1)

                        # simple elastic method, Eq D.19 - D.21 in ENDF 102. This has numerical issues, however:
                        # U_nn = numpy.exp(-2*1j*phi[:,0]) * (1 + numpy.sum( 1j*P*spin['neutronWidth'] /
                        #     (spin['energy']-E-1j*totalWidth/2) , axis=1))
                        # elasticSum += gfactor * abs( (1-U_nn)**2 )

                        # Instead of the above, use the following from RECENT:
                        elasticTerm1 = numpy.sum( totalWidth/2 * commonFactor, axis=1 )
                        elasticTerm2 = numpy.sum( dE * commonFactor, axis=1 )
                        sin2ps = 2*numpy.sin(phi[:,0])*numpy.cos(phi[:,0])
                        sinps2 = 2*numpy.sin(phi[:,0])**2
                        elasticSum += gfactor * ((sinps2-elasticTerm1)**2 + (sin2ps+elasticTerm2)**2)
            if True: #addMissingGfactor:
                elasticSum += 2*numpy.sin(phi[:,0])**2 * 2*self.missingGfactor[l]

        # get common factor 'beta' as a row vector:
        beta = numpy.pi / self.k(E)[:,0]**2
        capture = beta * captureSum
        elastic = beta * elasticSum
        fission = beta * fissionSum
        if self.removeNegatives:
            for reaction in (capture,elastic,fission):
                reaction[ reaction<=0 ] = 0
        total = elastic + capture + fission
        nonelastic = capture + fission
        return {'total':total, 'elastic':elastic, 'capture':capture, 'fission':fission, 'nonelastic':nonelastic}


###### Reich_Moore and R-Matrix Limited ######

# some helper functions:

def getR_S(E, Eres, captureWidth, widths, penetrabilities, externals):
    """
    Both versions of Reich_Moore formalisms (LRF=3 and 7) rely on building the complex matrix 'R'.
    Here we break it up into the real component R and the imaginary component S,
    which represent symmetric and anti-symmetric scattering respectively

    matrix elements R[i,j] and S[i,j] =

    (summed over resonances) partialWidth[i]*partialWidth[j] * coefficient/(dE**2+captureWidth**2),

    for the ith/jth open channel. For S, the coefficient is 'dE', for R 'captureWidth'
    and partialWidth[i] = widths[i] * penetrabilities[i]

        ( widths[i] is a row vector of resonance widths, and penetrabilities[i] is a column vector
        of the penetrability for each incident energy )

    After computing R and S, invert using invertMatrices to find RI and SI such that (I+R+jS)*(I+RI+jSI) = I
    where I is the identity matrix, and j = sqrt(-1)

    Incident energy dependence appears in both in E and the penetrabilities
    """
    NE = len(E)         # number of energies
    dim = len(widths)   # dimension of matrix at each energy
    try:
        # use wrapped c version if available:
        from . import getScatteringMatrices
        # convert any constant width/penetrability data to array w/correct dimensions:
        nRes = len(Eres)
        for i1 in range(dim):
            if len(widths[i1]) != nRes:
                widths[i1] = numpy.ones((1,nRes)) * widths[i1][0,0]
            if len(penetrabilities[i1]) != NE:
                penetrabilities[i1] = numpy.ones((NE,1)) * penetrabilities[i1][0,0]
        R, S = getScatteringMatrices.getScatteringMatrices( E, Eres.copy(), captureWidth,
                numpy.array(widths), numpy.array(penetrabilities) )
    except ImportError:
        # can't import c version, so use numpy (memory-intensive)
        # define some common factors.
        # These are all 2d arrays with len(E) rows, len(Eres) columns:
        dE = Eres - E
        DEN = 1/(dE**2 + captureWidth**2)
        captOverDEN = captureWidth*DEN
        dEoverDEN = dE*DEN
        del dE, DEN

        R = numpy.zeros((NE,dim,dim))
        S = numpy.zeros((NE,dim,dim))
        for i in range(dim):
            for j in range(i,dim):
                width_ij = widths[i]*widths[j]
                p_ij = penetrabilities[i]*penetrabilities[j]
                R[:,i,j] = p_ij.T * numpy.sum( width_ij * captOverDEN , axis=1)
                S[:,i,j] = p_ij.T * numpy.sum( width_ij * dEoverDEN , axis=1)
                # symmetrize:
                R[:,j,i] = R[:,i,j]
                S[:,j,i] = S[:,i,j]
        del captOverDEN, dEoverDEN

    for idx in range(len(externals)):
        p_ii = penetrabilities[idx]**2
        S[:, idx, idx] += (p_ii * externals[idx][0]).flatten()
        R[:, idx, idx] += (p_ii * externals[idx][1]).flatten()

    return R, S


def invertMatrices( R, S ):
    """
    find RI and SI such that (I+R+jS)*(I+RI+jSI) = I,
    where I is the identity matrix, and j = sqrt(-1)

    for more info, see comments for subroutine FROBNS3 in recent.f
    """
    assert R.shape == S.shape
    NE,dim,dim = R.shape

    if dim==1:
        # only have neutron width (and 'eliminated' capture width).
        # Can obtain quicker solution by manually inverting:
        DET = (R+1)**2 + S**2
        SI = -S/DET
        RI = -(R*(R+1)+S**2)/DET
        return RI,SI
    else:
        # have competition and/or fission, must invert matrices:
        def dotproduct(*args):
            # helper function: dotproduct(A,B,C) == numpy.dot(A, numpy.dot(B,C))
            import functools
            def quickdot(a,b):
                # take dot product of two Nx(MxM) arrays (containing MxM matrices):
                result = numpy.zeros( a.shape )
                dim = a.shape[1]
                for i in range(dim):
                    for j in range(dim):
                        result[:,i,j] = numpy.sum( a[:,i,:] * b[:,:,j], axis=1 )
                return result
            return functools.reduce( quickdot, args )

        def vectorInv(arr):
            # invert all MxM matrices in Nx(MxM) array
            # FIXME are special cases for 2x2 and 3x3 matrices still useful? numpy.linalg.inv may end up faster
            dim = arr.shape[1]
            if dim==2:
                arrinv = numpy.zeros(arr.shape)
                det = (arr[:,0,0]*arr[:,1,1]-arr[:,0,1]*arr[:,1,0])
                arrinv[:,0,0] = arr[:,1,1]/det
                arrinv[:,0,1] = -arr[:,0,1]/det
                arrinv[:,1,0] = -arr[:,1,0]/det
                arrinv[:,1,1] = arr[:,0,0]/det
                return arrinv
            elif dim==3:
                arrinv = numpy.zeros(arr.shape)
                det = (arr[:,0,0]*(arr[:,1,1]*arr[:,2,2]-arr[:,1,2]*arr[:,2,1])
                    -arr[:,0,1]*(arr[:,1,0]*arr[:,2,2]-arr[:,1,2]*arr[:,2,0])
                    +arr[:,0,2]*(arr[:,1,0]*arr[:,2,1]-arr[:,1,1]*arr[:,2,0]) )
                arrinv[:,0,0] = (arr[:,1,1]*arr[:,2,2]-arr[:,1,2]*arr[:,2,1])/det
                arrinv[:,0,1] = (arr[:,0,2]*arr[:,2,1]-arr[:,1,0]*arr[:,2,2])/det
                arrinv[:,0,2] = (arr[:,0,1]*arr[:,1,2]-arr[:,0,2]*arr[:,1,1])/det
                arrinv[:,1,0] = (arr[:,1,2]*arr[:,2,0]-arr[:,1,0]*arr[:,2,2])/det
                arrinv[:,1,1] = (arr[:,0,0]*arr[:,2,2]-arr[:,0,2]*arr[:,2,0])/det
                arrinv[:,1,2] = (arr[:,0,2]*arr[:,1,0]-arr[:,0,0]*arr[:,1,2])/det
                arrinv[:,2,0] = (arr[:,1,0]*arr[:,2,1]-arr[:,1,1]*arr[:,2,0])/det
                arrinv[:,2,1] = (arr[:,0,1]*arr[:,2,0]-arr[:,0,0]*arr[:,2,1])/det
                arrinv[:,2,2] = (arr[:,0,0]*arr[:,1,1]-arr[:,0,1]*arr[:,1,0])/det
                return arrinv
            else:
                return numpy.linalg.inv(arr)

        identity = numpy.zeros( (NE, dim, dim) )
        for i in range(dim):
            identity[:,i,i] = 1

        invRplusI = vectorInv( R + identity )
        SRS = dotproduct( S, invRplusI, S )
        RI = dotproduct( -vectorInv( R+identity + SRS ), R+SRS )
        SI = dotproduct( -invRplusI, S, RI+identity )
        return RI, SI


#### R-Matrix Limited ####
class RMatrixLimitedcrossSection(RRBaseClass):
    """
    extended Reich_Moore (LRF=7 in ENDF)
    Here, resonances are sorted primarily by J: within each 'spin group', total J is conserved
    One or more competitive channels may be used in this case.
    Also, each resonance may have contributions from multiple l-waves
    """

    def __init__(self, RMatrixForm, enableAngDists=False, verbose=True, **kw):
        parent = RMatrixForm.ancestor
        RRBaseClass.__init__(self, RMatrixForm, parent.domainMin, parent.domainMax, **kw)
        self.verbose = verbose
        if self.verbose:
            print ("From %f to %f %s, reconstructing using RMatrix (approximation %s)" %
                    (self.lowerBound,self.upperBound, self.energyUnit,self.RR.approximation))

        # for RML, need to know what residual nuclei are produced by each channel:
        thresholds = []
        for rreac in self.RR.resonanceReactions:
            reaction = rreac.link.link
            # elastic, capture or competitive?
            if reaction == self.reactionSuite.getReaction('capture'): rreac.tag = 'capture'
            elif reaction == self.reactionSuite.getReaction('elastic'): rreac.tag = 'elastic'
            elif 'fission' in rreac.label: rreac.tag = rreac.label
            else: rreac.tag = 'competitive'

            if rreac.Q is not None:
                Q = rreac.Q.getConstantAs(self.energyUnit)
            else:
                try:
                    Q = reaction.getQ(self.energyUnit, final=False)
                except ValueError as ex:
                    if 'fission' in rreac.label:
                        # fission may have energy-dependent Q-value. Use minimum Q to compute Xi
                        Qform = reaction.outputChannel.Q.evaluated
                        Q = PQUModule.PQU(Qform.rangeMin, Qform.rangeUnit).getValueAs(self.energyUnit)
                    else:
                        raise ex
            # adjust Q value for residual in excited state:
            for particle in reaction.outputChannel.products:
                if( hasattr( particle, 'getLevelAsFloat' ) ) : Q -= particle.getLevelAsFloat( self.energyUnit )
            # Xi = reaction threshold in the lab frame.
            targetToProjectileMassRatio = self.targetMass_amu / self.projectileMass_amu
            Xi = - ( (targetToProjectileMassRatio+1) / targetToProjectileMassRatio ) * Q
            particles = reaction.outputChannel.products
            rreac.reactionInfo = { 'Xi':Xi, 'particles': particles }
            if Xi > 0:
                thresholds.append(Xi)

            if rreac.boundaryConditionValue is not None:
                raise NotImplementedError("Reconstruction with numeric boundary condition")

        self.haveExternalRMatrix = False
        for sg in self.RR.spinGroups:
            sg.energy = numpy.array( sg.resonanceParameters.table.getColumn('energy',self.energyUnit) )

            for channel in sg.channels:
                resonanceReaction = self.RR.resonanceReactions[ channel.resonanceReaction ]
                column = sg.resonanceParameters.table.columns[ channel.columnIndex ]
                column.tag = resonanceReaction.tag
                if channel.externalRMatrix is not None:
                    self.haveExternalRMatrix = True

        # for energy grid generation:
        energies, totalWidths = [],[]
        for sg in self.RR.spinGroups:
            energies.extend( sg.resonanceParameters.table.getColumn('energy',self.energyUnit) )
            widths = [sg.resonanceParameters.table.getColumn( col.name, self.energyUnit ) for col in
                    sg.resonanceParameters.table.columns if col.name != 'energy']
            totalWidths.extend( numpy.sum( widths, axis=0 ) )
        zipped = sorted(zip(energies,totalWidths))
        self._energies, self._widths = list( zip(*zipped) )
        self._thresholds = thresholds

        if enableAngDists:
            self.setResonanceParametersByChannel()

    def penetrationFactorByChannel(self, c, Ein):
        if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]: return 1.0
        rho = self.rho(Ein-c.Xi, c)
        if c.channelClass == CPCHANNEL:
            pA, pB = self.particlePairs[c].reactionInfo['particles']
            eta = self.eta(Ein-c.Xi, pA, pB)
            return getCoulombWavefunctions.coulombPenetrationFactor(c.l,rho,eta)
        return self.penetrationFactor(c.l, rho)

    def phiByChannel(self, c, Ein):
        if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]: return 0.0
        rho = self.rhohat(Ein-c.Xi, c)
        if c.channelClass == CPCHANNEL:
            pA, pB = self.particlePairs[c].reactionInfo['particles']
            eta = self.eta(Ein-c.Xi, pA, pB)
            return getCoulombWavefunctions.coulombPhi(c.l,rho,eta)
        return self.phi( c.l, self.rho(Ein,c) )

    def shiftFactorByChannel(self, c, Ein):
        if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]: return 0.0
        rho = self.rho(Ein-c.Xi, c)
        if c.channelClass == CPCHANNEL:
            pA, pB = self.particlePairs[c].reactionInfo['particles']
            eta = self.eta(Ein-c.Xi, pA, pB)
            return getCoulombWavefunctions.coulombShiftFactor(c.l,rho,eta)
        return self.shiftFactor(c.l, rho)

    def isElastic(self, reactionDesignator):
        rs = self.RR.rootAncestor
        return rs.getReaction(reactionDesignator) == rs.getReaction('elastic')

    def getLMax(self, maxLmax=10):
        """
        LMax is determined by the behavior of the Blatt-Biedenharn Zbar coefficients.  Inside each one, there is
        a Racah coefficient and a Clebsh-Gordon coefficient.  The CG coefficient looks like this::

                ( l1 l2  L )
                (  0  0  0 )

        So, this means two things.  First, the CG coeff (and hence Zbar) will be zero if l1+l2+L=odd.
        Second, the maximum value of L will be l1max+l2max.  Hence, Lmax=2*lmax.
        """
        return min( max( [ 2*c.l for c in self.channels ] ), maxLmax )

    def setResonanceParametersByChannel( self, multipleSScheme='ENDF', useReichMooreApproximation=False, Ein=None, warnOnly=False ):
        """
        Reorganize member data into channels
        :param multipleSScheme:  ignored, kept so has same signature as overridden function
        :param useReichMooreApproximation:  ignored, kept so has same signature as overridden function
        :return:
        """
        self.nResonances = sum( [ len( sg.resonanceParameters.table ) for sg in self.RR.spinGroups ] )
        self.JMax = max( [sg.spin for sg in self.RR.spinGroups] )
        self.channelConstantsBc = []
        self.particlePairs=collections.OrderedDict()
        self.scatteringRadiiTrue = None
        self.scatteringRadiiEffective = None
        # FIXME: overrides should not be necessary anymore, ResonanceReaction / Channel classes now have logic
        # for looking up scatteringRadius / hardSphereRadius using ancestry
        self.overrides=collections.OrderedDict()
        allowedSs = getAllowedTotalSpins( self.projectileSpin, self.targetSpin, useFactor2Trick=False )
        warnings = []

        # The unified table of resonance parameter data
        channelDict = collections.OrderedDict()

        def addOrUpdateDict(theDict, theKey, theValue):
            if theKey not in theDict:
                theDict[theKey] = {}
            if theValue is not None:
                theDict[theKey][theValue[0]] = theValue[1]

        # Pack the channels
        channelIndex=0
        for sg in self.RR.spinGroups:
            ERs = sg.resonanceParameters.table.getColumn( sg.resonanceParameters.table.columns[0].name, self.energyUnit )
            for channel in sg.channels:
                pp = self.RR.resonanceReactions[ channel.resonanceReaction ]
                if not (pp.eliminated or pp.isFission()):

                    Ia, Ib = self.getParticleSpins(pp.label)
                    Pia, Pib = self.getParticleParities(pp.label)

                    # Check validity of channel spin
                    if channel.channelSpin not in getAllowedTotalSpins( Ia, Ib, useFactor2Trick=False ):
                        if warnOnly:
                            from fudge import warning
                            warnings.append(warning.InvalidSpinCombination( Ia, Ib, channel.channelSpin, channel.label))
                        else:
                            raise ValueError( 'Invalid spin combination: cannot couple Ia = %s and Ib = %s up to S = %s for channel "%s"' %
                                    (str(Ia), str(Ib), str(channel.channelSpin), channel.label))
                    # Check the validity of total angular momenta
                    if sg.spin not in getAllowedTotalSpins( channel.L, channel.channelSpin, useFactor2Trick=False ):
                        if warnOnly:
                            from fudge import warning
                            warnings.append(warning.InvalidAngularMomentaCombination(channel.L, channel.channelSpin, sg.spin, channel.label))
                        else:
                            raise ValueError( 'Invalid spin combination: cannot couple L = %s and S = %s up to J = %s for channel "%s"' %
                                    (str(channel.L), str(channel.channelSpin), str(sg.spin), channel.label ) )
                    # Check the validity of the channel parity
                    if Pia is not None and Pib is not None:
                        #print(Pia,Pib)
                        expectedPi=Pia*Pib*pow(-1,channel.L)
                        if expectedPi != int(sg.parity):
                            if warnOnly:
                                from fudge import warning
                                warnings.append(
                                    warning.InvalidParity(channel.L, channel.channelSpin, sg.spin, sg.parity, channel.label))
                            else:
                                raise ValueError(
                                    'Invalid parity: cannot couple L = %s and S = %s up to J = %s for channel "%s"' %
                                    (str(channel.L), str(channel.channelSpin), str(sg.spin), channel.label))
                    gfact = (2.0 * abs(sg.spin) + 1) / ((2 * float(Ia) + 1) * (2 * float(Ib) + 1))

                else:
                    gfact = 1.0

                # Construct the channel designator
                particles = list(map(str, pp.reactionInfo['particles']))
                ps0 = particles[0]
                ps1 = particles[1] if len(particles) > 1 else None
                if 'photon' in [ps0,ps1]:  channelClass = GAMMACHANNEL
                elif 'gamma' in [ps0,ps1]: channelClass = GAMMACHANNEL
                elif 'n'  in [ps0,ps1]:    channelClass = NEUTRONCHANNEL
                elif pp.isFission():       channelClass = FISSIONCHANNEL
                else:                      channelClass = CPCHANNEL
                if ps1 is None:
                    assert channelClass == FISSIONCHANNEL
                c=ChannelDesignator(
                    index=channelIndex,
                    reaction=pp.label,
                    l=channel.L,
                    s=channel.channelSpin,
                    J=sg.spin,
                    gfact=gfact,
                    particleA=ps0,
                    particleB=ps1,
                    Xi=pp.reactionInfo['Xi'],
                    isElastic=self.isElastic(pp.label),
                    channelClass=channelClass,
                    useRelativistic=self.RR.relativisticKinematics,
                    eliminated=pp.eliminated )

                if channel.externalRMatrix:
                    c.externalRMatrix = channel.externalRMatrix
                channelIndex+=1

                # Don't put closed channels onto lists
                if Ein is not None and not numpy.all(c.is_open(Ein)):
                    if numpy.any(c.is_open(Ein)):
                        print("WARNING: Rethink your grid!  Your grid straddles a threshold at %s %s for channel %s, l=%i, s=%i, J=%i."
                                % (str(c.Xi), self.energyUnit, c.reaction, c.l, c.s, c.J))
                    continue

                # Save the particle pair, helps finding them later
                self.particlePairs[c]=pp

                # Log the channels
                addOrUpdateDict(channelDict, c, None)  # This catches potential scattering only channels
                if channel.scatteringRadius or channel.hardSphereRadius:
                    self.overrides[c] = channel
                else:
                    self.overrides[c] = None

                # Now do the actual resonance packing
                for iR, width in enumerate(sg.resonanceParameters.table[:,channel.columnIndex]):
                    addOrUpdateDict(channelDict, c, (self._energies.index(ERs[iR]), width))

        # Set up the channel->resonance mappings
        self.allChannels = channelDict
        self.channels = collections.OrderedDict()  # just the kept ones
        self.eliminatedChannels = collections.OrderedDict()
        for c in self.allChannels:
            if c.eliminated:    self.eliminatedChannels[ c ] = self.allChannels[ c ]
            else:               self.channels[ c ] = self.allChannels[ c ]
        self.nChannels = len( self.channels )
        self.identityMatrix = numpy.identity( self.nChannels, dtype=complex )

        if warnOnly: return warnings

    def resetResonanceParametersByChannel(self, multipleSScheme='ENDF', useReichMooreApproximation=False, Ein=None ):
        return self.setResonanceParametersByChannel( multipleSScheme=multipleSScheme, useReichMooreApproximation=useReichMooreApproximation, Ein=Ein )

    def getChannelConstantsBc(self):
        """
        For ENDF's Reich-Moore, should be :math:`B_c = -\\ell`
        where :math:`\\ell` is the channel angular momentum, but the ENDF manual says nothing about it.

        There is a per-channel parameter BCH that we will interpret as :math:`B_c`
        """
        Bc=[]
        for c in self.channels:
            pp=self.particlePairs[c]
            if hasattr(pp,'boundaryCondition'):  #FIXME: this variable disappeared!
                Bc.append(pp.boundaryCondition)
            else:
                Bc.append(0.0)
        return Bc

    def getAPByChannel(self, c, trueOrEffective='true'):
        if trueOrEffective=='true': return self.getChannelScatteringRadiiTrue()[c]
        elif trueOrEffective=='effective': return self.getChannelScatteringRadiiEffective()[c]
        else: raise ValueError( "trueOrEffective must be 'true' or 'effective'")

    def getChannelScatteringRadiiTrue(self):
        APT = collections.OrderedDict()
        for c in self.channels:
            if self.RR.calculateChannelRadius:
                a = 0.123 * self.targetMass_amu**( 1. / 3. ) + 0.08
            elif self.overrides[c] is not None and self.overrides[c].scatteringRadius is not None:
                a = self.overrides[c].scatteringRadius.getValueAs('10*fm')
            else:
                a = self.particlePairs[c].getScatteringRadius().getValueAs('10*fm')
            APT[c] = a
        return APT

    def getChannelScatteringRadiiEffective(self):
        APE = collections.OrderedDict()
        for c in self.channels:
            if self.overrides[c] is not None and self.overrides[c].hardSphereRadius is not None:
                a = self.overrides[c].hardSphereRadius.getValueAs('10*fm')
            else:
                a = self.particlePairs[c].getHardSphereRadius().getValueAs('10*fm')
            APE[c] = a
        return APE

    def rho(self, Ein, c):
        """
        Compute :math:`\\rho_c(E) = a_c * k_c(E)`, using the true scattering radius.
        ENDF uses it for calculating shift and penetrabilities.

        :param Ein: incident energy in the lab frame (shifted by a threshold, if appropriate)
        :param c: the channel designator
        :return: the value of rho (dimensionless)
        """
        if self.scatteringRadiiTrue is None: self.scatteringRadiiTrue = self.getChannelScatteringRadiiTrue()
        pA,pB = self.particlePairs[c].reactionInfo['particles']
        return self.k_competitive(Ein, pA, pB) * self.scatteringRadiiTrue[c] # dimensionless

    def rhohat(self, Ein, c):
        """
        Compute :math:`\\hat{\\rho}_c(E) = a_c * k_c(E)`, using the effective scattering radius
        ENDF uses it for calculating the phase (but in truth, there should be no effective scattering radius).
        (Caleb uses self.k below, but I think it should be self.k_competitive for the sake of consistency)

        :param Ein: incident energy in the lab frame (shifted by a threshold, if appropriate)
        :param c: the channel designator
        :return: the value of rho (dimensionless)
        """
        if self.scatteringRadiiEffective is None: self.scatteringRadiiEffective = self.getChannelScatteringRadiiEffective()
        pA,pB = self.particlePairs[c].reactionInfo['particles']
        return self.k_competitive(Ein, pA, pB) * self.scatteringRadiiEffective[c] # dimensionless

    def omega(self, eta, L):
        if L > 0: return sum( [numpy.arctan2(eta,n) for n in range(1,L+1)] )
        return numpy.zeros_like(eta)

    def getL0Matrix(self, Ein):
        """
        Get the :math:`L^0` matrix of Froehner, :math:`{\\bf L^0}_{cc'} = \\delta_{cc'} (L_c-B_c)`
        where :math:`L_c = S_c + i P_c`

        """
        if not hasattr(self,'channelConstantsBc') or self.channelConstantsBc == []:
            self.channelConstantsBc = self.getChannelConstantsBc()
        L0 = numpy.zeros( (len(Ein), self.nChannels, self.nChannels), dtype=complex )

        for ic,c in enumerate(self.channels):
            if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]:
                shift = 0.0
                penet = 1.0
            elif c.channelClass == CPCHANNEL:
                rho = self.rho(Ein-c.Xi, c)
                pA, pB = self.particlePairs[c].reactionInfo['particles']
                eta = self.eta(Ein-c.Xi, pA, pB)
                shift = getCoulombWavefunctions.coulombShiftFactor(c.l,rho,eta)
                penet = getCoulombWavefunctions.coulombPenetrationFactor(c.l,rho,eta)
            else:
                rho = self.rho(Ein-c.Xi, c)
                penet = self.penetrationFactor(c.l, rho)
                shift = self.shiftFactor(c.l, rho)

            # ENDF appears to set Bc = Sc most of the time.  When that's not the case,
            # Bc is set to something non zero, so in that case, compute L0 correctly
            if self.channelConstantsBc[ ic ] != 0.0:
                L0[:,ic,ic] = (shift + 1j*penet - self.channelConstantsBc[ ic ]).flatten()
            else:
                L0[:,ic,ic] = (1j*penet).flatten()
        return L0

    def getRMatrix(self, Ein):
        """
        The R matrix in the Reich-Moore approximation is
        :math:`R_{cc'}=\\sum_\\lambda{\\frac{\\gamma_{\\lambda c}\\gamma_{\\lambda c'}}{E_\\lambda - E - i\\Gamma_{\\lambda\\gamma}/2}}`

        :param Ein:
        :type Ein: numpy.array(type=float)
        :return:
        :rtype:
        """
        R = numpy.zeros( (len(Ein), self.nChannels, self.nChannels), dtype=complex )

        # Loop through all resonances
        for iR, ER in enumerate( self._energies ):

            # Extract the gamma width for the first gamma channel that has this resonance
            gamWidth = 0.0
            for cg in self.eliminatedChannels:
                if iR in self.eliminatedChannels[ cg ]:
                    gamWidth = self.eliminatedChannels[ cg ][ iR ]
                    if gamWidth != 0.0: break

            # Precompute the reduced widths
            redWidth = []
            for ic,c in enumerate(self.channels):
                if iR in self.channels[c] and self.channels[c][iR] != 0.0:
                    width = self.channels[c][iR]
                    shiftedER=numpy.array([abs(ER-c.Xi)])
                    rho = self.rho(shiftedER, c)
                    if c.channelClass == NEUTRONCHANNEL:
                        pen = self.penetrationFactor( c.l, rho )
                    elif c.channelClass == CPCHANNEL:
                        pA, pB = self.particlePairs[c].reactionInfo['particles']
                        eta = self.eta(shiftedER, pA, pB)
                        pen = getCoulombWavefunctions.coulombPenetrationFactor(c.l, rho, eta)
                        if numpy.isnan(pen):
                            if VERBOSE: print( iR, ER-c.Xi, c.l, rho, eta, width )
                            raise ValueError('pen=%s for channel %s and resonance #%i, but L0[%i,%i]=%s '%(str(pen),str(c),iR,ic,ic,str(self.getL0Matrix(Ein-c.Xi)[:,ic,ic])))
                    else: pen = 1.0
                    if pen != 0.0:
                        redWidth.append( numpy.copysign( numpy.sqrt(numpy.abs(width/2.0/pen)), width )[0] )
                    else:
                        if VERBOSE: print( iR, ER-c.Xi, c, rho, eta, width )
                        redWidth.append( 0.0 )
                else:   redWidth.append( 0.0 )

            # Loop through all channels to accumulate the R Matrix elements
            for ic1, c1 in enumerate( self.channels ):
                if not iR in self.channels[ c1 ]: continue
                for ic2, c2 in enumerate( self.channels ):
                    if ic2 > ic1: break     # matrix is symmetric
                    if not iR in self.channels[ c2 ]: continue
                    dR = (redWidth[ ic1 ] * redWidth[ ic2 ] / ( ER-Ein - 1j * gamWidth / 2.0 ))
                    if any( numpy.isnan(dR) ):
                        if VERBOSE: print( redWidth )
                        raise ValueError('nan in R-matrix for channels %s and %s '%(str(c1),str(c2)))
                    R[:, ic1, ic2] += dR[:,0]
                    if ic1 != ic2:
                        R[:, ic2, ic1] += dR[:,0]

        for ic1, c1 in enumerate(self.channels):
            if c1.externalRMatrix:
                external_real, external_imag = c1.externalRMatrix.evaluate(Ein.flatten())
                R[:, ic1, ic1] += external_real + 1j * external_imag

        return R

    def getScatteringMatrixT(self, Ein, useTabulatedScatteringRadius=True, enableExtraCoulombPhase=False):
        """

        :param Ein:
        :type Ein: numpy.array(type=float)
        :param useTabulatedScatteringRadius:
        :type useTabulatedScatteringRadius: bool
        :param enableExtraCoulombPhase:
        :type enableExtraCoulombPhase: bool
        :return:
        :rtype:
        """
        # Make sure all channels are open at all requested energies
        if not numpy.all([c.is_open(Ein) for c in self.channels]):
            raise ValueError("One or more channels are not open on all energies in requested energy grid")

        T = numpy.zeros( ( len(Ein), self.nChannels, self.nChannels ), dtype=complex )
        U = self.getScatteringMatrixU( Ein,
                                       useTabulatedScatteringRadius=useTabulatedScatteringRadius,
                                       enableExtraCoulombPhase=enableExtraCoulombPhase )
        for ic1, c1 in enumerate( self.channels ):
            pA, pB = self.particlePairs[c1].reactionInfo['particles']
            eta = self.eta(Ein-c1.Xi, pA, pB)
            wc  = self.omega(eta,c1.l)
            for ic2 in range( self.nChannels ):
                if enableExtraCoulombPhase:
                    if ic1==ic2: eTwoIWcDeltacc=numpy.exp(2j*wc).flatten()
                    else:        eTwoIWcDeltacc=0.0j
                    T[:, ic1, ic2] = eTwoIWcDeltacc  - U[:, ic1, ic2]
                else:
                    T[:, ic1, ic2] = float(ic1==ic2) - U[:, ic1, ic2]
        return T

    def getEiPhis(self, Ein, useTabulatedScatteringRadius=None, enableExtraCoulombPhase=False):
        """
        The phase factor for the collision matrix, :math:`\\Omega_c=e^{\\omega_c-\\varphi_c}`

        :param Ein:
        :type Ein: numpy.array(type=float)
        :param useTabulatedScatteringRadius:
        :type useTabulatedScatteringRadius: bool
        :param enableExtraCoulombPhase:
        :type enableExtraCoulombPhase: bool
        :return:
        :rtype:
        """
        # Precompute phase factor
        eiphis = []
        for ic, c in enumerate( self.channels ):
            if c.channelClass in [FISSIONCHANNEL, GAMMACHANNEL]:
                eiphis.append( numpy.ones( Ein.shape, dtype=complex ) )
            elif c.channelClass == CPCHANNEL:
                rho = self.rhohat(Ein-c.Xi, c)
                pA, pB = self.particlePairs[c].reactionInfo['particles']
                eta = self.eta(Ein-c.Xi, pA, pB)
                phic = getCoulombWavefunctions.coulombPhi( c.l, rho, eta ).flatten()
                if numpy.any(numpy.isnan(phic)):
                    raise ValueError( 'phi is NaN for channel %s: %s ' % (str(c),str( list( zip(Ein.flatten()-c.Xi,phic.flatten() ) ) ) ) )
                if enableExtraCoulombPhase:
                    wc = self.omega(eta,c.l).flatten()
                    eiphis.append( numpy.exp( 1j * wc.flatten() - 1j * phic.flatten() ) )
                else:
                    eiphis.append( numpy.exp( -1j * phic.flatten() ) )
            else:
                rho = self.rhohat(Ein-c.Xi, c)
                phic = self.phi( c.l, rho )
                eiphis.append( numpy.exp( -1j * phic.flatten() ) )
        return eiphis

    def k_competitive(self, Ex, pA, pB):
        """
        Calculate k for any 2-body output channel.

        Note that if pA and pB are target and neutron, this reduces to self.k(E) as defined above
        in the ResonanceReconstructionBaseClass

        :param Ex: incident energy - Xi (Xi is the lab frame reaction threshold)
        :param pA: particle A
        :param pB: particle B
        :return: k in b**-1/2
        """

        particleA = self.RR.PoPs[pA.pid]
        pA_mass = particleA.getMass( 'amu' )
        particleB = self.RR.PoPs[pB.pid]
        pB_mass = particleB.getMass( 'amu' )
        targ_mass = self.targetMass_amu

        return self.kConstant1 * numpy.sqrt( pA_mass * pB_mass / (pB_mass + pA_mass) / self.neutronMass_amu *
                targ_mass / (targ_mass + self.projectileMass_amu) ) * numpy.sqrt(abs(Ex))

    def eta(self, Ex, pA, pB):
        """
        eta, the  Sommerfeld parameter, given by :math:`\\eta = Z_A Z_B m_{red} \\alpha / ( \\hbar c k )`

        for competitive channels with 2 charged particles, parameter eta is used to find penetrability.
        eta is given in eq D.79 of ENDF manual and $e^2$ is the fine structure constant $\alpha$ and
        $m_{red}$ is the reduced mass.  eta is dimensionless.

        :param Ex: The incident energy
        :param pA: particle A
        :param pB: particle B
        :return:   the Sommerfeld parameter [dimensionless]
        """

        particleA = self.RR.PoPs[pA.pid]
        pA_z = chemicalElementMiscPoPsModule.ZAInfo( particleA )[0]

        particleB = self.RR.PoPs[pB.pid]
        pB_z = chemicalElementMiscPoPsModule.ZAInfo( particleB )[0]

        if pA_z * pB_z == 0: return numpy.zeros_like(Ex)

        pA_mass = particleA.getMass( 'amu' )
        pB_mass = particleB.getMass( 'amu' )
        mn = self.projectileMass_amu
        eSq_over_hbarSq = 3.4746085579272e-1    # ?, should be b**-1/2
        k_comp = self.k_competitive(Ex, pA, pB) # in b**-1/2
        eta = (eSq_over_hbarSq * pA_z * pB_z * (pA_mass/mn) * (pB_mass/(pB_mass+pA_mass)) / k_comp)
        return eta

    @blockwise
    def getCrossSection(self, E):
        """
        Return reconstructed cross sections at incident energy E
        :param E: incident energy (unit = self.energyUnit).  May be scalar or vector
        :return: dictionary with cross sections for each reaction.
            Cross sections are either scalar or vector, depending on type of parameter 'E'
        """
        elasticSum = 0
        absorbtionSum = 0
        fissionSum = 0
        nCompetitive = len( [ch for ch in self.RR.resonanceReactions if ch.tag=='competitive'] )
        competitiveSum = [0,] * nCompetitive

        eliminatedReaction = [rr for rr in self.RR.resonanceReactions if rr.eliminated]
        if self.RR.approximation is resolvedModule.RMatrix.Approximation.ReichMoore and len(eliminatedReaction) != 1:
            raise TypeError("Only 1 reaction can be eliminated in Reich-Moore approximation!")
        haveEliminated = False  # only set to true if at least one non-zero width is present

        for sg in self.RR.spinGroups:
            j = sg.spin
            gfact = (2*j+1.0)/((2*self.projectileSpin+1)*(2*self.targetSpin+1))

            widths = []
            penetrabilities = []
            phis = []    # phi is l-dependant
            Xis = []
            thresholdIndices = []
            chanIds = []
            eliminatedWidth = numpy.zeros( sg.resonanceParameters.table.nRows )
            externals = []

            for chan in sg.channels:
                resonanceReaction = self.RR.resonanceReactions[chan.resonanceReaction]
                chanWidths = numpy.array(sg.resonanceParameters.table.getColumn(chan.columnIndex, self.energyUnit))

                if resonanceReaction.eliminated:
                    eliminatedWidth += chanWidths / 2
                    if any(eliminatedWidth): haveEliminated = True
                    continue
                elif resonanceReaction.isFission():
                    # penetrability treated as 1 everywhere
                    reducedWidths = numpy.sqrt( abs(chanWidths) / 2 )
                    reducedWidths[ chanWidths<0 ] *= -1
                    widths.append( reducedWidths )
                    penetrabilities.append( numpy.array([[1]]) )
                    chanIds.append( chan.resonanceReaction )
                    continue
                else:
                    l = chan.L
                    pA,pB = resonanceReaction.reactionInfo['particles']
                    Xi = resonanceReaction.reactionInfo['Xi']
                    Xis.append( Xi )
                    chanIds.append( chan.resonanceReaction )

                    if Xi>0:
                        # find index of first incident energy above threshold:
                        import bisect
                        thresholdIndex = bisect.bisect( E, Xi )
                    else:
                        thresholdIndex = 0
                    thresholdIndices.append( thresholdIndex )

                    # channel radius:
                    def rho( Ex ):
                        if self.RR.calculateChannelRadius:
                            a = 0.123 * self.targetMass_amu**( 1. / 3. ) + 0.08
                        else:
                            a = chan.getScatteringRadius().getValueAs('10*fm')
                        return self.k_competitive(Ex, pA, pB) * a

                    # rhohat, for calculating phase shift phi:
                    rhohat = chan.getHardSphereRadius().getValueAs('10*fm') * self.k(E)
                    phinow = self.phi(l, rhohat)
                    phinow[ phinow/rhohat < 1e-6 ] = 0  # see subroutine 'facphi' in RECENT
                    phis.append( phinow )

                    # penetrability:
                    Ex1 = abs(sg.energy-Xi)         # evaluated at resonances
                    Ex2 = E-Xi; Ex2[Ex2<0] = 0      # evaluated at each incident energy
                    eta1 = self.eta(Ex1, pA, pB)
                    if resonanceReaction.eliminated or resonanceReaction.isFission():
                        penetrabilityAtResonances = numpy.ones_like( Ex1 )
                        penetrabilityAtEin = numpy.ones_like( Ex2 )
                    elif any(eta1 > 0):
                        # output channel has two charged particles, need Coulomb penetrability:
                        penetrabilityAtResonances = getCoulombWavefunctions.coulombPenetrationFactor(l, rho(abs(Ex1)), eta1)
                        if thresholdIndex > 0:  # threshold reaction
                            eta2 = numpy.zeros_like( Ex2 )
                            eta2[ thresholdIndex: ] = self.eta( Ex2[ thresholdIndex: ], pA, pB )

                            penetrabilityAtEin = numpy.zeros_like( Ex2 )
                            penetrabilityAtEin[ thresholdIndex: ] = numpy.sqrt(
                                    getCoulombWavefunctions.coulombPenetrationFactor(l, rho(Ex2[ thresholdIndex: ]),
                                    eta2[ thresholdIndex: ] ) )

                        else:
                            eta2 = self.eta(Ex2, pA, pB)
                            penetrabilityAtEin = numpy.sqrt( getCoulombWavefunctions.coulombPenetrationFactor(l,rho(Ex2), eta2) )
                    else:
                        # no Coulomb contribution:
                        penetrabilityAtResonances = self.penetrationFactor(l, rho(abs(Ex1)))
                        penetrabilityAtEin = numpy.sqrt( self.penetrationFactor(l, rho(Ex2)) )

                    reducedWidths = numpy.sqrt( abs(chanWidths) / (2 * penetrabilityAtResonances) )
                    reducedWidths[ chanWidths<0 ] *= -1
                    widths.append( reducedWidths )
                    penetrabilities.append( penetrabilityAtEin )

                if chan.externalRMatrix:
                    externals.append(chan.externalRMatrix.evaluate(E))
                else:
                    externals.append((E*0, E*0))

            # are there any threshold reactions (negative Q-values)?
            # If so, we must break up the calculation above/below each threshold
            if any(thresholdIndices):
                if thresholdIndices != sorted(thresholdIndices):
                    # sort by increasing threshold
                    (thresholdIndices, chanIds, widths, penetrabilities, externals) = list(zip(*sorted(
                        zip(thresholdIndices,chanIds,widths,penetrabilities,externals),
                        key=lambda tmp: tmp[0]  # sort only by threshold index (not by id or width)
                    )))
                assert thresholdIndices[-1] <= len(E)
                thresholdIndices = list(thresholdIndices) + [len(E)]

                RI = numpy.zeros((len(E), len(widths), len(widths)))
                SI = numpy.zeros(RI.shape)

                threshSet = sorted(set(thresholdIndices))
                for i1 in range(len(threshSet)-1):
                    low = threshSet[i1]
                    high = threshSet[i1+1]
                    nOpen = len([th for th in thresholdIndices if th-low <= 0])
                    RI_now, SI_now = invertMatrices(
                            *getR_S(E[low:high], sg.energy, eliminatedWidth, widths[:nOpen],
                                    [p[low:high] for p in penetrabilities[:nOpen]],
                                    [(extR[low:high], extI[low:high]) for (extR, extI) in externals[:nOpen]])
                            )
                    if nOpen==1:  # numpy insists on same number of dimensions for copy:
                        RI_now.shape = SI_now.shape = (len(RI_now),1,1)

                    RI[low:high, :nOpen, :nOpen] = RI_now
                    SI[low:high, :nOpen, :nOpen] = SI_now

            else:   # no threshold reactions
                RI,SI = invertMatrices(
                        *getR_S(E, sg.energy, eliminatedWidth, widths, penetrabilities, externals)
                        )

            # reconstruct:
            chanIds = numpy.array( chanIds )
            elasID = [tmp.label for tmp in self.RR.resonanceReactions if tmp.tag=='elastic'][0]
            competitiveIDs = [tmp.link.link.label
                              for tmp in self.RR.resonanceReactions if tmp.tag=='competitive']
            for chan in self.RR.resonanceReactions:
                if chan.eliminated:
                    continue
                # which matrix elements correspond to this channel?
                thisChanIds = numpy.where( chanIds==chan.label )[0]
                if len(thisChanIds)==0:
                    continue
                id1, id2 = min(thisChanIds), max(thisChanIds)+1
                assert len(thisChanIds) == id2-id1
                if chan.tag == 'elastic':
                    elas = 0; absorb = 0
                    # diagonal component
                    for i in range(id1,id2):
                        phi = phis[i] [:,0] # convert to row vector
                        sinsqr = numpy.sin(phi)**2
                        sincos = numpy.sin(phi)*numpy.cos(phi)
                        elas += ( sinsqr * (sinsqr+2*RI[:,i,i]) + sincos * (sincos + 2*SI[:,i,i]) )
                        absorb += -RI[:,i,i]
                    # add cross-terms:
                    elas += numpy.sum(numpy.sum( RI[:,id1:id2,id1:id2]**2 + SI[:,id1:id2,id1:id2]**2,
                            axis=1), axis=1)
                    absorb -= numpy.sum(numpy.sum( RI[:,id1:id2,id1:id2]**2 + SI[:,id1:id2,id1:id2]**2,
                            axis=1), axis=1)
                    elasticSum += gfact * 4 * elas
                    absorbtionSum += gfact * 4 * absorb

                else: # competitive or fission:
                    # sum of cross-terms between this channel and elastic:
                    elasticIds = numpy.where( chanIds==elasID )[0]
                    if len(elasticIds) == 0:
                        continue
                    el1, el2 = min(elasticIds), max(elasticIds)+1
                    comp = numpy.sum(numpy.sum(SI[:,id1:id2,el1:el2]**2 + RI[:,id1:id2,el1:el2]**2,
                            axis=1),axis=1)
                    if chan.isFission():
                        fissionSum += gfact * 4 * comp
                    else:
                        # may have more than one competitive: add to correct channel
                        competitiveSum[ competitiveIDs.index(chan.link.link.label)] += gfact * 4 * comp

        # get common factor 'beta' as a row vector:
        beta = numpy.pi / self.k(E)[:,0]**2
        elastic = beta * elasticSum
        competitive = [ beta * comp for comp in competitiveSum ]
        fission = beta * fissionSum
        if haveEliminated:
            absorbtion = beta * absorbtionSum
        else:
            absorbtion = sum(competitive) + fission # more numerically stable in absence of eliminated channel
        total = elastic + absorbtion

        retDict = {'total':total, 'elastic':elastic, 'fission':fission, 'nonelastic':absorbtion}
        competitiveNames = [ch.link.link.label for ch in self.RR.resonanceReactions if ch.tag == 'competitive']
        for key,val in zip(competitiveNames,competitive):
            retDict[ key ] = val
        if haveEliminated:
            retDict[eliminatedReaction[0].tag] = absorbtion - fission - sum(competitive)
        return retDict


##### unresolved resonance region. Only one formalism here: #####
class URRcrossSection(ResonanceReconstructionBaseClass):

    def __init__(self, URRform, verbose=True, **kw):
        """
        Initialize the URR management class

        :param URRform: The tabulatedWidths class to be reconstructed
        :param verbose: duh, verbosity flag
        :param kw: Python dict of keywords, currently none defined
        :return: None
        """

        ResonanceReconstructionBaseClass.__init__(self, URRform)
        self.verbose = verbose

        # energy boundaries for this region:
        self.URR = URRform
        parent = URRform.ancestor
        self.lowerBound = parent.domainMin
        self.upperBound = parent.domainMax
        if self.verbose:
            print ("Unresolved from %f to %f %s" % (self.lowerBound,self.upperBound,self.energyUnit))
        self.levelSpacings=None
        self.averageWidths=None
        self.DOFs=None
        self.egrid=None

    def getWidthsAndSpacings(self):
        """
        URR tables give us the average resonance parameters, the average resonance spacing
        and the number of degrees of freedom, assuming that the widths are distributed by a chi^2
        probability density function with the given number of degrees of freedom.

        This function sets several member data as dicts (indexed by L & J) of data and interpolable functions::

            - levelSpacings::

                ..math::
                    D_{L,J}(E)

            - averageWidths, assigned simple labels like 'elastic' and 'capture' for convenience::

                ..math::
                    \\Gamma_{L,J,c}

              there is a wrinkle here. The elastic 'width' is, according to the ENDF manual, a "reduced width" GN0.
              This is misleading since basic science literature (and Froehner in JEFF-18) defines the reduced width
              as::

                ..math::
                    \\Gamma_{L,J,n} = 2 P_L \\gamma_{L,J,c}^2

              However, ENDF's reduced neutron width is wildly different (see Eq. (D.99) in ENDF-102)::

                ..math::
                    \\Gamma{L,J,n} = \\Gamma_{L,J,n}^0 \\sqrt{E} P_L/\\rho

            - DOFs: degrees of freedom with same labels as averageWidths

            - reactionLabels: which reactions correspond to 'elastic', 'capture' and 'fission'

        :return: None
        """

        self.levelSpacings = collections.OrderedDict()  # FIXME why OrderedDict? Does order matter?
        self.averageWidths = collections.OrderedDict()
        self.DOFs = collections.OrderedDict()
        self.reactionLabels = collections.OrderedDict()

        reactionFinder = self.URR.findAttributeInAncestry('getReaction')
        for reaction in self.URR.resonanceReactions:
            found = False
            for simpleLabel in ('elastic','capture','fission'):
                if reaction.link.link is reactionFinder(simpleLabel):
                    self.reactionLabels[reaction.label] = simpleLabel
                    found = True
            if not found:
                raise NotImplementedError("Reaction '%s' not yet supported in unresolved region" % reaction.label)

        for L in self.URR.Ls:
            l = L.L
            for J in L.Js:
                j = J.J

                self.levelSpacings[(l,j)]=J.levelSpacing.data
                self.averageWidths[(l,j)]={}
                self.DOFs[(l,j)]={}

                # Set the widths
                for wid in J.widths:
                    if wid.resonanceReaction == 'competitive':
                        simpleLabel = wid.resonanceReaction
                    else:
                        simpleLabel = self.reactionLabels[ wid.resonanceReaction ]
                    self.averageWidths[(l,j)][simpleLabel] = wid.data  # FIXME: average neutron width is funky ENDF "reduced width", we should fix it!
                    self.DOFs[(l, j)][simpleLabel] = wid.degreesOfFreedom

    def rho(self, E):
        """
        The dimensionless parameter rho::

            ..math::
                \\rho = k(E) a

        We always calculate channel radius for unresolved region according to the ENDF manual

        :param E: the incident neutron energy
        :return: rho
        """
        a = 0.123 * self.targetMass_amu**( 1. / 3. ) + 0.08
        return self.k(E) * a

    def getFluctuationIntegrals(self, E, widths, DOF):
        """
        From subroutine GNRL3 in RECENT. If possible, this will be replaced
        with more basic approach (rather than using lookup table)... not finding
        appropriate equations right now

        Comments from GNRL3 sourcecode::

              Calculate unresolved resonance fluctuation function
              (original coding from AVERAGE4 by Mulki Bhat)
              (new weighting scheme from MC^2-II)

              This routine has been modified to calculate elastic, capture
              and fission fluctuation functions all during one call (as
              opposed to the original version which calculated each reaction
              separately).

              GNX, GGX, GFX and GXX are the widths for elastic, capture,
              fission and competition. MUN, MUF and MUX are the number of
              degrees of freedom for elastic, fission and competition (infinite
              number of degrees assumed for capture). RN, RC and RF are the
              calculated fluctuation integrals for elastic, capture and fission

              The number of degrees of freedom for each distribution (elastic,
              fission or competition) may be 1 to 4. If the number of degrees
              of freedom for any distribution is less than 1 or more than 4
              it will be treated as an infinite number of degrees of freedom
              (which infers that the widths are not distributed, but are rather
              all equal to the average value). This last case is simulated by
              defining an additional 10 point quadrature in which the weight
              for one point is 1.0 and the weight for all other points is zero.
              for the one point of weight 1.0 the average width will be used.

              :type E: Incident energy, may be scalar or numpy array
              :type widths: dictionary of numpy arrays,
                containing average widths for each channel (for a specific L/J combination)
              :type DOF: dictionary of degrees of freedom for specific L/J

              :returns tuple(numpy.array) with fluctuation integrals for elastic, capture and fission
        """
        XX = [[ 3.0013465e-03,1.3219203e-02,1.0004488e-03,1.3219203e-02,1.0e+0],
            [7.8592886e-02,7.2349624e-02,2.6197629e-02,7.2349624e-02,0.0e+0],
            [4.3282415e-01,1.9089473e-01,1.4427472e-01,1.9089473e-01,0.0e+0],
            [1.3345267e+00,3.9528842e-01,4.4484223e-01,3.9528842e-01,0.0e+0],
            [3.0481846e+00,7.4083443e-01,1.0160615e+00,7.4083443e-01,0.0e+0],
            [5.8263198e+00,1.3498293e+00,1.9421066e+00,1.3498293e+00,0.0e+0],
            [9.9452656e+00,2.5297983e+00,3.3150885e+00,2.5297983e+00,0.0e+0],
            [1.5782128e+01,5.2384894e+00,5.2607092e+00,5.2384894e+00,0.0e+0],
            [2.3996824e+01,1.3821772e+01,7.9989414e+00,1.3821772e+01,0.0e+0],
            [3.6216208e+01,7.5647525e+01,1.2072069e+01,7.5647525e+01,0.0e+0],]
        WW = [[1.1120413e-01,3.3773418e-02,3.3376214e-04,1.7623788e-03,1.0e+0],
            [2.3546798e-01,7.9932171e-02,1.8506108e-02,2.1517749e-02,0.0e+0],
            [2.8440987e-01,1.2835937e-01,1.2309946e-01,8.0979849e-02,0.0e+0],
            [2.2419127e-01,1.7652616e-01,2.9918923e-01,1.8797998e-01,0.0e+0],
            [1.0967668e-01,2.1347043e-01,3.3431475e-01,3.0156335e-01,0.0e+0],
            [3.0493789e-02,2.1154965e-01,1.7766657e-01,2.9616091e-01,0.0e+0],
            [4.2930874e-03,1.3365186e-01,4.2695894e-02,1.0775649e-01,0.0e+0],
            [2.5827047e-04,2.2630659e-02,4.0760575e-03,2.5171914e-03,0.0e+0],
            [4.9031965e-06,1.6313638e-05,1.1766115e-04,8.9630388e-10,0.0e+0],
            [1.4079206e-08,0.0000000e+00,5.0989546e-07,0.0000000e+00,0.0e+0],]

        if numpy.all(widths['capture']<=0):
            return 0,0,0
        shape = widths['elastic'].shape
        RE,RC,RF = numpy.zeros(shape), numpy.zeros(shape), numpy.zeros(shape)

        # We may store non-integer DOF, but we just end up converting to int:
        MUN = int( DOF['elastic'] )
        if MUN<1 or MUN>4: MUN=5
        MUF = int( DOF.get('fission',0) )
        if MUF<1 or MUF>5: MUF=5
        MUX=int( DOF.get('competitive',0) )
        if MUX<1 or MUX>4: MUX=5
        for j in range(10):
            xj = XX[j][MUN-1]
            wj = WW[j][MUN-1]
            if numpy.any(widths.get('fission',[])) and numpy.any(widths.get('competitive',[])):
                effj = widths['elastic']*xj+widths['capture']
                for k in range(10):
                    xk = XX[k][MUF-1]
                    wk = WW[k][MUF-1]
                    effjk = effj+widths['fission']*xk
                    for i in range(10):
                        xi = XX[i][MUX-1]
                        wi = WW[i][MUX-1]
                        factor = wi*wk*wj*xj/(effjk+widths['competitive']*xi)
                        RE += xj*factor
                        RC += factor
                        RF += xk*factor
            elif numpy.any(widths.get('fission',[])):
                effj = widths['elastic']*xj+widths['capture']
                for k in range(10):
                    xk = XX[k][MUF-1]
                    wk = WW[k][MUF-1]
                    factor = wk*wj*xj/(effj+widths['fission']*xk)
                    RE += xj*factor
                    RC += factor
                    RF += xk*factor
            elif numpy.any(widths.get('competitive',[])):
                effj = widths['elastic']*xj + widths['capture']
                for k in range(10):
                    xk = XX[k][MUX-1]
                    wk = WW[k][MUX-1]
                    factor = wk*wj*xj/(effj+widths['competitive']*xk)
                    RE += xj*factor
                    RC += factor
            else:   # only elastic and capture widths:
                factor = wj * xj / ( widths['elastic'] * xj + widths['capture'] )
                RE += xj*factor
                RC += factor
        for arr in (RE,RC,RF):
            arr[ widths['elastic']<=0 ] = 0
        return RE,RC,RF

    def generateEnergyGrid(self, interpolateWidths=True):
        """
        Get the energy grid for reconstructing unresolved resonances. Usually this is just the energy grid
        chosen by the evaluator for storing energy-dependent widths.

        ENDF has been revised (with VIII) to interpolate on widths rather than reconstructed cross sections.
        For discussion see Red's rant in D.E. Cullen "A Short History of ENDF/B Unresolved Resonance Parameters",
        LLNL Report LLNL-TR-461199, ENDF Report ENDF-369 (2010)).

        :param interpolateWidths: if True, interpolate the average widths
        :return: a tuple containing the grid and the interpolate flag
        """
        egrid = set()
        for L in self.URR.Ls:
            for J in L.Js:
                egrid.update( J.levelSpacing.data.domainGrid )  # FIXME what about units?
                for width in J.widths:
                    egrid.update( width.data.domainGrid )
        if min(egrid) > self.lowerBound or max(egrid) < self.upperBound:
            raise ValueError("unresolved region parameters don't span stated energy range:" +
                    " (%e - %e) vs (%e - %e)" % (min(egrid), max(egrid), self.lowerBound, self.upperBound) )
        egrid.update( [self.lowerBound, self.upperBound] )
        egrid = sorted(egrid)
        lidx = egrid.index(self.lowerBound)
        uidx = egrid.index(self.upperBound)
        egrid = numpy.array( egrid[lidx:uidx+1] )

        # Fill in any large gaps in original egrid:
        energyGridGaps = numpy.argwhere( egrid[1:] / egrid[:-1] > 3 ).flatten()
        if len( energyGridGaps ):
            e_set = set( egrid )
            fixGrid = numpy.array( [1.0, 1.25, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.2, 8.5] )
            for idx in energyGridGaps:
                lowend, highend = egrid[idx], egrid[idx+1]
                lowPow, highPow = int(math.floor(math.log10(lowend))), int(math.ceil(math.log10(highend)))
                lowend *= 1 + 1e-12                      # Kludge to handle cases where 3.9999999999999996e-05 should be 4e-5.
                for i in range( lowPow, highPow + 1 ):
                    points_to_add = [point for point in 10**i * fixGrid if lowend < point < highend]
                    e_set.update( points_to_add )
            interpolateWidths=True
            egrid = numpy.array( sorted( e_set ) )

        return egrid.tolist(), interpolateWidths

    @blockwise
    def getCrossSection(self, E):
        """
        Return unresolved region reconstructed cross sections at incident energy E
        :param E: incident energy (unit = self.energyUnit).  May be scalar or vector
        :return: dictionary with cross sections for each reaction.
            Cross sections are either scalar or vector, depending on type of parameter 'E'
        """

        if getattr(self, 'averageWidths', None) is None:  # do some pre-sorting for convenience
            self.getWidthsAndSpacings()

        # Integrate to get the average cross section, the ENDF way.  Assumes resonances are SLBW ones (ugh)
        captureSum = 0
        elasticSum = 0
        fissionSum = 0
        rho = self.rho(E)
        rhohat = self.URR.getScatteringRadius().getValueAs('10*fm', E) * self.k(E)
        for L in self.URR.Ls:
            l = L.L
            phi = self.phi(l,rhohat)
            for J in L.Js:
                j = J.J
                gfactor = (2.0*abs(j)+1)/((2*self.projectileSpin+1)*(2*self.targetSpin+1))

                VL = self.DOFs[(l,j)]['elastic'] * self.penetrationFactor(l,rho) / rho
                # Save all the widths in a new widths container to simplify coding in the fluctuating integral widget
                widths = {}
                for label in self.averageWidths[(l, j)]:
                    function = self.averageWidths[(l, j)][label]
                    widths[label] = numpy.array([function.evaluate(energy) for energy in E.flatten()]
                                                ).reshape(-1, 1)    # convert to column vector

                widths['elastic'] *= self.energyFactor * VL*numpy.sqrt(E)   # convert reduced width to 'regular' width

                RE, RC, RF = self.getFluctuationIntegrals(E, widths, self.DOFs[(l, j)])

                # common factor for all reactions:
                levelSpacing = numpy.array(
                    [self.levelSpacings[(l, j)].evaluate(energy) for energy in E.flatten()]).reshape(-1, 1)
                comfac = 2 * numpy.pi * gfactor * widths['elastic'] / levelSpacing

                elasticSum += ((RE * widths['elastic'] - 2 * numpy.sin(phi)**2) * comfac)[:, 0]
                captureSum += (RC * widths['capture'] * comfac)[:, 0]
                if 'fission' in widths:
                    fissionSum += (RF * widths['fission'] * comfac)[:, 0]

            elasticSum += 4*(2*l+1)*numpy.sin(phi[:, 0])**2

        # common factor 'beta' as row vector:
        beta = numpy.pi / self.k(E)[:, 0]**2
        capture = beta * captureSum
        elastic = beta * elasticSum
        fission = beta * fissionSum
        for reaction in (capture, elastic, fission):
            reaction[reaction <= 0] = 0
        total = elastic + capture + fission
        nonelastic = capture + fission

        xscs = {'total': total, 'elastic': elastic, 'capture': capture, 'fission': fission, 'nonelastic': nonelastic}
        return xscs

    def getTransmissionCoefficients(self, skipFission=True, method='weakCoupling'):
        """
        Compute the effective transmission coefficient in the resolved resonance regime.
        The method argument controls the schemes to do this.

        In each case, define the variable tau as::

            ..math::
                \\tau_c = \\pi\\overline{\\Gamma}_c/D

        These are the options:

            * `method=='weakCoupling'`:  This is the traditional weak coupling limit::

                ..math::
                    T_c=2\\tau_c

            * `method=='1stOrder'`: This is the same as the weak coupling limit.

            * `method=='2ndOrder'`: This is the 2nd order correction to the weak coupling limit::

                ..math::
                    T_c=2\\tau_c(1-\\tau_c)

            * `method=='sumRule'`: Use Moldauer's sum rule to extract the transmission coefficients directly from
              the RRR tables P.A. Moldauer Phys. Rev. Lett. 19, 1047-1048 (1967).  The equation is::

                ..math::
                    T_c=2\\tau_c\\left[\\sqrt{1+\\tau_c^2}-\\tau_c\\right]

            * `method=='opticalModel'`: Use Moldauer's 'optical model' transmission coefficient::

                ..math::
                    T_c=1-\\exp{(-2\\tau_c)}

            * `method`==`SPRT`: Use the

                ..math::
                    T_c=\\frac{2\\tau_c}{ (1 + \\tau_c/2)^2 + (P_c R^\\inf_c)^2 }

        :param skipFission: flag to skip fission, what else?
        :param method: method to compute the transmission coefficient
        :return: a dictionary of results, sorted by channel.  Tc entries are all PQU's or None if it cannot be computed
        """
        allowedMethods=["weakCoupling", '1stOrder', '2ndOrder', 'sumRule', 'opticalModel', 'SPRT']
        if method not in allowedMethods:
            raise ValueError('Transmission coefficient calculation method must be one of '+', '.join(allowedMethods))
        # Initialize the reduced width factors for the elastic channel
        redWidthFactor={}
        if not self.averageWidths: self.getWidthsAndSpacings()
        for lj in self.averageWidths:
            if lj[0] not in redWidthFactor:
                redWidthFactor[lj[0]]=XYs1dModule.XYs1d.createFromFunction(
                    XYs1dModule.XYs1d.defaultAxes(
                        labelsUnits={
                            XYs1dModule.yAxisIndex : ( 'gamma' , '' ),
                            XYs1dModule.xAxisIndex : ( 'Ex', 'eV' ) }),
                    self.averageWidths[lj]['elastic'].domain(),
                    lambda E,nope: math.sqrt(E)*self.penetrationFactor( lj[0], self.rho(E) )/self.rho(E),
                    {},
                    1e-6,
                    100)

        # Now compute the Tc's
        Tc={}
        channelClass={'capture':GAMMACHANNEL, 'elastic':NEUTRONCHANNEL, 'fission':FISSIONCHANNEL}
        for rxn in channelClass.keys():
            for lj in self.averageWidths:
                if rxn=='elastic':
                    tau=math.pi*redWidthFactor[lj[0]]*self.averageWidths[lj][rxn]/self.levelSpacings[lj]
                elif rxn=='fission':
                    if skipFission: continue
                    if rxn not in self.averageWidths[lj]: continue
                    tau=math.pi*self.averageWidths[lj][rxn]/self.levelSpacings[lj]
                else:
                    tau=math.pi*self.averageWidths[lj][rxn]/self.levelSpacings[lj]
                c=ChannelDesignator(lj[0], lj[1], rxn, len(Tc), int(2.0*abs(lj[0]-lj[1])), gfact=None,
                                    particleA=None, particleB=None, isElastic=(rxn=='elastic'),
                                    channelClass=channelClass[rxn], useRelativistic=False, eliminated=False)
                if method in ["weakCoupling", '1stOrder']:
                    Tc[c] = 2.0 * tau
                elif method=='2ndOrder':
                    Tc[c] = 2.0 * tau * (1.0 - tau)
                elif method=="opticalModel":
                    Tc[c] = tau.applyFunction(lambda x, par: 1.0 - math.exp(-2.0 * tau.evaluate(x)), None) #FIXME: "tau.evaluate(x)" should only be "x", but applyFunction() is broken as is exp()
                elif method == 'sumRule':
                    Tc[c] = tau.applyFunction(lambda x, par: 2.0 * tau.evaluate(x) * (math.sqrt(1.0 + tau.evaluate(x) * tau.evaluate(x)) - tau.evaluate(x)), None)
                else: #method=="SPRT"
                    Tc[c] = 2.0*tau/(1.0+tau/2.0)/(1.0+tau/2.0)

                # Fix axis label, still says "Gamma"
                Tc[c].axes[0].label='Tc(rxn=%s, L=%i, J=%s)'%(c.reaction, c.l, str(c.J))
        return Tc
