# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import time
import math
import numpy
import bisect

from fudge.core.math.fudgemath import RoundToSigFigs

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import gridded as griddedModule
from xData import XYs1d as XYs1dModule

from fudge.productData.distributions import angularEnergy as angularEnergyModule
from fudge.productData.distributions import energyAngular as energyAngularModule

from . import base as baseModule
from . import incoherentInelastic as incoherentInelasticModule

debug = False   # set True to disable multiprocessing
import multiprocessing
multiprocessing.set_start_method("fork")

class NeedShortCollisionTime( Exception ):

    pass

def process( self, style, energyMin, energyMax, temperature, kwargs ):
    """
    Convert S(alpha,beta) into double-differential cross section, store results under associated style
    Uses multiprocessing to parallelize over incident energy unless debug == True.

    @param self: incoherentInelastic.Form instance
    @param style: styles.Heated instance
    @param energyMin: minimum incident energy for double-differential cross section
    @param energyMax: maximum incident energy "  "
    @param temperature: desired temperature, in units of kwargs["temperatureUnit"]
    @param kwargs: processing options
    @return:
    """

    verbose = kwargs.get( 'verbosity', 0 )
    tolerance = 0.001

    # create energyAngular distribution by default
    order = kwargs.get("TNSL_order", "energy-angle")
    assert order in ("energy-angle", "angle-energy")

    primaryScatterer = self.getPrimaryScatterer()
    secondaryAtoms = [atom for atom in self.scatteringAtoms if atom is not primaryScatterer]

    if primaryScatterer.distinctScatteringKernel is not None:
        raise NotImplementedError("Processing with distinct scattering kernel present")
    if not isinstance(primaryScatterer.selfScatteringKernel.kernel, griddedModule.Gridded3d):
        raise NotImplementedError("Processing for self-scattering kernel of type %s" %
                type(primaryScatterer.selfScatteringKernel.kernel))

    gridded3d = primaryScatterer.selfScatteringKernel.kernel

    betas = numpy.array( gridded3d.axes[2].values )
    beta_interpolation = gridded3d.axes[2].interpolation

    alphas = numpy.array( gridded3d.axes[1].values )
    alpha_interpolation = gridded3d.axes[1].interpolation
    alpha_min = alphas[0]  # use power-law extrapolation down to alpha_min, then switch to constant S(alpha)
    if alpha_min > 1e-10:
        alpha_min = max(1e-10, 1e-5 * alpha_min)

    Sab_arr = gridded3d.array.constructArray( )

    temperatureUnit = gridded3d.axes[3].unit
    temperature = PQUModule.PQU( temperature, kwargs['temperatureUnit'] ).getValueAs( temperatureUnit )
    temperatures = gridded3d.axes[3].values
    temperatureIndex = bisect.bisect( temperatures, temperature )
    if( temperatureIndex == len( temperatures ) ) :
        temperatureIndex -= 1
    else :
        if( abs( temperature - temperatures[temperatureIndex-1] ) < abs( temperatures[temperatureIndex] - temperature ) ) :
            temperatureIndex -= 1

    # S_alpha_beta does not support temperature interpolation. Override temperature with the closest match in the file,
    # raise ValueError if the match isn't within 3%, and print a warning if the match is within 3% but not very close.
    newTemperature = temperatures[temperatureIndex]
    if not 0.97 < temperature / newTemperature < 1.03:
        raise ValueError(f"Requested TNSL temperature {temperature}{temperatureUnit} is not close to any temperature in the evaluation. " +
                         f"Closest available temperature = {newTemperature}{temperatureUnit}")
    if( not( numpy.allclose( temperatures[temperatureIndex], temperature, atol = 1e-12 ) ) ) :
        print( f"    WARNING: requested temperature {temperature}{temperatureUnit} not found, using {newTemperature}{temperatureUnit} instead" )

    temperature = temperatures[temperatureIndex]
    temperatureFactor = temperature / PQUModule.PQU( 0.0253, "eV/k" ).getValueAs( gridded3d.axes[3].unit )
    incidentEnergyUnit = kwargs['incidentEnergyUnit']
    kT = PQUModule.PQU( temperature, "k * %s" % gridded3d.axes[3].unit ).getValueAs( incidentEnergyUnit )     # k * temperature (i.e., an energy).

    Sab_arr = Sab_arr[temperatureIndex,:,:]                 # Table at desired temperature.

    liquid_small_beta_coef, alpha_min_liquid, log_S_0 = 0, 0, 0
    beta0 = betas.searchsorted(0)
    Sab_beta0 = Sab_arr[beta0]
    if Sab_beta0[0] > Sab_beta0[1]:
        # special case for liquids where S increases as alpha and beta go to 0
        log_S_0 = math.log(Sab_arr[beta0,0])
        alpha_min_liquid = alpha_min * 100
        liquid_small_beta_coef = (log_S_0 - math.log(Sab_arr[beta0+1,0])) * alphas[0]/betas[1]**2

    # Convert 2-d Sab array into a set of interpolation functions, one for each value of beta. Also store flags for when we need to switch to SCT approximation
    Sab = {}
    SCTflags = {}
    truncatedSab = []
    for betaIndex, beta in enumerate( betas ) :
        Sab_atBeta = Sab_arr[betaIndex,:]
        if( all( Sab_atBeta == 0 ) ) :
            SCTflags[beta] = 1.1 * alphas[-1]               # Always use SCT for this beta.
        elif( any( Sab_atBeta == 0 ) ) :
            nonzero = numpy.nonzero( Sab_atBeta )[0][0]     # First non-zero element.
            if nonzero > 0:
                SCTflags[beta] = alphas[nonzero]
                Sab_atBeta = Sab_atBeta[nonzero:]

            zeroes = Sab_atBeta == 0
            if( any( zeroes ) ) :
                truncatedSab.append( ( beta, sum( zeroes ) ) )
                Sab_atBeta[ zeroes ] = 1e-90                # In case log-lin or log-log interpolation

            Sab[beta] = XYs1dModule.XYs1d(list(zip(alphas[nonzero:], Sab_atBeta)), interpolation=alpha_interpolation, axes=axesModule.Axes(2))
        else :
            Sab[beta] = XYs1dModule.XYs1d(list(zip(alphas, Sab_arr[betaIndex,:])), interpolation=alpha_interpolation, axes=axesModule.Axes(2))

    if( truncatedSab ) :
        betasWithZeroes, nzeroes = list( zip( *truncatedSab ) )
        if( verbose > 1 ) : print( "    WARNING: %d zeroes (in %d betas) were changed to 1e-90" % ( sum( nzeroes ), len( betasWithZeroes ) ) )

    # Generate parameters for extrapolating to smaller alpha (unless we already switched to SCT).
    # Per Jesse Holmes' email, for small alpha S = C * alpha^N.  Compute C and N for each beta.
    alpha_extrapolation = {}
    for beta in Sab :
        if( beta in SCTflags ) : continue
        xys = Sab[beta]
        a1, s1 = list( map( math.log, xys[0] ) )
        a2, s2 = list( map( math.log, xys[1] ) )
        C = ( s1 * a2 - s2 * a1 ) / ( a2 - a1 )
        N = ( s1 - C ) / a1
        alpha_extrapolation[beta] = ( C, N )

    neutronMassAMU = self.findAttributeInAncestry('PoPs')['n'].mass.float('amu')
    principalAWR = primaryScatterer.mass.value / neutronMassAMU
    prefactor = primaryScatterer.numberPerMolecule * primaryScatterer.boundAtomCrossSection.getValueAs('b') / ( 2.0 * kT )

#<editor-fold desc="helper functions" defaultstate="collapsed">
    def getS_alpha_beta( alpha, beta ) :
        """
        Find Sab, either by interpolating on the grid or extrapolating below the alpha grid.
        Raises NeedShortCollisionTime if SCT approximation should be used.
        """

        if alpha <= 0: return 0

        if beta > betas[-1] or alpha > alphas[-1]: raise NeedShortCollisionTime()
        if self.asymmetric and beta < betas[0]: raise NeedShortCollisionTime()

        beta_index = bisect.bisect( betas, beta ) - 1

        beta1 = betas[beta_index]
        if beta1 in SCTflags and alpha < SCTflags[beta1]: raise NeedShortCollisionTime()
        if alpha < alphas[0]:
            if beta1 == 0 and alpha < alpha_min_liquid and liquid_small_beta_coef != 0:
                S1 = log_S_0 + math.log(alphas[0] / 1e-6) / 2 - liquid_small_beta_coef * beta ** 2 / 1e-6
                return math.exp(S1) # don't try to interpolate along beta

            C1, N1 = alpha_extrapolation[beta1]
            if alpha < alpha_min:
                S1 = math.exp(C1) * alpha_min**N1
            else :
                S1 = math.exp(C1) * alpha**N1
        else :
            S1 = Sab[beta1].evaluate(alpha)

        if beta1 == beta: return S1

        beta2 = betas[beta_index+1]
        if beta2 in SCTflags and alpha < SCTflags[beta2]: raise NeedShortCollisionTime()
        if alpha < alphas[0]:
            if beta2 == 0 and alpha < alpha_min_liquid and liquid_small_beta_coef != 0:
                S2 = log_S_0 + math.log(alphas[0] / 1e-6) / 2 - liquid_small_beta_coef * beta ** 2 / 1e-6
                return math.exp(S2)

            C2, N2 = alpha_extrapolation[beta2]
            if alpha < alpha_min:
                S2 = math.exp( C2 ) * alpha_min**N2
            else :
                S2 = math.exp( C2 ) * alpha**N2
        else :
            S2 = Sab[beta2].evaluate( alpha )

        xys = XYs1dModule.XYs1d([[beta1, S1], [beta2, S2]], interpolation=beta_interpolation, axes=axesModule.Axes(2))
        return xys.evaluate(beta)

    def evaluate(energy, mu, eprime):
        """
        Get double-differential cross section
        :param energy: incident energy
        :param mu: cos(scattering angle)
        :param eprime: outgoing energy
        :return: sigma(E,mu,Eprime)
        """

        beta = (eprime - energy) / kT
        alpha = (eprime + energy - 2.0 * mu * math.sqrt(energy * eprime)) / (principalAWR * kT)
        SCT = False

        alpha_lookup = alpha
        beta_lookup = beta
        if not self.asymmetric: beta_lookup = abs(beta_lookup)                 # S(-beta) = S(beta)
        if self.calculatedAtThermal:
            beta_lookup *= temperatureFactor
            alpha_lookup *= temperatureFactor

        try:
            Sab_term = getS_alpha_beta(alpha_lookup, beta_lookup) * math.exp(-0.5 * beta)
        except NeedShortCollisionTime:
            SCT = True
            Sab_term = primaryScatterer.shortCollisionTime(alpha, beta, temperature)

        result = prefactor * math.sqrt(eprime / energy) * Sab_term

        for atom in secondaryAtoms:
            if isinstance(atom.selfScatteringKernel.kernel, incoherentInelasticModule.SCTApproximation):
                # this atom only contributes when SCT is used
                if not SCT:
                    atom_Sab_term = 0
                else:
                    atom_Sab_term = atom.shortCollisionTime(alpha, beta, temperature)
                result += atom.numberPerMolecule * atom.boundAtomCrossSection.getValueAs('b') / (
                        2.0 * kT) * math.sqrt(eprime / energy) * atom_Sab_term

        return result

    def generate_eprime_grid(energy, significant_digits=None):
        """
        Generate initial outgoing energy grid for specified incident energy.
        Start with log-spaced grid of 4 points / decade from 1e-9 eV to 10 eV,
        then add points for each tabulated $\beta$.

        @param energy: incident energy
        @param significant_digits: round to desired sigfigs if specified
        @return:
        """
        offset = math.log10(PQUModule.PQU(1, 'eV').getValueAs(incidentEnergyUnit))
        eprimes = set(numpy.logspace(offset-9, offset+1, 10*4+1))
        eprimes.update([val for val in betas * kT + energy if val >= 0])
        if not self.asymmetric:
            eprimes.update([val for val in -betas * kT + energy if val >= 0])
        eprimes.add(0)
        eprimes = sorted(eprimes)
        if significant_digits is not None:
            eprimes = sorted(set(RoundToSigFigs(eprimes, sigfigs=significant_digits)))
        return eprimes

    def bisectionTest1D(spectrum, x1, y1, x2, y2, evaluator, epsilon, accuracy, level, max_depth):
        """
        Add points to 1-d function as needed to reach desired accuracy.
        Points are recursively added to the 'spectrum' list between (x1,y1) and (x2,y2), using the 'evaluator'
        function to determine the yMid value at xMid.
        Stop bisecting once either target accuracy or max_depth is reached.
        """

        if level > max_depth: return
        xMid = 0.5 * (x1 + x2)
        yMid = evaluator(xMid)
        yEstimate = 0.5 * (y1 + y2)
        if abs(yEstimate - yMid) <= accuracy * yMid:
            return
        if (xMid - x1) > (2 * epsilon * x1):
            bisectionTest1D(spectrum, x1, y1, xMid, yMid, evaluator, epsilon, accuracy, level + 1, max_depth)
        spectrum.append([xMid, yMid])
        if (x2 - xMid) > (2 * epsilon * xMid):
            bisectionTest1D(spectrum, xMid, yMid, x2, y2, evaluator, epsilon, accuracy, level + 1, max_depth)

    def bisectionTest2D(surface, x1, y1, x2, y2, evaluator, epsilon, accuracy, level, max_depth):
        """
        Similar to 1d version but adds an outgoing spectrum + spectrum integral at each value x
        """

        if level > max_depth: return
        xMid = RoundToSigFigs(0.5 * (x1 + x2), 12)
        if xMid in (x1, x2):
            return
        yMid, yMidSpectrum = evaluator(xMid)
        yEstimate = 0.5 * (y1 + y2)
        if abs(yEstimate - yMid) <= accuracy * yMid:
            return

        if (xMid - x1) > (2 * epsilon * x1):
            bisectionTest2D(surface, x1, y1, xMid, yMid, evaluator, epsilon, accuracy, level + 1, max_depth)
        surface.append([xMid, yMid, yMidSpectrum])
        if (x2 - xMid) > (2 * epsilon * xMid):
            bisectionTest2D(surface, xMid, yMid, x2, y2, evaluator, epsilon, accuracy, level + 1, max_depth)

    # <editor-fold desc="Evaluate angular-energy distribution" defaultstate="collapsed">
    def get_energy_spectrum_at_mu( energy, mu, eprimes, accuracy = 0.001, epsilon = 1e-6 ):
        """
        Generate appropriate outgoing energy grid at given incident energy / mu, and return spectrum.

        :param energy: incident energy
        :param mu: cos(scattering angle)
        :param eprimes: initial outgoing energy grid (grid depends on energy but not mu)
        :param accuracy, epsilon: parameters for refining interpolation grid
        :return: list [[eprime0, sigma0],[eprime1,sigma1]...]
        """

        initialGrid = [ [ ep, evaluate( energy, mu, ep ) ] for ep in eprimes ]
        yMax = max( [ v[1] for v in initialGrid ] )
        yCutoff = 10**math.floor( math.log10( yMax ) - 10 )
        if mu > 0.9999:
            yCutoff *= 1e-10 # allow more orders of magnitude due to strong forward-scattering peak
        while( initialGrid[-2][-1] < yCutoff ) : initialGrid.pop( )
        while( initialGrid[1][-1] < yCutoff ) : initialGrid.pop( 0 )

        spectrum = [ initialGrid.pop( 0 ) ]
        x1, y1 = spectrum[0]
        for x2, y2 in initialGrid:
            if x2 - x1 > 2 * epsilon * x1:
                bisectionTest1D(spectrum, x1, y1, x2, y2, evaluator=lambda eprime: evaluate(energy, mu, eprime),
                                epsilon=epsilon, accuracy=accuracy, level=0, max_depth=3)
            spectrum.append( [ x2, y2 ] )
            x1, y1 = x2, y2

        if( mu == 1 ):      # check for zeros near E==Eprime, remove if found
            badIndices = []
            for idx, (x,y) in enumerate( spectrum ):
                if y==0 and abs(x-energy)/energy < 1e-6:
                    badIndices.append(idx)
            for badIndex in badIndices[::-1]:
                spectrum.pop(badIndex)

        spectrum = angularEnergyModule.XYs1d(data=spectrum, outerDomainValue=mu, axes=axesModule.Axes(2))
        return spectrum.integrate(), spectrum

    def getAngularEnergyDistribution( energy, distributionAxes, accuracy = 0.001, epsilon = 1e-6, significant_digits = None ):
        """
        Generate appropriate mu grid at an incident energy, return P(mu) and outgoing energy distribution for each mu

        :param energy: incident energy
        :param accuracy, epsilon: parameters for refining interpolation grid
        :return list [[mu0, P(mu0|E), P(E'|E,mu0)], [mu1, P(mu1|E), P(E'|E,mu1)], ...]
        """
        mulist = [-1.0, -0.5, 0, 0.5, 0.75, 0.95, 0.99, 0.999999, 1.0]  # careful, mu=1 may diverge!
        eprimes = generate_eprime_grid(energy, significant_digits=significant_digits)

        initialGrid = []
        for mu in mulist:
            integral, spectrum = get_energy_spectrum_at_mu(energy, mu, eprimes)
            initialGrid.append([mu, integral, spectrum])

        surface = [initialGrid.pop(0)]
        x1, y1, y1Spec = surface[0]
        for x2, y2, y2Spec in initialGrid:
            bisectionTest2D(surface, x1, y1, x2, y2,
                            evaluator=lambda mu: get_energy_spectrum_at_mu(energy, mu, eprimes),
                            epsilon=epsilon, accuracy=accuracy, level=0, max_depth=5)
            surface.append([x2, y2, y2Spec])
            x1, y1 = x2, y2

        angularEnergy2d = angularEnergyModule.XYs2d( axes = distributionAxes, outerDomainValue = energy )
        for mu, integral, POfEPrime in surface:
            POfEPrime.axes = distributionAxes
            angularEnergy2d.append( POfEPrime )

        integral = angularEnergy2d.integrate()
        # FIXME thin after taking the integral but before normalizing?
        angularEnergy2d.scaleDependent(1/integral, insitu=True)
        return integral, angularEnergy2d
    # </editor-fold>

    # <editor-fold desc="Evaluate energy-angular distribution" defaultstate="collapsed">
    def get_mu_spectrum_at_eprime(energy, eprime, accuracy=0.001, epsilon=1e-6):
        """
        Return outgoing angular spectrum at given incident / outgoing energy.

        :param energy: incident energy
        :param eprime: outgoing energy
        :param accuracy, epsilon: parameters for refining interpolation grid
        :return: list [[mu0, sigma0],[mu1,sigma1]...]
        """

        mu_list = [-1, -0.5, 0, 0.5, 1]
        if 0.1 < eprime / energy < 10:
            # need more forward-peaked distribution near E == Eprime
            mu_list = [-1.0, -0.5, 0, 0.5, 0.75, 0.95, 0.99, 0.999999, 1.0]  # careful, mu=1 may diverge!

        initialGrid = [[mu, evaluate(energy, mu, eprime)] for mu in mu_list]

        spectrum = [initialGrid.pop(0)]
        x1, y1 = spectrum[0]
        for x2, y2 in initialGrid:
            if x2 - x1 > 2 * epsilon * x1:
                bisectionTest1D(spectrum, x1, y1, x2, y2, evaluator=lambda mu: evaluate(energy, mu, eprime),
                                epsilon=epsilon, accuracy=accuracy, level=0, max_depth=3)
            spectrum.append([x2, y2])
            x1, y1 = x2, y2

        spectrum = energyAngularModule.XYs1d(spectrum, outerDomainValue=eprime, axes=axesModule.Axes(2))
        return spectrum.integrate(), spectrum

    def getEnergyAngularDistribution(energy, distributionAxes, accuracy=0.01, epsilon=1e-6, significant_digits=None):
        """
        Generate appropriate outgoing energy grid at an incident energy, return P(eprime) and outgoing mu distribution for each eprime

        :param energy: incident energy
        :param accuracy, epsilon: parameters for refining interpolation grid
        :param epsilon:
        :param significant_digits:
        :return list [[eprime0, P(eprime0|E), P(mu|E,eprime0)], [eprime1, P(eprime1|E), P(mu'|E,eprime1)], ...]
        """
        eprimes = generate_eprime_grid(energy, significant_digits=significant_digits)

        initialGrid = []
        for eprime in eprimes:
            integral, spectrum = get_mu_spectrum_at_eprime(energy, eprime)
            initialGrid.append([eprime, integral, spectrum])

        yMax = max([v[1] for v in initialGrid])
        yCutoff = 10**math.floor(math.log10(yMax) - 15)
        while initialGrid[-2][1] < yCutoff: initialGrid.pop()
        while initialGrid[1][1] < yCutoff: initialGrid.pop(0)

        surface = [initialGrid.pop(0)]
        x1, y1, y1Spec = surface[0]
        for x2, y2, y2Spec in initialGrid:
            bisectionTest2D(surface, x1, y1, x2, y2, evaluator=lambda eprime: get_mu_spectrum_at_eprime(energy, eprime),
                            epsilon=epsilon, accuracy=0.01, level=0, max_depth=3)
            surface.append([x2, y2, y2Spec])
            x1, y1 = x2, y2

        energyAngular2d = energyAngularModule.XYs2d( axes = distributionAxes, outerDomainValue = energy )
        for eprime, integral, POfMu in surface:
            POfMu.axes = distributionAxes
            energyAngular2d.append( POfMu )

        integral = energyAngular2d.integrate()
        # FIXME thin after taking the integral but before normalizing?
        energyAngular2d.scaleDependent(1/integral, insitu=True)
        return integral, energyAngular2d
    # </editor-fold>

    def computeAtEnergy(energy, verbose=0, order="energy-angle"):
        """
        Compute outgoing angle/energy distribution and cross section at given incident energy.
        """
        if( debug and ( verbose > 1 ) ) :
            print("Compute TNSL at incident energy %g" % energy)

        if order == "energy-angle":
            evaluator = getEnergyAngularDistribution
            distributionAxes = energyAngularModule.defaultAxes( kwargs['incidentEnergyUnit'] )
        elif order == "angle-energy":
            evaluator = getAngularEnergyDistribution
            distributionAxes = angularEnergyModule.defaultAxes( kwargs['incidentEnergyUnit'] )
        else:
            raise AssertionError("Order must be 'energy-angle' or 'angle-energy'")

        return evaluator(energy, distributionAxes, significant_digits=11)

    def computeAtEnergies(energies):
        """
        Return (energy, crossSection, distribution) for each incident energy in list.
        Uses multiprocessing for faster results.
        """
        if debug:
            results = [(E, *computeAtEnergy(E, verbose, order=order)) for E in energies]
        else:
            # compute incident energies in parallel
            from multiprocessing import Process, Queue
            from fudge.defaults import numTasks
            processes = []
            resultsQueue = Queue()

            def enqueueResult( energy, verbose ):
                try:
                    results = (energy, *computeAtEnergy( energy, verbose, order=order ))
                    resultsQueue.put( results )
                    sys.exit( 0 )
                except Exception as err:
                    print(err)
                    sys.exit( 1 )

            nEnergies = len(energies)
            index = 0
            def addProcess( energy, processes ):
                process = Process( target=enqueueResult, args = (energy, verbose) )
                process.start()
                processes.append( [index, energy, process] )

            for pid in range( min(numTasks, nEnergies) ):
                addProcess( energies[index], processes )
                index += 1

            results = []
            while True:
                while not resultsQueue.empty():
                    results.append( resultsQueue.get() )

                nextProcesses = []
                for i2, energy, process in processes:
                    if process.is_alive():
                        nextProcesses.append( [i2, energy, process] )
                    else:
                        if process.exitcode != 0:
                            raise Exception("Error encountered in subprocess!")
                        if index < nEnergies:
                            addProcess( energies[index], nextProcesses )
                            index += 1
                            process.join()

                if len(results) == nEnergies:
                    break

                time.sleep(0.1)
                processes = nextProcesses

            results.sort()
        return results
#</editor-fold>

    # incident energy grid:
    energies = numpy.logspace( math.log10( energyMin ), math.log10( energyMax ), 180 ).tolist()

    if verbose:
        print("  Initial calculation at %d incident energies" % len(energies))
    results = computeAtEnergies(energies)

    # add incident energy points as needed (only for cross section):
    for idx in range(10):
        x,y,spectra = zip(*results)

        # FIXME following logic is nearly identical to checkInterpolation in reconstructResonances.py.
        # Move to common location?
        x = numpy.array(x)
        y = numpy.array(y)
        m = (y[2:]-y[:-2])/(x[2:]-x[:-2])
        delta_x = x[1:-1] - x[:-2]
        b = y[:-2]

        interpolated = numpy.zeros_like(x)
        interpolated[1:-1] = m*delta_x+b
        interpolated[0] = y[0]
        interpolated[-1] = y[-1]

        delt = interpolated / y
        mask = (delt>1+tolerance) + (delt<1-tolerance)

        badindices = numpy.arange(len(mask))[mask]
        if len(badindices)==0:
            break

        mask = numpy.zeros(len(x), dtype=bool)
        mask[badindices] = True
        midpoints = (x[:-1]+x[1:])/2
        energies_needed = sorted(set(list(midpoints[mask[:-1]]) + list(midpoints[mask[1:]])))
        if verbose:
            print("  Iteration %d, adding %d points" % (idx, len(energies_needed)))

        newResults = computeAtEnergies(energies_needed)
        results += newResults
        results.sort()

    # unpack results:
    crossSection = []

    if order == "energy-angle":
        energyAngularAxes = energyAngularModule.defaultAxes( kwargs['incidentEnergyUnit'] )
        energyAngular3d = energyAngularModule.XYs3d(axes=energyAngularAxes, interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)

        for energy, xsc, doublediff in results:
            crossSection.append( [energy, xsc] )
            if energy in energies:  # keep file size down by only adding energy grid to distribution
                energyAngular3d.append( doublediff )

        distribution = energyAngularModule.Form( style.label, self.productFrame, energyAngular3d )
    else:
        angularEnergyAxes = angularEnergyModule.defaultAxes( kwargs['incidentEnergyUnit'] )
        angularEnergy3d = angularEnergyModule.XYs3d(axes=angularEnergyAxes)

        for energy, xsc, doublediff in results:
            crossSection.append( [energy, xsc] )
            if energy in energies:  # keep file size down by only adding energy grid to distribution
                angularEnergy3d.append( doublediff )

        distribution = angularEnergyModule.Form( style.label, self.productFrame, angularEnergy3d )


    kwargs['multiplicity'] = XYs1dModule.XYs1d(data=[[energies[0], 1.0], [energies[-1], 1.0]], axes=axesModule.Axes(2))
    kwargs['energyAccuracy'] = kwargs['accuracy']
    kwargs['momentumAccuracy'] = kwargs['accuracy']
    kwargs['projectileMass'] = kwargs['neutronMass']
    # FIXME targetMass must be set for calculateAverageProductData even though it's not used for lab frame distributions
    kwargs['targetMass'] = None
    averageProductEnergy, averageProductMomentum = distribution.calculateAverageProductData( style, **kwargs )
    averageProductEnergy = averageProductEnergy[0]
    averageProductMomentum = averageProductMomentum[0]

    crossSection, averageProductEnergy, averageProductMomentum = baseModule.processedToForms( style.label, crossSection, averageProductEnergy, averageProductMomentum, kwargs )

    crossSection /= primaryScatterer.numberPerMolecule     # TNSL cross section is for the whole molecule, but we need it per atom

    return( crossSection, averageProductEnergy, averageProductMomentum, distribution )
