Processing GNDS files
=====================

The ``processProtare`` script (in fudge/bin/processProtare.py) is the main driver for processing
nuclear data evaluations. This script can be used to reconstruct resonances (if necessary),
heat reaction cross sections to various temperatures, generate group-averaged cross sections and transfer matrices,
and generate average product energy and momentum deposition data.

The name ``processProtare`` includes the short-hand 'ProTarE', which stands for "Projectile, Target and Evaluation".
Each GNDS ``reactionSuite`` contains this combination: an evaluation of all reactions involving a combination of
projectile and target.

Processed data generated during processProtare are stored internally in Fudge data structures, and once processing
is complete the resulting file is written back out as a GNDS/XML file.

Options for processing
----------------------

processProtare supports many command-line options. Perhaps the most useful is the '-h' or '--help' option::

        >python processProtare.py -h
        usage: processProtare.py [-h] [--pid PID] [--tid TID] [--noLazyParse]
                         [--exit EXIT] [--tag TAG] [-o OUTPUT]
                         [--writeConvertedUnits] [--energyUnit ENERGYUNIT]
                         [--CoulombPlusNuclearMuCutOff COULOMBPLUSNUCLEARMUCUTOFF]
                         [--doNotAddNuclearPlusInterference] [-t TEMPERATURES]
                         [--temperatureUnit TEMPERATUREUNIT]
                         [--legendreMax LEGENDREMAX] [-mg] [-up]
                         [--skipMultiGroupSums] [-mc] [--fluxFile FLUXFILE]
                         [--fluxID FLUXID] [--groupFile GROUPFILE] [-g GID]
                         [--formatVersion {1.10,2.0}] [--cullProcessedData]
                         [--threads THREADS] [-v]
                         [--reconAccuracy RECONACCURACY] [--restart]
                         [--preProcessLabel PREPROCESSLABEL]
                         [--prefixRecon PREFIXRECON]
                         [--prefixMuCutoff PREFIXMUCUTOFF]
                         [--prefixAEP PREFIXAEP] [--prefixHeated PREFIXHEATED]
                         [--prefixMultiGroup PREFIXMULTIGROUP]
                         [--prefixMonteCarlo PREFIXMONTECARLO]
                         [--prefixUpScatter PREFIXUPSCATTER]
                         mapOrProtareFileName [outputFile]

This script processes all data in a GNDS file as requested by input arguments. In addition to entering arguments on the comannd
line, arguments can also be read from an input file by specifying an input file on the command line. The character '@' must prefix
the input file's name. For example, to include the arguments in a file named pp.input execute this script as::

    processProtare.py -mc @pp.input -vvv -up eval.xml proc.xml

...

Sample usage
------------

processProtare takes a GNDS file as input. It also currently requires a LLNL-specific file called 'bdfls', which contains a database of
group boundaries and fluxes.

After translating an ENDF file into GNDS, you can do simple processing by:

::

        > python processProtare.py GNDS_file_name.xml -t 2.586e-8

This performs basic processing: Doppler-broadens cross sections to 300 K (2.586e-8 MeV/k) and computes average
energy and momentum deposited into reaction products. It writes a new processed file with the extension ".proc.xml" (i.e., if
the original file was named n-001_H_001.xml, this will output n-001_H_001.proc.xml).
This processing adds two new 'style' elements to the GNDS file: ``heated`` corresponding to the Doppler broadening, and ``averageProductData``
corresponding to the energy/momentum calculations.

Below are some examples of more advanced processing options.

Heating cross sections to multiple temperatures:

::

        > python processProtare.py GNDS_file_name.xml -t 300 -t 600 -t 900
                --temperatureUnit K 

Using eV instead of the default MeV for energies:

::

        > python processProtare.py GNDS_file_name.xml --energyUnit eV

Heating to two temperatures and generating cumulative distribution functions (CDFs) for faster sampling in Monte Carlo transport codes:

::

        > python processProtare.py GNDS_file_name.xml -t 300 -t 600
                --temperatureUnit K --energyUnit eV -mc
        
Generating multi-group cross sections and transfer matrices (requires that Merced be compiled):

::

        > python processProtare.py GNDS_file_name.xml --groupFile groups.xml --fluxFile fluxes.xml
                -t 300 -t 600 --temperatureUnit K --energyUnit eV -mg
                --fluxID LLNL_fid_1 -g n=LLNL_gid_4 -g photon=LLNL_gid_4

The 'LLNL_gid_4' group structure is defined in groups.xml, and 'LLNL_fid_1' is defined in fluxes.xml.
The resulting processed GNDS file will contain grouped cross sections and transfer matrices at both 300 and 600 K.


To-do list
----------

processProtare is not yet complete. Some functionality yet to be implemented:

- generating probability tables from unresolved resonance parameters
- reconstructing detailed angular distributions from resonance parameters
- processing covariances
- and likely other features
