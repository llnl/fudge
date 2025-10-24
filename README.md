# FUDGE (For Updating Data and Generating Evaluations)

**FUDGE** is a nuclear data package built around Python that supports reading, viewing, modifying,
writing and processing nuclear data in the **GNDS** format.

The current release of **FUDGE** focuses on supporting version 2.0 of the Generalized Nuclear Database Structure (**GNDS**).
The **FUDGE** package includes tools to translate other nuclear data formats to and from **GNDS**,
plus tools for testing, visualizing, and processing **GNDS** data.

# Contents of this README:

- Getting Started
- Basic Use (translating ENDF files)
- Getting Help

## Getting Started:

Installing **FUDGE** requires Python (version 3.9 or higher) and NumPy (version 1.15 or higher).

Optional packages matplotlib and PyQT5 are also recommended to support plotting. 

### The following installation options are currently supported:

- Installation via `pip install`: This is ideal for use in virtual environments. The 
      following steps are recommended to `pip install` FUDGE:

    - Create a new virtual environment if desired. For example, using Anaconda:
 
          conda create --name fudge python=3.10 numpy matplotlib PyQT5
          conda activate fudge

      Or, using venv:

          python3 -m venv fudge
          source fudge/bin/activate
          pip install wheel
          pip install setuptools
          pip install numpy
          pip install matplotlib
          pip install PyQT5

      Some users have reported problems with installing PyQT5 with conda on Windows. It may be omitted from the
      above `conda create` command and, if required, be installed via `pip install` after activating the environment.

    - Install FUDGE:

          # Install latest available version:
          pip install git+https://github.com/LLNL/fudge.git

      Or,

          # Install a tagged release:
          pip install git+https://github.com/LLNL/fudge.git@6.11.0


- Installation by cloning the git repository and building with the unix `make` command: 
  This is the typical mode for active FUDGE maintenance and development.
  The following steps are recommended:

    - Ensure that NumPy (version 1.15 or greater) is installed

    - Clone FUDGE in the current directory: 
     
        git clone https://github.com/LLNL/fudge.git
            # or using SSH (requires creating a github account and registering an ssh key):
            git clone git@github.com:LLNL/fudge.git
     
    - Build FUDGE:

          cd fudge; make -s

    - [Optional] Set environment variables. For general use of the FUDGE package, 
          some changes should be made to your computer's environment. The following lines 
          make the required change. Note that <path_to_FUDGE> indicates the path to the 
          directory containing this README.md file.

      - on unix with the bash, ksh, or similar shell, put the following line in the 
          .bashrc or equivalent file:

            export PYTHONPATH=$PYTHONPATH:<path_to_FUDGE>

      - on csh, tcsh or similar shell into the appropriate shell file:

            setenv PYTHONPATH $PYTHONPATH:<path_to_FUDGE>

### Notes for Windows users:

We recommend installing and using FUDGE with Windows Subsystem for Linux (WSL).
When using WSL, be sure to install in a Linux-managed directory, e.g. /home/* instead of /mnt/c/*.

FUDGE can also be installed using the MinGW environment, but some portions of the code are optimized for
unix-style process management and perform best on WSL.
Compiling C extensions on MinGW may require an extra step:
if command 'cc.exe' is not available, locate a C compiler and set environment variable CC to that compiler
before running 'make' or 'pip install'. For example,
```
set CC=gcc
pip install ...
```

When installing with MinGW, the PYTHONPATH should be amended in the registry (see for example <http://www.support.tabs3.com/main/R10463.htm>)

Please let us know if you run into trouble installing or using FUDGE on Windows!

### Differences between FUDGE installed via `pip install` vs. via `make`:

When FUDGE is installed via `pip install`, several executable scripts are installed into 
the virtual environment 'bin' directory and can be used directly. The same behavior can be achieved
using the `make` installation method by adding <path_to_FUDGE>/bin and <path_to_FUDGE>/brownies/bin to the PATH.


## Basic Use:

After installing **FUDGE**, a user may be interested in translating an **ENDF-6** file into a **GNDS** file and vice versa.
Two tools for translating **GNDS** files to and from **ENDF-6** files are included in the *brownies/bin* directory.

The Python script **endf2gnds.py** translations an **ENDF-6** file to a **GNDS** file. 
This script writes a **GNDS** file containing the **reactionSuite** node and, if covariance data are 
present in the ENDF-6 file, its also writes a **GNDS** file containing the **covarianceSuite** node.

For example, the following command translating the **ENDF/B-VIII.0** file 
n-001_H_001.endf and generates n-001_H_001.xml and n-001_H_001-covar.xml:
```
python3 -m brownies.bin.endf2gnds.py n-001_H_001.endf n-001_H_001.xml

# or, if FUDGE was installed via pip:
endf2gnds.py n-001_H_001.endf n-001_H_001.xml
```

Some ENDF-6 files fail to translate unless the `--continuumSpectraFix` option is included when calling endf2gnds.
This option addresses a problem with unnormalizeable continuum photon spectra. It has no effect unless an evaluation
contains an unnormalizeable spectrum so can be used safely for all translations. Sample use:
```
endf2gnds.py <originalFile.endf> <newFile.xml> --continuumSpectraFix
```

The Python script **gnds2endf.py** translates a **GNDS** file to an **ENDF-6** file. This
script will look for a corresponding covariance file and, if present, add its data to the **ENDF-6** 
result. For example, the following translates the *n-001_H_001.xml* file generated in the **endf2gnds.py** command above 
back into an **ENDF-6** file:
```
python3 -m brownies.bin.gnds2endf.py n-001_H_001.xml

# or,
gnds2endf.py n-001_H_001.xml
```

**FUDGE** scripts use the **argparse.py** module to parse input arguments. To get *help* on a script, run the script
with the *-h* option. For example, to get help on the usage of the *bin/processProtare.py* script execute:
```
<path_to_FUDGE>/bin/processProtare.py -h
```

There are more scripts in the *bin* directory that a user may find useful. The *fudgeScripts.py* script provides an
overview of the other scripts included with FUDGE:

```
>fudgeScripts.py

Scripts in FUDGE:
    GNDSType.py                  - This script prints the GNDS type of each file listed.
    ZA_Info.py                   - For each argument entered, which must be an isotope name specified by either its ZA (1000 * Z + A) or its PoPs id, this script prints that argument, its ZA and PoPs
                                   id.
    ...
```

One can also run **FUDGE** interactively from the Python command prompt.
The following Python commands show how to read in the *n-001_H_001.xml* generated by the example above, print a list
of its reactions and plot each reaction's cross section on a single plot:
```
from fudge import GNDS_file     # FUDGE must be in the PYTHONPATH.
n_H1 = GNDS_file.read( 'n-001_H_001.xml' )
for reaction in n_H1.reactions: print( reaction )

crossSections = []
for reaction in n_H1.reactions:
    crossSection = reaction.crossSection.toPointwise_withLinearXYs( lowerEps = 1e-7 )
    crossSection.plotLegendKey = str( reaction )
    crossSections.append( crossSection )

crossSection.multiPlot( crossSections, rangeMin = 1e-5, xylog = 3,
        title = 'n + U230', xLabel = 'Neutron energy [eV]', yLabel = 'Cross section [b]' )
```

## Note regarding HDF5 support:

**FUDGE** imports the h5py module in some places to read and write HDF5 files. In some circumstances h5py fails to import
if MPI-related environment variables are set. This issue can cause the Python interpreter to crash. The best fix at the
moment is to ensure that environment variable `PMI_FD` is not set before running anything requiring h5py.

## Getting Help:

**FUDGE** documentation is currently undergoing an update but some documentation can be found by executing the following
**make** command in the **FUDGE** directory (note this only works after cloning FUDGE, not with pip install):
```
make -s docs
```
To view the documentation, open the file *doc/sphinx/_build/html/index.html* inside the **FUDGE** directory.

For other questions about the **FUDGE**, please feel free to contact the maintainers at
[mattoon1@llnl.gov](mailto:mattoon1@llnl.gov), [beck6@llnl.gov](mailto:beck6@llnl.gov),
or [gert1@llnl.gov](mailto:gert1@llnl.gov).

# License

**FUDGE** is distributed under the terms of the 3-clause BSD license.

All new contributions must be made under the 3-clause BSD license.

See LICENSE, COPYRIGHT and NOTICE for details.

SPDX-Licence-Identifier: BSD-3-Clause

This package includes several components:  
LLNL-CODE-683960	(FUDGE)

LLNL-CODE-770134	(numericalFunctions  
LLNL-CODE-771182	(statusMessageReporting)  
LLNL-CODE-725546	(Merced)

**FUDGE** is a product of the Nuclear Data and Theory Group at Lawrence Livermore National Laboratory (LLNL).

This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory
under Contract DE-AC52-07NA27344.
