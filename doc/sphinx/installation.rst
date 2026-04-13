.. _installation:

Installation
============

The basic use of the FUDGE package requires the following:

* Git
* Python (versions 3.7 or higher) general-purpose programming language (with the numpy package included)
* C and C++ compilers for the translation of the compute-intensive modules.

While there are numerous options available for the installation of FUDGE, only the following are currently supported:

* Installation via `pip install` in a virtual environment or normal computing environmet that has the necessary installation privileges. This requires the following steps:

(i)  Ensure that NumPy (version 1.15 or higher) and wheels are installed in your Python environment

(ii) Run the command

    ::

        pip install git+https://github.com/LLNL/fudge.git

* Installation via Makefiles which is our typical development mode for active FUDGE code maintenance and improvements. The following steps are recommended:

(i) Ensure that NumPy (version 1.15 or higher) is installed

(ii) Use the following command to clone FUDGE to the current folder (users at LLNL may wish to clone the latest development version from CZ gitlab instead):

    ::

        git clone git+https://github.com/LLNL/fudge.git

(iii) Build FUDGE with the following command:

    ::

        cd fudge; make -s


[Optional] Extra python packages
--------------------------------
  
For more advanced use, Fudge depends on these additional packages:

* matplotlib for general plotting.
* PyQt5 and Qt v5 for interactive plotting. The underlying plotting is still done with matplotlib with the Qt toolkit and bindings providing the graphical
  user interface to interact and modify some of the common matplotlib plot options.
* scipy and h5py are required for more advance usage including procssing data libraries.


[Optional] Setting Environment Variables (not required after installing with pip)
---------------------------------------------------------------------------------

For general use of the fudge package, some changes should be made to your computer's environment. The following lines make the required change. 
Note that <path_to_fudge> indicates the path to root FUDGE folder:

* On unix with bash, ksh, etc, put this line in the .bashrc or equivalent:

    ::

        export PYTHONPATH=$PYTHONPATH:<path_to_fudge>

* On csh, tcsh or similar, do:
    
    ::
    
        setenv PYTHONPATH $PYTHONPATH:<path_to_fudge>

* On Windows, the environment variable should be added to the registry (see for example http://www.support.tabs3.com/main/R10463.htm)


Troubleshooting
---------------

Users on OS X may have trouble compiling Merced (the C++ code responsible for generating transfer matrices
for deterministic transport). Merced uses OpenMP threading, which is currently not supported by the default
clang compiler on OS X.
This problem is associated with an error message like the following:

::

    error: unsupported option '-fopenmp'

You may safely ignore this error message unless you wish to use the deterministic processing capabilities in Fudge.
If you encounter this error and do need to use Merced, you will need to install a compiler with OpenMP support,
then specify which compiler to use when building Merced.

For example, the g++-mp-4.9 compiler (available from Macports) supports openmp

::

   sudo port install gcc-9
   # then
   make merced CXX=g++-mp-4.9

Windows users may also have trouble with installing Merced. The solution here appears to be to install
VisualStudio, and specify the 'CL' compiler while building Merced:

::

   make merced CXX=CL

If you encounter other trouble with installing, building extensions or setting environment
variables, please let us know (email: mattoon1@llnl.gov).

