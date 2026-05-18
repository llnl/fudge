######################
    FUDGE and GNDS
######################

:Release: |version|
:Date: |today|

Comments:

    * BRB: Should we move xData/documentation to LUPY?
    * BRB: In the isotopicAbundances stuff, should we change the attribute "evalution" to "library"?
        Oops, library is use in the iam node. Do we need a library specifier in the iam node.
    * BRB: Can we simplify the name IsotopicAbundancesByChemicalElement?
    * BRB: I used the node isotopicAbundancesByChemicalElement twice. Maybe one of them should be renamed.
    * BRB: methods like __init__ are not automatically include by sphinx. I think we should include some
        of them (see https://stackoverflow.com/questions/5599254/how-to-use-sphinxs-autodoc-to-document-a-classs-init-self-method).
    * BRB: In GNDS specs, I to not see PoPs/fissionFragmentData stuff.

Welcome to **FUDGE** (For Updating Data and Generating Evaluations), a nuclear data platform built with Python that permits
reading, viewing, modifying, writing and processing of nuclear data. FUDGE is a product of the Nuclear Data and
Theory group (NDT) at Lawrence Livermore National Laboratory (LLNL).

For installation instructions please see the file README.md or visit https://github.com/llnl/fudge.

This release of FUDGE expands support for processing and visualizing data in the Generalized Nuclear Data (GNDS) format.
GNDS is intended to modernize how nuclear data are stored and used, eventually replacing the legacy ENDF-6 format.
This version of FUDGE supports the latest versions of the GNDS standard, versions 2.0 and 2.1, as well as older versions '2.0.LLNL_4' and '1.10'.
The package includes tools to translate other formats to and from GNDS, plus tools
for testing, visualizing and processing GNDS-formatted data.
Python code within the *fudge* directory only support GNDS formatted data. Python code to support lecagy formats 
(e.g., ENDF-6 and ACE) reside in the *brownies* directory.

FUDGE requires Python 3 (version 3.7 and higher).
See the installation notes for more details on working with multiple versions of Python.

FUDGE is comprised of several sub-packages (e.g., **PoPs**, **xData**). 
Users are encouraged to check out the :doc:`sub-packages` documentation before continuing.
Also see the tutorials in :doc:`tutorial/index`. 

If you just want to see what a GNDS file looks like, visit the ENDF/B-VIII library pages. The library
is available in both ENDF-6 and `GNDS format <http://www.nndc.bnl.gov/endf/b8.0/gndsfiles.html>`_.

Contents
========

.. toctree::
    :maxdepth: 1

    downloading
    installation
    sub-packages
    gnds
    Tutorial <tutorial/index>
    externals
    reporting-bugs
    todo
    LUPY <LUPY/index.rst>
    pqu <pqu/modules.rst>
    xData <xData/index.rst>
    PoPs <PoPs/index.rst>
    Isotopic abundance <isotopicAbundances/index.rst>
    thanks

.. toctree::
   :hidden:
   
   The FUDGE Package <fudge/index>
   glossary 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`glossary`
* :ref:`search`
