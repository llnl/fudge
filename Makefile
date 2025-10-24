# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL := /bin/bash

SUBMODULES = xData PoPs numericalFunctions brownies Merced crossSectionAdjustForHeatedTarget pqu

.PHONY: default build inplace bin rebuild-test-data clean realclean build \
        crossSectionAdjustForHeatedTarget numericalFunctions fudgeVersion

PYTHON = python3

default: all

all: fudgeVersion inplace2 bin

fudgeVersion:
	$(PYTHON) getFudgeVersion.py

build:
	$(PYTHON) setup.py --quiet build

inplace: inplace2
	@echo 'INFO: This target does not build bin'

inplace2: crossSectionAdjustForHeatedTarget numericalFunctions
	cd fudge/processing/resonances; $(MAKE) PYTHON=$(PYTHON)

# bin/Makefile calls both Merced/Makefile and fudge/processing/deterministic/upscatter/Makefile and copy executables
bin:
	cd bin; $(MAKE)

clean:
	rm -rf build
	find . -name "*.pyc" -exec rm -f {} \;
	find . -name "*.so" -exec rm -f {} \;
	cd crossSectionAdjustForHeatedTarget; $(MAKE) realclean
	cd numericalFunctions; $(MAKE) realclean
	cd bin; $(MAKE) realclean
	cd fudge/processing/deterministic/upscatter; $(MAKE) realclean
	cd PoPs; $(MAKE) realclean
	cd pqu; $(MAKE) realclean
	cd doc/sphinx; $(MAKE) realclean

realclean: clean
	rm -rf Out
	cd fudge; $(MAKE) realclean
	cd Merced; $(MAKE) realclean

pyclean:
	find . -name "*.pyc" -exec rm -f {} \;

tar:
	fileName=`git describe`; \
	rm -rf $$fileName; \
	mkdir $$fileName; \
	absolutePath=`cd $$fileName; pwd`; \
	git archive --format=tar HEAD | (cd $$fileName && tar -xf -); \
	git submodule foreach --recursive "git archive --prefix=\$$displaypath/ --format=tar HEAD | (cd $$absolutePath && tar -xf -)"; \
	find $$fileName -iname ".git*" -exec rm {} \; ; \
	tar -cf $${fileName}.tar $$fileName

docs:
	cd doc/sphinx; $(MAKE) $(MODE) html; cd ../..

rebuild-test-data:
	cd fudge/covariances/test; $(PYTHON) rebuild_test_data.py
	cd fudge/processing/resonances/test; $(PYTHON) rebuild_test_data.py

testPipInstall: rebuild-test-data check-pqu check-nf check-fudge
	$(PYTHON) -c "import fudge; print(f'FUDGE VERSION: {fudge.__version__}')"

check: rebuild-test-data check-pqu check-nf check-fudge

check-pqu:
	echo ===== check-pqu =====
	cd pqu/Check; $(PYTHON) check.py

check-nf:
	echo ===== check-nf =====
	cd numericalFunctions/nf_specialFunctions/Python/Test/UnitTesting/; $(MAKE) PYTHON=$(PYTHON) check

check-PoPs:
	echo ===== check-PoPs =====
	cd PoPs/Test/; $(MAKE) check

check-brownies: # BROKEN ?
	echo ===== check-brownies =====
	cd brownies; $(MAKE) PYTHON=$(PYTHON) check

check-heat: # BROKEN ... Need to change from bdfls and gnuplot to GNDS and new PyQt5 based plot methods
	echo ===== check-heat =====
	$(PYTHON) crossSectionAdjustForHeatedTarget/Python/Test/t.py

FUDGECOVTESTFILES =   fudge/covariances/test/test_base.py                           fudge/covariances/test/test_mixed.py \
                      fudge/covariances/test/test_summed.py                         fudge/covariances/test/test_covarianceSuite.py

check-cov: rebuild-test-data
	echo ===== check-cov =====
	for testFile in $(FUDGECOVTESTFILES); do echo ; \
		echo \>\>\> TESTING $$testFile ; \
		echo ; $(PYTHON) $$testFile; done

FUDGETESTFILES =    fudge/core/math/test/testFudgeMath.py                           fudge/core/math/test/test_linearAlgebra.py \
                    fudge/processing/resonances/test/test_getScatteringMatrices.py  brownies/legacy/endl/test/test_endlProject.py \
                    fudge/productData/distributions/test/__init__.py                fudge/reactionData/test/test_crossSection.py \
                    fudge/covariances/test/test_base.py                             fudge/covariances/test/test_mixed.py \
                    fudge/covariances/test/test_summed.py                           fudge/covariances/test/test_covarianceSuite.py \
                    fudge/processing/resonances/test/test_reconstructResonances.py  fudge/processing/resonances/test/test_makeUnresolvedProbabilityTables.py \
                    xData/test/test_XYs.py


check-merced:
	echo ===== check-merced =====
	mercedPath=$(shell pwd)/bin/merced; cd Merced/TestSuite; \
	num_threads=`python -c "import os; print(os.cpu_count())"` ;\
	num_cases=`find $(MERCEDTESTFOLDERS) -iname "in.*" | wc -l`;\
	currentCase=0;\
	for subdir in $(MERCEDTESTFOLDERS); do \
	  cd $$subdir;\
	  for file in in.*; do \
	    currentCase=$$((currentCase+1)) ;\
	    rm -f utfil $${file/in/new} $${file/in./}.info;\
	    $$mercedPath -num_threads $$num_threads $$file &> $${file/in./}.info;\
	    if [[ `cmp $${file/in/out} utfil` ]]; then \
          echo "  " $${subdir}/$$file output differs from baseline;\
        else\
          echo -ne "  $$currentCase of $$num_cases\r";\
	    fi;\
	    mv utfil $${file/in/new};\
	  done;\
	  cd ../;\
	done

MERCEDTESTFOLDERS = Compton ENDFLegendre ENDFdoubleDiff GeneralEvaporation Kalbach Legendre \
                    Legendre2Body Madland Maxwell Watt coherent doubleDifferential \
                    evaporation isotropic phaseSpace twoBody two_step uncorrelated weights

check-fudge: rebuild-test-data
	echo ===== check-fudge =====
	for testFile in $(FUDGETESTFILES); do echo ; \
		echo \>\>\> TESTING $$testFile ; \
		$(PYTHON) $$testFile; done

crossSectionAdjustForHeatedTarget:
	if [[ -d crossSectionAdjustForHeatedTarget/build ]]; then rm -rf crossSectionAdjustForHeatedTarget/build; fi
	cd crossSectionAdjustForHeatedTarget; $(PYTHON) setup.py --quiet build 
	find crossSectionAdjustForHeatedTarget/build -iname "*crossSectionAdjustForHeatedTarget*" \
	  -ipath "*build/lib*/crossSectionAdjustForHeatedTarget*/*crossSectionAdjustForHeatedTarget*" \
	  -exec cp {} crossSectionAdjustForHeatedTarget/lib \;

numericalFunctions:
	export PYTHONPATH=${PYTHONPATH}:`pwd`; cd numericalFunctions; $(PYTHON) setup.py --quiet build
	find numericalFunctions/build -ipath "numericalFunctions/build/lib*/numericalFunctions/*" ! -iname "__init__.py" \
	  -exec cp {} numericalFunctions/lib \;
