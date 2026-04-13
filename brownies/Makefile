# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
PYTHON = python3

.PHONY: default clean realclean

default:

clean:
	find . -name "*.pyc" -exec rm -f {} \;

realclean: clean

FUDGETESTFILES = \
    BNL/utilities/test/test_c4.py \
    BNL/utilities/test/test_bins.py \
    BNL/inter/test/test_datatables.py \
    BNL/inter/test/test_metrics.py \
    BNL/inter/test/test_report.py \
    BNL/plot_evaluation/test/test_plotio.py \
    BNL/plot_evaluation/test/test_plotstyles.py \
    BNL/plot_evaluation/test/test_exfor2endf.py \
    BNL/plot_evaluation/test/test__init__.py \
    BNL/restools/test/test_level_analysis.py  \
    BNL/restools/test/test_level_density.py   \
    BNL/restools/test/test_level_generator.py

rebuild-test-data:
	cd BNL/inter/test/; $(PYTHON) rebuild_test_data.py

check: rebuild-test-data
	for testFile in $(FUDGETESTFILES); \
	do echo ; \
	echo =======================================================================; \
	echo \>\>\> TESTING $$testFile ; \
	echo =======================================================================; \
	echo ; \
	$(PYTHON) $$testFile; \
	done
