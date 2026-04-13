# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

.PHONY: clean realclean build

default: build

build:
	$(PYTHON) -m pip install --quiet . --no-build-isolation
	find build -iname "_get*" -ipath "*build/lib*/_get*" -exec cp {} . \;

clean:
	rm -rf build
	rm -f *.pyc

realclean: clean
	cd test; $(MAKE) realclean
