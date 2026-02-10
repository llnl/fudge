# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

subDirs = Doc Src TestSuite bin

.PHONY: default bin clean realclean check doc doSubDirs prefixDeploy

default: bin

bin:
	cd bin; $(MAKE)

prefixDeploy: bin
ifdef PREFIX
	cp bin/merced $(PREFIX)/bin
endif

clean:
	$(MAKE) doTarget=clean doSubDirs

realclean:
	$(MAKE) doTarget=realclean doSubDirs

check:

doc:
	cd Doc; $(MAKE) doc

doSubDirs:
	SAVED_PWD=`pwd`; \
	for subDir in $(subDirs); do cd $$SAVED_PWD/$$subDir; $(MAKE) $(doTarget); done
