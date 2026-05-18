# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

SMR_PATH = ../../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../../..

ptwXYCFLAGS = -lnfut_utilities -lptwXY_utilities \
			  -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = groupTestAll1 integration integrationXY integrationWithWeight_x integrationWithWeight_sqrt_x integration2 runningIntegration \
	groupTwoFunction groupTwoFunction2 equalProbableBins groupFourFunctions

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

groupTestAll1: groupTestAll1.c
	$(CC) -o groupTestAll1 groupTestAll1.c $(localCFLAGS) -lm

integration: integration.c
	$(CC) -o integration integration.c $(localCFLAGS) -lm

integration2: integration2.c
	$(CC) -o integration2 integration2.c $(localCFLAGS) -lm

integrationXY: integrationXY.c
	$(CC) -o integrationXY integrationXY.c $(localCFLAGS) -lm

integrationWithWeight_x: integrationWithWeight_x.c
	$(CC) -o integrationWithWeight_x integrationWithWeight_x.c $(localCFLAGS) -lm

integrationWithWeight_sqrt_x: integrationWithWeight_sqrt_x.c
	$(CC) -o integrationWithWeight_sqrt_x integrationWithWeight_sqrt_x.c $(localCFLAGS) -lm

runningIntegration: runningIntegration.c
	$(CC) -o runningIntegration runningIntegration.c $(localCFLAGS) -lm

groupTwoFunction: groupTwoFunction.c
	$(CC) -o groupTwoFunction groupTwoFunction.c $(localCFLAGS) -lm

groupTwoFunction2: groupTwoFunction2.c
	$(CC) -o groupTwoFunction2 groupTwoFunction2.c $(localCFLAGS) -lm

equalProbableBins: equalProbableBins.c
	$(CC) -o equalProbableBins equalProbableBins.c $(localCFLAGS) -lm

groupFourFunctions: groupFourFunctions.c
	$(CC) -o groupFourFunctions groupFourFunctions.c $(localCFLAGS) -lm

check: default
	for f in $(TARGETS); do ./$$f $(CHECKOPTIONS); done
	./groupTwoFunction2 -a $(CHECKOPTIONS)

clean:
	rm -f *.dat v

realclean: clean
	rm -fr $(TARGETS) *.dSYM
