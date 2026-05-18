# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

.PHONY: all pdf crossSectionAdjustForHeatedTarget.pdf clean realclean

all: pdf

pdf: crossSectionAdjustForHeatedTarget.pdf

crossSectionAdjustForHeatedTarget.pdf: crossSectionAdjustForHeatedTarget.tex
	pdflatex crossSectionAdjustForHeatedTarget.tex
	pdflatex crossSectionAdjustForHeatedTarget.tex

clean:
	rm -f *.log *.aux *.dvi *.pdf

realclean: clean
