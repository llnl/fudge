# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

.PHONY: default doc pdf clean realclean merced.pdf

default:

doc: merced.pdf

pdf: merced.pdf

merced.pdf: merced.tex assertion.tex doubleDiffFormula.tex \
  doubleDiffTable.tex fig4-1.tex fig4-2.tex fig4-3.tex fig4-4.tex \
  fig4-5.tex fig5-1.tex fig5-2.tex fig8-1.tex fig11-1.tex fig11-2.tex \
  figA-1.tex gammas.tex interpolate.tex isotropic.tex legendre-cm.tex \
  legendre-lab.tex relativity.tex summary.tex transfer.tex \
  two_body.tex uncorrelated.tex usage.tex
	pdflatex merced
	pdflatex merced
	pdflatex merced

clean:
	rm -f *.log *.aux *.toc *.out *.pdf

realclean: clean
