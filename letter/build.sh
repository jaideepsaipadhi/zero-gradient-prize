#!/bin/sh
# Build both versions of the letter. Run from letter/.
set -e
pdflatex -interaction=nonstopmode letter >/dev/null
bibtex letter
pdflatex -interaction=nonstopmode letter >/dev/null
pdflatex -interaction=nonstopmode letter
J='\def\FIVEP{}\input{letter}'
pdflatex -interaction=nonstopmode -jobname=letter_5p "$J" >/dev/null
bibtex letter_5p
pdflatex -interaction=nonstopmode -jobname=letter_5p "$J" >/dev/null
pdflatex -interaction=nonstopmode -jobname=letter_5p "$J"
