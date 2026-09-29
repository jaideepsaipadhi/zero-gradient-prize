#!/bin/sh
# Build the JSC submission PDF (run from jsc/).
set -e
python3 flatten.py
pdflatex -interaction=nonstopmode main >/dev/null || true
bibtex main
pdflatex -interaction=nonstopmode main >/dev/null || true
pdflatex -interaction=nonstopmode main >/dev/null || true
grep -n '^!' main.log || echo "no LaTeX errors"
grep -i 'undefined' main.log || echo "no undefined references/citations"
grep -c 'Overfull' main.log || true
grep -o 'Output written.*' main.log
