#!/bin/bash

set -e

jupyter nbconvert --to script 'generatePEOgromacs.ipynb'

Nseg=3
newline='Nseg = '$Nseg
linetoreplace=$(cat generatePEOgromacs.py | grep 'Nseg =')
sed -i '/'"$linetoreplace"'/c\'"$newline" generatePEOgromacs.py

python3 generatePEOgromacs.py
rm generatePEOgromacs.py

mv peo.itp ff/
