#!/bin/bash

set -e

jupyter nbconvert --to script 'generatePEGgromacs.ipynb'

Nseg=2
newline='Nseg = '$Nseg
linetoreplace=$(cat generatePEGgromacs.py | grep 'Nseg =')
sed -i '/'"$linetoreplace"'/c\'"$newline" generatePEGgromacs.py

python3 generatePEGgromacs.py
rm generatePEGgromacs.py
