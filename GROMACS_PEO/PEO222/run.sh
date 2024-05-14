#!/bin/bash

export GMX_MAXBACKUP=-1

gmx=gmx # /home/simon/Softwares/gromacs/build-gpu/bin/gmx

${gmx} grompp -f input/em.mdp -p topol.top -o em -pp em -po em 
${gmx} mdrun -deffnm em -v -nt 4 -pin on
cp em.gro conf.gro

${gmx} grompp -f input/nvt.mdp -p topol.top -o nvt -pp nvt -po nvt
${gmx} mdrun -deffnm nvt -v -nt 4 -pin on
cp nvt.gro conf.gro
