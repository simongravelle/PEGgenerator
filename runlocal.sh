#!/bin/bash

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f input/em.mdp -o em -pp em -po em
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm em
mv em.gro conf.gro

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f input/nvt.mdp -o nvt -pp nvt -po nvt
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm nvt
