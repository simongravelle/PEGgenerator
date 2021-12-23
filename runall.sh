gmx grompp -f input/em.mdp -o em -pp em -po em
gmx mdrun -v -deffnm em

mv em.gro conf.gro

gmx grompp -f input/nvt.mdp -o nvt -pp nvt -po nvt
gmx mdrun -v -deffnm nvt
