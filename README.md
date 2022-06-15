# Generate PEG topology for GROMACS and LAMMPS

<img align="right" width="30%" src="PEG.png">

### Description

Python scripts for generating a single all-atom PEG molecule with a number Nseg 
of monomer. The charmm36 force field is used, and the output is readable by GROMACS
and LAMMPS.

### How to (for GROMACS):

Generate the configuration files using generatePEGgromacs.ipynb, or simply by exebuting generatepeg.sh:

```
    sh generatepeg_gromacs.sh
```
You can vary the number of monomer by changing the value of Nseg. Then, minimise the energy of the PEG using gromacs:

```
    gmx grompp -f input/em.mdp -o em -pp em -po em
    gmx mdrun -v -deffnm em
    mv em.gro conf.gro
```
Eventually, relaxe the PEG molecule in the NVT ensemble:
```
    gmx grompp -f input/nvt.mdp -o nvt -pp nvt -po nvt
    gmx mdrun -v -deffnm nvt
```

### Output (GROMACS)

Pre-equilibrated topology files are given here : 
* [PEG200](PEG200/)
* [PEG600](PEG600/)
* [PEG1200](PEG1200/) 
 
This [video](https://www.youtube.com/watch?v=FkFdO58UdOA) has been made using 
the PEG molecule generated with this script.

### Contact

Feel free to contact me by email if you have inquiries. You can find contact 
details on my [personal page](https://simongravelle.github.io/).
