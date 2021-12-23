# Generate PEG topology for GROMACS

### Description

Python script for generating a single all-atom PEG molecule with a number Nseg of monomer. The charmm36 force field is used, and the output is readable by GROMACS.

### How to:

Generate the configuration files using generatePEGgromacs.ipynb, or simply by exebuting generatepeg.sh:

```
    sh generatepeg.sh
```
You can vary the number of polymer by changing the value of Nseg. Then, minimise the energy of the PEG using gromacs:

```
    gmx grompp -f innput/em.mdp -o em -pp em -po em
    gmx mdrun -v -deffnm em
```

### Output

This [video](https://www.youtube.com/watch?v=8ldIHP175TI) has been made with this script.

### Contact

Feel free to contact me by email if you have inquiries. You can find contact details on my [personal page](https://simongravelle.github.io/).
