

import numpy as np
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt 

from utilities import check_input, convert_coordinates, load_molecule_templates, \
    build_peg, analyze_atom_types, center_molecule, check_box_boundaries, \
    validate_atoms, estimate_molar_mass, build_cc_bonds, build_co_bonds, \
    build_ch_bonds, build_oh_bonds, combine_bonds, \
    check_missing_ch_bonds, calculate_angles, check_angle_counts, calculate_dihedrals, \
    check_dihedral_counts, validate_dihedrals, combine_angles, combine_dihedrals, \
    write_lammps_data, check_dihedral_atom_ids, write_gro

### number of segments for the PEG molecules

nseg = 1
spacing_nm = 0.2 # initial distance between 2 monomers in nm
spacing_A = spacing_nm*10
### load end-groups (head) and monomere up and down
ENDPATCH_FILE = "../DATA/endpatch_CH3.dat"
MONOMER_FILE = "../DATA/monomer.dat"

CONVERSION_TABLE = {1: "CC32A", 2: "OC30A", 3: "HCA2", 4: "OC311", 5: "HCP1",}

endpatch, monomer = load_molecule_templates(ENDPATCH_FILE, MONOMER_FILE)

check_input("monomer.dat", monomer, expected_atoms=7)
check_input("endpatch_CH3.dat", endpatch, expected_atoms=8)

# Place the atoms and build the peg
atoms = build_peg(endpatch, monomer, nseg, spacing_nm)
atoms = convert_coordinates(atoms)

# Validate atom count
expected_atoms = (2 * len(endpatch) + nseg * len(monomer))
assert atoms.shape == (expected_atoms, 7)

# Center the molecule 
atoms = center_molecule(atoms)          
assert np.allclose(atoms[:, 4:7].mean(axis=0), 0.0, atol=1e-10)
x, y, z = atoms[:,4], atoms[:,5], atoms[:,6]

### initial box size
Lx, Ly, Lz = np.max(x)*2+10, np.max(y)*2+10, np.max(z)*2+10

check_box_boundaries(atoms, Lx, Ly, Lz)
validate_atoms(atoms, nseg, monomer, endpatch)

car, hyd, oxy = analyze_atom_types(atoms, CONVERSION_TABLE)

molmass = estimate_molar_mass(atoms)

write_gro("peg.gro", atoms, CONVERSION_TABLE, box=(Lx, Ly, Lz))


# Build bonds
cc_bonds, cc_bond_types = build_cc_bonds(car, cutoff=spacing_A)
assert len(cc_bonds) == nseg+1
if "OH" in ENDPATCH_FILE:
    print("Warning : will not assign the right atom type of C-O-H terminated PEG")
co_bonds, co_bond_types = build_co_bonds(car, oxy, cutoff = 2)
if "CH3" in ENDPATCH_FILE:
    assert len(co_bonds) == len(oxy)*2
ch_bonds, ch_bond_types = build_ch_bonds(car, hyd, cutoff = 1.5)
if "CH3" in ENDPATCH_FILE:
    assert len(ch_bonds) == len(hyd)
oh_bonds, oh_bond_types = build_oh_bonds(oxy, hyd)

bonds = combine_bonds((cc_bonds, cc_bond_types),
                      (co_bonds, co_bond_types),
                      (ch_bonds, ch_bond_types),
                      (oh_bonds, oh_bond_types))

angles, angle_types = calculate_angles(atoms, bonds, CONVERSION_TABLE)
assert np.sum(angle_types == 8) == nseg + 2 
assert np.sum(angle_types == 3) + np.sum(angle_types == 7) == nseg*2 + 2

# Build dihedrals
dihedrals, dihedral_types = calculate_dihedrals(atoms, angles, CONVERSION_TABLE)
assert (np.sum(dihedral_types == 4) + np.sum(dihedral_types == 5))/2 == nseg + 1

check_dihedral_counts(atoms, dihedrals, CONVERSION_TABLE)
validate_dihedrals(atoms, dihedrals, dihedral_types)

# check angles and dihedrals
angles = combine_angles(angles, angle_types)
dihedrals = combine_dihedrals(dihedrals, dihedral_types)
check_dihedral_atom_ids(atoms, dihedrals)

masses = {1: 12.011, 2: 15.999, 3: 1.008, 4: 15.999, 5: 1.008}
write_lammps_data("single-peo.data", atoms, bonds, angles, dihedrals, Lx, Ly, Lz, masses)
