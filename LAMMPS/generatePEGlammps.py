import numpy as np
import warnings
warnings.filterwarnings('ignore')

from utilities import check_input, convert_coordinates, load_molecule_templates, \
    build_peg, analyze_atom_types, center_molecule, check_box_boundaries, \
    validate_atoms, estimate_molar_mass, build_cc_bonds, build_co_bonds, \
    build_ch_bonds, build_oh_bonds, combine_bonds, check_bond_counts, \
    check_missing_ch_bonds, calculate_angles, check_angle_counts, calculate_dihedrals, \
    check_dihedral_counts, validate_dihedrals, combine_angles, combine_dihedrals, \
    write_lammps_data, check_dihedral_atom_ids, write_gro

### number of segments for the PEG molecules

nseg = 1
spacing = 0.28 # distance between 2 monomers

### initial box size
Lx, Ly, Lz = 200, 200, 200

### load end-groups (head) and monomere up and down
ENDPATCH_FILE = "../DATA/endpatch_CH3.dat"
MONOMER_FILE = "../DATA/monomer.dat"

CONVERSION_TABLE = {1: "CC32A", 2: "OC30A", 3: "HCA2", 4: "OC311", 5: "HCP1",}

endpatch, monomer = load_molecule_templates(ENDPATCH_FILE, MONOMER_FILE)

check_input("monomer.dat", monomer, expected_atoms=7)
check_input("endpatch_CH3.dat", endpatch, expected_atoms=8)

atoms = build_peg(endpatch, monomer, nseg, spacing)

# Convert coordinates
atoms = convert_coordinates(atoms)

# Validate atom count
expected_atoms = (2 * len(endpatch) + nseg * len(monomer))
assert atoms.shape == (expected_atoms, 7)

# Center the molecule 
atoms = center_molecule(atoms)

assert np.allclose(atoms[:, 4:7].mean(axis=0), 0.0, atol=1e-10)

print("\n[OK] Molecule centered at the origin")

check_box_boundaries(atoms, Lx, Ly, Lz)

validate_atoms(atoms, nseg, monomer, endpatch)

car, hyd, oxy = analyze_atom_types(atoms, CONVERSION_TABLE)

molmass = estimate_molar_mass(atoms)

write_gro("peg.gro", atoms, CONVERSION_TABLE, box=(Lx, Ly, Lz))

# Build bonds
cc_bonds, cc_bond_types = build_cc_bonds(car)
co_bonds, co_bond_types = build_co_bonds(car, oxy)
ch_bonds, ch_bond_types = build_ch_bonds(car, hyd)
oh_bonds, oh_bond_types = build_oh_bonds(oxy, hyd)

bonds = combine_bonds((cc_bonds, cc_bond_types),
                      (co_bonds, co_bond_types),
                      (ch_bonds, ch_bond_types),
                      (oh_bonds, oh_bond_types))

check_bond_counts(atoms, bonds)

check_missing_ch_bonds(car, hyd, ch_bonds)

# Build angles
angles, angle_types = calculate_angles(atoms, bonds, CONVERSION_TABLE)
check_angle_counts(atoms, angles, CONVERSION_TABLE)

# Build dihedrals
dihedrals, dihedral_types = calculate_dihedrals(atoms, angles, CONVERSION_TABLE)

check_dihedral_counts(atoms, dihedrals, CONVERSION_TABLE)
validate_dihedrals(atoms, dihedrals, dihedral_types)

# check angles and dihedrals
angles = combine_angles(angles, angle_types)
dihedrals = combine_dihedrals(dihedrals, dihedral_types)
check_dihedral_atom_ids(atoms, dihedrals)

masses = {1: 12.011, 2: 15.999, 3: 1.008, 4: 15.999, 5: 1.008}
write_lammps_data("single-peo.data", atoms, bonds, angles, dihedrals, Lx, Ly, Lz, masses)
