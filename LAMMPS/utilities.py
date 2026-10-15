

import numpy as np


def check_input(name, data, expected_atoms=None):
    assert data.ndim == 2 and data.shape[1] == 6, (
        f"{name}: expected 6 columns, got {data.shape}"
    )
    assert np.isfinite(data).all(), f"{name}: NaN or infinite values"
    assert np.all(data[:, 0] == np.arange(1, len(data) + 1)), (
        f"{name}: atom IDs must be consecutive, starting at 1"
    )
    assert np.all(data[:, 1] == data[:, 1].astype(int)), (
        f"{name}: atom types must be integers"
    )
    assert np.all(np.isin(data[:, 1].astype(int), [1, 2, 3, 4, 5])), (
        f"{name}: unknown atom type"
    )
    if expected_atoms is not None:
        assert len(data) == expected_atoms, (
            f"{name}: expected {expected_atoms} atoms, got {len(data)}"
        )
    print(f"[OK] {name}: {len(data)} atoms, valid input")


def convert_coordinates(atoms, factor=10.0):
    """Convert the x, y and z columns from nm to Angstroms."""

    atoms = atoms.copy()
    atoms[:, 3:6] *= factor

    return atoms


def load_molecule_templates(endpatch_file, monomer_file):
    """Load and validate the end-patch and monomer templates."""

    endpatch = np.loadtxt(endpatch_file, ndmin=2)
    monomer = np.loadtxt(monomer_file, ndmin=2)

    check_input(
        "monomer.dat", monomer, expected_atoms=7
    )
    check_input(
        "endpatch_CH3.dat", endpatch, expected_atoms=8
    )

    # Each row must contain:
    # local ID, atom type, x, y, z, additional property
    assert endpatch.shape[1] == 6
    assert monomer.shape[1] == 6

    return endpatch, monomer

def build_peg(endpatch, monomer, nseg, spacing):
    """Build one PEG molecule from two end patches and nseg monomers."""

    if nseg < 1:
        raise ValueError("nseg must be at least 1.")

    if spacing <= 0:
        raise ValueError("spacing must be positive.")

    n_atoms = 2 * len(endpatch) + nseg * len(monomer)
    atoms = np.zeros((n_atoms, 7), dtype=float)
    cursor = 0

    # Columns: atom ID, molecule ID, atom type, charge, x, y, z

    # First end patch
    for m in endpatch:
        atoms[cursor] = (
            cursor + 1, 1, m[1], m[2], -m[3]-spacing, -m[4], -m[5])
        cursor += 1

    # Repeated monomers
    for seg in range(nseg):
        for m in monomer:
            atoms[cursor] = (cursor + 1, 1, m[1], m[2], m[3] + seg * spacing * 2, m[4], m[5])
            cursor += 1

    # Second end patch
    end_offset = (2 * (nseg - 1)) * spacing + spacing

    for m in endpatch:
        atoms[cursor] = (cursor + 1, 1, m[1], m[2], m[3]+end_offset, m[4], m[5])
        cursor += 1

    assert cursor == n_atoms
    return atoms

def analyze_atom_types(atoms, conversion_table):
    """Print atom counts and return useful atom selections."""

    atom_types = atoms[:, 2].astype(int)

    print(f"\nTotal atoms: {len(atoms)}")

    for atom_type, label in conversion_table.items():
        count = np.count_nonzero(atom_types == atom_type)
        print(f"{label}: {count} atoms")

    # Atom selections
    car = atoms[atom_types == 1]

    hyd = atoms[
        (atom_types == 3) | (atom_types == 5)
    ]

    oxy = atoms[
        (atom_types == 2) | (atom_types == 4)
    ]

    return car, hyd, oxy

def center_molecule(atoms):
    """Center the molecule at the origin using its geometric centroid."""

    atoms = atoms.copy()

    # Center x, y and z coordinates only
    coordinates = atoms[:, 3:6]
    center = coordinates.mean(axis=0)

    atoms[:, 3:6] -= center

    return atoms

def check_box_boundaries(atoms, Lx, Ly, Lz):
    """Check whether all atoms lie inside the simulation box."""

    box_lengths = [Lx, Ly, Lz]
    axis_names = ["x", "y", "z"]

    coordinates = atoms[:, 3:6]

    for col, length, axis in zip(range(3), box_lengths, axis_names):
        lo = coordinates[:, col].min()
        hi = coordinates[:, col].max()

        print(f"[CHECK] {axis}: [{lo:.3f}, {hi:.3f}] Å")

        if lo < -length / 2 or hi > length / 2:
            raise ValueError(
                f"Atoms outside the {axis}-box: "
                f"[{lo:.3f}, {hi:.3f}] Å, "
                f"box limits are [{-length/2:.3f}, "
                f"{length/2:.3f}] Å"
            )

    print("\n[OK] All atoms are inside the simulation box")

def validate_atoms(atoms, Nseg, monomer, endpatch):
    """Validate atom count, IDs, molecule IDs, types and coordinates."""

    expected_atoms = Nseg * len(monomer) + 2 * len(endpatch)

    assert atoms.ndim == 2 and atoms.shape[1] == 7, f"Expected an (N, 7) array, got {atoms.shape}"
    assert len(atoms) == expected_atoms, f"Expected {expected_atoms} atoms, got {len(atoms)}"
    assert np.isfinite

def estimate_molar_mass(atoms):
    """Estimate PEG molar mass in g/mol from atom types."""

    atom_types = atoms[:, 2].astype(int)
    car = np.count_nonzero(atom_types == 1)
    oxy = np.count_nonzero(np.isin(atom_types, [2, 4]))
    hyd = np.count_nonzero(np.isin(atom_types, [3, 5]))
    molmass = car * 12 + oxy * 16 + hyd * 1

    print(f"\nPEG molecule with {molmass} g/mol")
    return molmass

def build_cc_bonds(car):
    """Build carbon-carbon bonds between consecutive carbon pairs."""
    bonds = []
    for idx in range(0, len(car) - 1, 2):
        idx0 = int(car[idx, 0])
        idx1 = int(car[idx + 1, 0])
        bonds.append(sorted([idx0, idx1]))
    bonds = np.array(bonds, dtype=int).reshape(-1, 2)
    bond_types = np.ones(len(bonds), dtype=int)
    print(f"\n[OK] Generated {len(bonds)} carbon-carbon bonds")
    return bonds, bond_types

def build_co_bonds(car, oxy, cutoff=4):
    """Build carbon-oxygen bonds based on the distance cutoff in Å."""
    bonds = []
    for c in car:
        for o in oxy:
            distance = np.linalg.norm(c[4:7] - o[4:7])
            if distance <= cutoff:
                bonds.append(sorted([int(c[0]), int(o[0])]))
    bonds = np.array(bonds, dtype=int).reshape(-1, 2)
    bond_types = np.full(len(bonds), 2, dtype=int)
    print(f"[OK] Generated {len(bonds)} carbon-oxygen bonds")
    return bonds, bond_types

def build_ch_bonds(car, hyd, cutoff=2.5):
    """Build carbon-hydrogen bonds based on the distance cutoff in Å."""
    bonds = []
    for c in car:
        for h in hyd:
            distance = np.linalg.norm(c[4:7] - h[4:7])
            if distance <= cutoff:
                bonds.append(sorted([int(c[0]), int(h[0])]))
    bonds = np.array(bonds, dtype=int).reshape(-1, 2)
    bond_types = np.full(len(bonds), 3, dtype=int)
    print(f"[OK] Generated {len(bonds)} carbon-hydrogen bonds")
    return bonds, bond_types

def build_oh_bonds(oxy, hyd, cutoff=1.5):
    """Build oxygen-hydrogen bonds based on the distance cutoff in Å."""
    bonds = []
    for o in oxy:
        for h in hyd:
            distance = np.linalg.norm(o[4:7] - h[4:7])
            if distance <= cutoff:
                bonds.append(sorted([int(o[0]), int(h[0])]))
    bonds = np.array(bonds, dtype=int).reshape(-1, 2)
    bond_types = np.full(len(bonds), 4, dtype=int)
    print(f"[OK] Generated {len(bonds)} oxygen-hydrogen bonds")
    return bonds, bond_types

def combine_bonds(*bond_groups):
    """Combine bond arrays and their corresponding bond types."""
    bonds = np.vstack([group[0] for group in bond_groups if len(group[0]) > 0]) if any(len(group[0]) > 0 for group in bond_groups) else np.empty((0, 2), dtype=int)
    bond_types = np.concatenate([group[1] for group in bond_groups if len(group[1]) > 0]) if any(len(group[1]) > 0 for group in bond_groups) else np.empty(0, dtype=int)
    bonds = np.column_stack((np.arange(1, len(bonds) + 1), bond_types, bonds))
    print(f"[OK] Combined {len(bonds)} bonds")
    return bonds

def check_bond_counts(atoms, bonds):
    """Compare generated bond counts with estimates for linear PEG with CH3 ending."""
    atom_types = atoms[:, 2].astype(int)
    n_car = np.count_nonzero(atom_types == 1)
    n_oxy = np.count_nonzero(np.isin(atom_types, [2, 4]))
    n_hyd = np.count_nonzero(np.isin(atom_types, [3, 5]))
    actual = {bond_type: np.count_nonzero(bonds[:, 1] == bond_type) for bond_type in [1, 2, 3, 4]}
    expected_cc = n_car // 2
    expected_co = n_oxy * 2
    expected_ch = n_hyd
    print(f"\nC-C bonds: expected approximately {expected_cc}, generated {actual[1]}")
    print(f"C-O bonds: expected approximately {expected_co}, generated {actual[2]}")
    print(f"C-H bonds: expected approximately {expected_ch}, generated {actual[3]}")

def check_missing_ch_bonds(car, hyd, ch_bonds):
    """Find hydrogen atoms that have no generated carbon-hydrogen bond."""
    bonded_hydrogens = set(ch_bonds[:, 1]) | set(ch_bonds[:, 0]) if len(ch_bonds) else set()
    missing = hyd[~np.isin(hyd[:, 0].astype(int), list(bonded_hydrogens))]
    if len(missing) > 0:
        print(f"\nHydrogens without C-H bonds: {len(missing)}")
        for h in missing:
            distances = np.linalg.norm(car[:, 4:7] - h[4:7], axis=1)
            nearest = np.argmin(distances)
            print(f"H atom {int(h[0])}, type {int(h[2])}: nearest C atom {int(car[nearest, 0])}, distance {distances[nearest]:.3f} Å")

def calculate_angles(atoms, bonds, conversion_table):
    """Calculate angles and assign force-field angle types."""
    angles = []
    angle_types = []
    atom_types = {int(atom[0]): int(atom[2]) for atom in atoms}
    bonded_a = np.concatenate((bonds[:, 2], bonds[:, 3]))

    for a in atoms:
        ida = int(a[0])
        tpa = int(a[2])
        if np.count_nonzero(bonded_a == ida) > 1:
            id_neighbors = np.unique(bonds[(bonds[:, 2] == ida) | (bonds[:, 3] == ida), 2:4])
            id_neighbors = id_neighbors[id_neighbors != ida]
            for idx, idb in enumerate(id_neighbors):
                for idc in id_neighbors[idx + 1:]:
                    tpb = atom_types[int(idb)]
                    tpc = atom_types[int(idc)]
                    angle_type = None

                    if (tpb == 3) and (tpa == 1) and (tpc == 3):
                        angle_type = 5
                    elif ((tpb == 3) and (tpa == 1) and (tpc == 4)) or ((tpb == 4) and (tpa == 1) and (tpc == 3)):
                        angle_type = 6
                    elif ((tpb == 3) and (tpa == 1) and (tpc == 1)) or ((tpb == 1) and (tpa == 1) and (tpc == 3)):
                        angle_type = 2
                    elif ((tpb == 4) and (tpa == 1) and (tpc == 1)) or ((tpb == 1) and (tpa == 1) and (tpc == 4)):
                        angle_type = 7
                    elif ((tpb == 1) and (tpa == 4) and (tpc == 5)) or ((tpb == 5) and (tpa == 4) and (tpc == 1)):
                        angle_type = 1
                    elif ((tpb == 1) and (tpa == 1) and (tpc == 2)) or ((tpb == 2) and (tpa == 1) and (tpc == 1)):
                        angle_type = 3
                    elif ((tpb == 3) and (tpa == 1) and (tpc == 2)) or ((tpb == 2) and (tpa == 1) and (tpc == 3)):
                        angle_type = 4
                    elif (tpb == 1) and (tpa == 2) and (tpc == 1):
                        angle_type = 8
                    else:
                        print(tpb, tpa, tpc)
                        print(conversion_table[tpb], conversion_table[tpa], conversion_table[tpc])
                        raise ValueError("Unknown angle type")

                    angles.append([int(idb), ida, int(idc)])
                    angle_types.append(angle_type)

    angles = np.array(angles, dtype=int).reshape(-1, 3)
    angle_types = np.array(angle_types, dtype=int)
    print(f"[OK] Generated {len(angles)} angles")
    return angles, angle_types

def check_angle_counts(atoms, angles, CONVERSION_TABLE):
    """Check angle counts and report the total number of angles."""
    atom_types = {int(atom[0]): int(atom[2]) for atom in atoms}
    angle_counts = {}
    for angle in angles:
        types = tuple(atom_types[int(atom_id)] for atom_id in angle)
        angle_counts[types] = angle_counts.get(types, 0) + 1
    print(f"[OK] Total angles: {len(angles)}")
    for types, count in sorted(angle_counts.items()):
        print(f"  {count:4d}  {' - '.join(CONVERSION_TABLE[t] for t in types)}")

def calculate_dihedrals(atoms, angles, CONVERSION_TABLE):
    """Generate dihedrals from pairs of angles sharing a central bond."""

    dihedrals = []
    dihedral_types = []

    atom_types = {
        int(atom[0]): int(atom[2])
        for atom in atoms
    }

    # Force-field dihedral type IDs.
    rules = {
        (1, 1, 2, 1): [1, 2, 3],
        (2, 1, 1, 2): [4, 5],
        (1, 1, 4, 5): [6, 7, 8],
        (4, 1, 1, 4): [9, 10, 11],
        (4, 1, 1, 2): [12, 13, 14],
    }

    print("Raw angles (atom IDs):")
    for angle in angles:
        ids = tuple(map(int, angle))
        names = tuple(
            CONVERSION_TABLE[int(atoms[atom_id - 1, 2])]
            for atom_id in ids
        )
        print(ids, names)

        seen = set()

        for angle1 in angles:
            i, j, k = map(int, angle1)

            for angle2 in angles:
                a, b, c = map(int, angle2)

                # Join angles across the bond j-k.
                # The shared bond can appear in either orientation.
                if a == j and b == k:
                    path = (i, j, k, c)
                elif b == j and a == k:
                    path = (i, j, k, c)
                elif b == k and c == j:
                    path = (i, j, k, a)
                elif c == k and b == j:
                    path = (i, j, k, a)
                else:
                    continue

                if len(set(path)) != 4:
                    continue

                canonical = min(path, path[::-1])

                if canonical in seen:
                    continue

                seen.add(canonical)

                types = tuple(atom_types[x] for x in canonical)
                type_names = tuple(CONVERSION_TABLE[t] for t in types)

                print("Candidate:", canonical, type_names)

                matched = rules.get(types)

                if matched is None:
                    matched = rules.get(types[::-1])

                if matched is None:
                    continue

                for dtype in matched:
                    dihedrals.append(canonical)
                    dihedral_types.append(dtype)

    dihedrals = np.array(dihedrals, dtype=int).reshape(-1, 4)
    dihedral_types = np.array(dihedral_types, dtype=int)

    print(f"[OK] Generated {len(dihedrals)} dihedrals")

    return dihedrals, dihedral_types

def check_dihedral_counts(atoms, dihedrals, CONVERSION_TABLE):
    """Check dihedral counts and report the total number of dihedrals."""
    atom_types = {int(atom[0]): int(atom[2]) for atom in atoms}
    dihedral_counts = {}
    for dihedral in dihedrals:
        types = tuple(atom_types[int(atom_id)] for atom_id in dihedral)
        dihedral_counts[types] = dihedral_counts.get(types, 0) + 1
    print(f"[OK] Total dihedrals: {len(dihedrals)}")
    for types, count in sorted(dihedral_counts.items()):
        print(f"  {count:4d}  {' - '.join(CONVERSION_TABLE[t] for t in types)}")

def validate_dihedrals(atoms, dihedrals, dihedral_types):
    """Validate dihedral atom IDs and type assignments."""
    atom_ids = set(atoms[:, 0].astype(int))
    assert dihedrals.ndim == 2 and dihedrals.shape[1] == 4
    assert len(dihedrals) == len(dihedral_types)
    assert np.all(np.isin(dihedrals, list(atom_ids)))
    assert np.all(dihedrals[:, 0] != dihedrals[:, 1])
    assert np.all(dihedrals[:, 1] != dihedrals[:, 2])
    assert np.all(dihedrals[:, 2] != dihedrals[:, 3])
    assert np.all(np.isfinite(dihedral_types))
    print("[OK] Dihedrals validated")

def combine_angles(angles, angle_types):
    """Combine angle IDs, types, and atom IDs into a LAMMPS array."""
    assert len(angles) == len(angle_types)
    angle_ids = np.arange(1, len(angles) + 1).reshape(-1, 1)
    combined_angles = np.column_stack((angle_ids, angle_types, angles))
    print(f"[OK] Combined {len(combined_angles)} angles")
    return combined_angles.astype(int)

def combine_dihedrals(dihedrals, dihedral_types):
    """Combine dihedral IDs, types, and atom IDs into a LAMMPS array."""
    assert len(dihedrals) == len(dihedral_types)
    dihedral_ids = np.arange(1, len(dihedrals) + 1).reshape(-1, 1)
    combined_dihedrals = np.column_stack((dihedral_ids, dihedral_types, dihedrals))
    print(f"[OK] Combined {len(combined_dihedrals)} dihedrals")
    return combined_dihedrals.astype(int)

def write_lammps_data(filename, atoms, bonds, angles, dihedrals,
                      Lx, Ly, Lz, masses):
    """Write molecular topology and coordinates to a LAMMPS data file."""

    with open(filename, "w") as f:
        f.write("PEG molecule\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write(f"{len(dihedrals)} dihedrals\n\n")

        f.write("5 atom types\n")
        f.write("5 bond types\n")
        f.write("8 angle types\n")
        f.write("14 dihedral types\n\n")

        f.write(f"{-Lx/2:.6f} {Lx/2:.6f} xlo xhi\n")
        f.write(f"{-Ly/2:.6f} {Ly/2:.6f} ylo yhi\n")
        f.write(f"{-Lz/2:.6f} {Lz/2:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for atom_type, mass in sorted(masses.items()):
            f.write(f"{atom_type} {mass:.6f}\n")

        f.write("\nAtoms # full\n\n")
        np.savetxt(
            f,
            atoms[:, [0, 1, 2, 3, 4, 5, 6]],
            fmt="%d %d %d %.6f %.6f %.6f %.6f"
        )

        f.write("\nBonds\n\n")
        np.savetxt(f, bonds, fmt="%d")

        f.write("\nAngles\n\n")
        np.savetxt(f, angles, fmt="%d")

        if len(dihedrals) > 0:
            f.write("\nDihedrals\n\n")
            np.savetxt(f, dihedrals, fmt="%d")

    print(f"[OK] LAMMPS data written to {filename}")

def check_dihedral_atom_ids(atoms, dihedrals):
    """Check that dihedral atom IDs exist and each quartet is valid."""
    atom_ids = set(atoms[:, 0].astype(int))
    for d in dihedrals:
        if not all(int(i) in atom_ids for i in d[2:6]):
            print(f"[ERROR] Invalid atom ID in dihedral: {d}")
        if len(set(d[2:6])) != 4:
            print(f"[ERROR] Repeated atom IDs in dihedral: {d}")
    print("[OK] Dihedral atom ID check complete")

def write_gro(filename, atoms, conversion_table, box, title="PEG molecule"):
    """Write atom coordinates to a GROMACS .gro file."""

    with open(filename, "w") as f:
        f.write(f"{title}\n")
        f.write(f"{len(atoms):5d}\n")

        for atom in atoms:
            atom_id = int(atom[0])
            atom_type = int(atom[2])
            x, y, z = atom[4:7] / 10.0
            atom_name = conversion_table[atom_type]

            f.write(
                f"{1:5d}"
                f"{'PEG':<5s}"
                f"{atom_name:>5s}"
                f"{atom_id % 100000:5d}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
            )

        bx, by, bz = box
        f.write(f"{bx / 10:10.5f}{by / 10:10.5f}{bz / 10:10.5f}\n")
