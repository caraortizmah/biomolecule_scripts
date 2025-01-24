import sys
import re
import numpy as np
from ase.io import read
from rdkit import Chem
#from rdkit.Chem import rdDistGeom, AllChem
from rdkit.Geometry import Point3D


# without ASE
def get_inertia_tensor(atoms, fullinfo=False):
    """
    Function that returns the inertia tensor
    of a molecule and its center of mass.
    Using this function is possible to
    get the moments of inertia and their axes and 
    to create distances between group of atoms along
    the main axe.
    """

    coords = atoms.get_positions()  # x, y, z coordinates
    
    # Take masses from the atoms
    masses = atoms.get_masses()
    
    # Calculate the center of mass, weighted by atomic masses
    center_of_mass = np.sum(masses[:, np.newaxis] *\
                            coords, axis=0) / np.sum(masses)

    # Translate coordinates to the center of mass
    translated_coords = coords - center_of_mass
    

    
    # Extract x, y, z coordinates into separate arrays
    x = translated_coords[:, 0]
    y = translated_coords[:, 1]
    z = translated_coords[:, 2]

    # Compute components of the inertia tensor using broadcasting
    Ixx = np.sum(masses * (y**2 + z**2))
    Iyy = np.sum(masses * (x**2 + z**2))
    Izz = np.sum(masses * (x**2 + y**2))
    Ixy = -np.sum(masses * (x * y))
    Ixz = -np.sum(masses * (x * z))
    Iyz = -np.sum(masses * (y * z))

    # Construct the inertia tensor matrix
    inertia_tensor = np.array([
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz]
    ])
    
    if not fullinfo:
        return inertia_tensor
    else:
        return inertia_tensor, center_of_mass

def create_AB_distance(atoms, delta_d, position_lim, armin=True):
    """
    Function that uses the moments of inertia
    of a molecule to create a separation between two group
    of atoms distributed in the eigenvector corresponding to 
    the smallest eigenvalue: the axis along changes of 
    inertia are lower; the main axis of separation.
    The distance separation is delta_d.
    position_lim is a number that indicate up to which atom 
    position will be considered as a group of atoms to be 
    translated along the main axis.
    Having armin True as default the smallest eigenvalue is 
    used from the principal moments of inertia.
    if armin is False, then the largest eigenvalue is used.
    """

    coords = atoms.get_positions()  # x, y, z coordinates

    # Get the intertia tensor and the center of the mass
    inertia_tensor = get_inertia_tensor(atoms)#, True)

    # Compute the principal moments of inertia (eigenvalues of the inertia tensor)
    principal_moments_of_inertia, axes = np.linalg.eig(inertia_tensor)
    
    if armin:
        # Main axis (smallest eigenvalue)
        main_axis = axes[:, np.argmin(principal_moments_of_inertia)]
    else:
        # Main axis (largest eigenvalue)
        main_axis = axes[:, np.argmax(principal_moments_of_inertia)]

    # Translate group A (using a position in a list of atoms as criteria)
    # along the main axis by delta_d
    translation_vector = delta_d * main_axis
    coords[0:position_lim] += translation_vector

    return coords

def create_AB_dist_guessinggroup(atoms, delta_d):
    """
    Function that uses the moments of inertia
    of a molecule to create a separation between two group
    of atoms distributed in the eigenvector corresponding to 
    the smallest eigenvalue (the main axis of separation).
    The distance separation is delta_d.
    The criteria of separation between group A and B will
    depend on the eigenvector that corresponds to the 
    smallest eigenvalue (axis along changes of inertia are lower).
    """

    coords = atoms.get_positions()  # x, y, z coordinates

    # Get the intertia tensor and the center of the mass
    inertia_tensor, com = get_inertia_tensor(atoms, True)

    # Compute the principal moments of inertia (eigenvalues of the inertia tensor)
    principal_moments_of_inertia, axes = np.linalg.eig(inertia_tensor)
    
    # Main axis (largest eigenvalue)
    main_axis = axes[:, np.argmin(principal_moments_of_inertia)]
    # Translate coordinates to the center of mass
    centered_coords = coords - com
    
    # Project coordinates onto the main axis
    projections = np.dot(centered_coords, main_axis)

    # Separate into two groups
    median_proj = np.median(projections)
    group_A_indices = np.where(projections < median_proj)[0]

    # Translate group A along the main axis by delta_d
    translation_vector = delta_d * main_axis
    centered_coords[group_A_indices] += translation_vector

    # Reconstruct the final coordinates
    coords_final = centered_coords + com

    return coords_final

def save_pdb_separation(pdb_file, delta_d, position_lim):
    """
    Call create_AB_distance() to create a new
    coordinates from the same structure separating
    a group of atoms, from the first on the pdb file
    until the position number atom_threshold_num to
    the rest of atoms by delta_d Angstroms.

    Args:
        pdb_file (str): PDB file.
        delta_d (float): Desired change in separation 
        distance (positive to increase,
         negative to decrease).
        atom_threshold_num (int): the atom position on
         the pdb file where represents one group of
         atoms, i.e. amino acid.
    
    Returns:
        output_pdb (str): The adjusted PDB file.
    """

    structure_atoms = read(pdb_file)

    # Create a RDKit molecule object for the new coordinates
    #  having a new distance separation
    mol1_1 = Chem.RWMol()
    for atom in structure_atoms:
        atom_num = int(atom.number)  # Atomic number from ASE
        mol1_1.AddAtom(Chem.Atom(atom_num))

    # Create another RDKit molecule having another coordinates
    mol1_2 = Chem.RWMol()
    for atom in structure_atoms:
        atom_num = int(atom.number)  # Atomic number from ASE
        mol1_2.AddAtom(Chem.Atom(atom_num))

    # Add 3D coordinates to RDKit molecule
    conf1_1 = Chem.Conformer(mol1_1.GetNumAtoms())
    for i, atom_coords in enumerate(\
        create_AB_distance(structure_atoms, \
                           delta_d, position_lim)):
        x, y, z = atom_coords  # Atom coordinates from the function
        conf1_1.SetAtomPosition(i, Point3D(x, y, z))
    mol1_1.AddConformer(conf1_1)

    conf1_2 = Chem.Conformer(mol1_2.GetNumAtoms())
    for i, atom_coords in enumerate(\
        create_AB_distance(structure_atoms, \
                           delta_d, position_lim, \
                            armin=False)):
        x, y, z = atom_coords  # Atom coordinates from the function
        conf1_2.SetAtomPosition(i, Point3D(x, y, z))
    mol1_2.AddConformer(conf1_2)

    # Save to a PDB file
    suffix="v1"
    new_pdb = re.sub(r"(\.pdb)$", f".{suffix}\\1", pdb_file)
    Chem.rdmolfiles.MolToPDBFile(mol1_1, new_pdb)

    # Save to a PDB file
    suffix="v2"
    new_pdb = re.sub(r"(\.pdb)$", f".{suffix}\\1", pdb_file)
    Chem.rdmolfiles.MolToPDBFile(mol1_2, new_pdb)

def main():

    # Check for the correct number of arguments
    if len(sys.argv) != 4:
        print("Usage: python create_distance_pairaa.py pdb_list.txt distance_separation nth_atom")
        sys.exit(1)

    # Parse command-line arguments
    filename = sys.argv[1]
    delta_d = float(sys.argv[2])
    position_lim = int(sys.argv[3])
    
    with open(filename, 'r') as file:
        for line in file:
            current_file = line.strip()
            save_pdb_separation(current_file, delta_d, position_lim)

if __name__ == "__main__":
    main()
