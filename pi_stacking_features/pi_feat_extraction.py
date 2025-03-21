import sys
import re
import numpy as np
from ase.io import read
from rdkit import Chem
#from rdkit.Chem import rdDistGeom, AllChem
from rdkit.Geometry import Point3D


def get_aromatic_dihedral(coords1, coords2):
    """
    Function that returns the dihedral angle
    between the two planes formed by the 
    two aromatic side chains of the two 
    atomatic rings.
    Each plane is built by using a set of
    three atom coordinates that belong to 
    the aromatic ring of each amino acid.
    """

    # Compute normal vectors of both planes
    normal1 = np.cross(coords1[1] - coords1[0], coords1[2] - coords1[0])
    normal2 = np.cross(coords2[1] - coords2[0], coords2[2] - coords2[0])

    # Normalize the normal vectors
    normal1 /= np.linalg.norm(normal1)
    normal2 /= np.linalg.norm(normal2)

    # Compute the angle between the normal vectors
    cos_theta = np.dot(normal1, normal2)
    angle = np.arccos(np.clip(cos_theta, -1.0, 1.0)) * 180 / np.pi

    return angle

def get_centroids(coords):
    """
    Calculate the centroid of the 3D coordinates
    returning an [x, y, z] mean array.
    """
    return np.mean(coords, axis=0)

def distance_centroids(coords1, coords2):
    """
    Function that calculate the distance between 
    two group of atoms.
    Each group of atoms contains three atoms, the minimum
    number points to form a plane.
    From each plane is calculated a centroid and the 
    euclidean distance between both of them is returned.
    """

    # Compute the centroid for both planes
    centroid1 = get_centroids(coords1)
    centroid2 = get_centroids(coords2)
    return np.linalg.norm(centroid2 - centroid1)

def translate_coord(coords, target_centroid):
    """
    Translate a group of coordinates to a given target
    centroid of other group of atoms.
    """
    
    # Calculate the centroid
    centroid = get_centroids(coords)
    # Translate coordinates
    return coords + (target_centroid - centroid)

def get_rmsd_shift_centroid(coords1, coords2):
    """
    Function that calculates the rmsd between two
    group of coordinates after doing a translation
    of the coordinates of the second group of atoms.
    """
    
    # Calculate centroid of group 1
    centroid1 = get_centroids(coords1)

    # Translate group 2 coordinates according to the
    # centroids difference between both groups
    translated_coords = translate_coord(coords2, centroid1)

    # Calculate the RMSD between two groups
    diff = coords1 - translated_coords
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd

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

#*
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
        # Translate group A (using a position in a list of atoms as criteria)
        # along the main axis by delta_d
        translation_vector = delta_d * main_axis
    else:
        # Main axis (largest eigenvalue)
        main_axis = axes[:, np.argmin(principal_moments_of_inertia)]
        main_axis2 = axes[:, np.argmax(principal_moments_of_inertia)]
        # Translate group A (using a position in a list of atoms as criteria)
        # along the main axis by delta_d
        translation_vector = delta_d * 2 \
            * ( main_axis - main_axis2 )


    coords[0:position_lim] += translation_vector

    return coords
#*
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

#get pi_stacking related features
def get_pi_stacking_features(pdb_file, atm_indx_planes, atm_indx_rings): 
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
        suffix (str): the suffix for the new pdb version
    
    Returns:
        output_pdb (str): The adjusted PDB file.
    """

    structure_atoms = read(pdb_file)

    # Extract only atom coordinates for the two planes
    coords_planes = structure_atoms.get_positions()[atm_indx_planes]

    # Extract only atom coordinates for the two bencenes
    coords_rings = structure_atoms.get_positions()[atm_indx_rings]

    """
    coords_planes are the 3D coordinates from the pdb file
    and atm_indx_planes is a 6-element array.
    The first three for one aromatic ring plane and the
    last three for the other plane.
    coords_rings are the 3D coordinates from the pdb file
    and atm_indx_rings is a 6-element array.
    The first three for one bencene ring plane and the
    last three for the other bencene ring.
    """
    
    # Get coordinates representing a 2D plane
    plane1, plane2 = coords_planes[0:3], coords_planes[-3:]

    # Get coordinates representing a reduced bencene ring
    ring1, ring2 = coords_rings[0:3], coords_rings[-3:]

    # Calculate the dihedral angle of two planes of the aromatic side chains
    dihedral_angle = get_aromatic_dihedral(plane1, plane2)

    # Calculate a distance between the two aromatic rings (approximation)
    pi_stacking_distance = distance_centroids(plane1, plane2)

    # Calculate displacement of the two aromatic side chains
    in_plane_displacement = get_rmsd_shift_centroid(ring1, ring2)

    # Get the intertia tensor and the center of the mass
    inertia_tensor, com = get_inertia_tensor(structure_atoms, True)

    # Compute the principal moments of inertia (eigenvalues of the inertia tensor)
    principal_moments_of_inertia, axes = np.linalg.eig(inertia_tensor)

    return pi_stacking_distance, dihedral_angle, in_plane_displacement,\
           principal_moments_of_inertia, com

def save_features(pdb_file, indx_planes, indx_rings):
    """
    """

    pi_features = get_pi_stacking_features(pdb_file, indx_planes, indx_rings)

    print (pdb_file, pi_features)

    # Create a RDKit molecule object for the new coordinates
    # creating an artificial plane with 3 atoms of the aromatic ring
    plane_1 = Chem.RWMol()
    # Extracting only 3 atoms using the indices
    for indx in atom_indices[0:3]:
        atom = structure_atoms[indx]
        atom_num = int(atom.number)
        plane_1.AddAtom(Chem.Atom(atom_num))

    plane_2 = Chem.RWMol()
    # Extracting only 3 atoms using of the indices
    for indx in atom_indices[-3:]:
        atom = structure_atoms[indx]
        atom_num = int(atom.number)
        plane_2.AddAtom(Chem.Atom(atom_num))

    #for atom in structure_atoms:
    #    atom_num = int(atom.number)  # Atomic number from ASE
    #    mol1_1.AddAtom(Chem.Atom(atom_num))

    # Create another RDKit molecule having another coordinates
    #mol1_2 = Chem.RWMol()
    #for atom in structure_atoms:
    #    atom_num = int(atom.number)  # Atomic number from ASE
    #    mol1_2.AddAtom(Chem.Atom(atom_num))

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
    new_pdb = re.sub(r"(\.pdb)$", f".v1.{suffix}\\1", pdb_file)
    Chem.rdmolfiles.MolToPDBFile(mol1_1, new_pdb)

    # Save to a PDB file
    new_pdb = re.sub(r"(\.pdb)$", f".v2.{suffix}\\1", pdb_file)
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
            save_pdb_separation(current_file, delta_d, position_lim, "i")
            save_pdb_separation(current_file, -delta_d, position_lim, "d")

if __name__ == "__main__":
    main()
