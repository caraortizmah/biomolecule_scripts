import sys
import numpy as np
from ase.io import read
import pandas as pd


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

def get_pi_stacking_features(pdb_file, atm_indx_planes, atm_indx_rings): 
    """
    Function that uses a pdb file and position atoms
    in that file as indices to extract structural 
    features related to (specifically) pi-stacking
    between atomatic amino acids.

    Args:
        pdb_file (str): PDB file.
        atm_indx_planes (numpy.ndarray): A 6-element array that 
         contains position numbers of the target atoms in the
         PDB file. The 3 first elements should correspond to 
         three atoms of the aromatic side chain of one amino acid
         and the last 3 elements of the array should correspond to
         other three atoms of another aromatic side chain of other
         amino acid.
        atm_indx_rings (numpy.ndarray): A 6-element array that 
         contains position numbers of the target atoms in the
         PDB file. The 3 first elements should correspond to 
         three atoms of the bencene of one amino acid and the 
         last 3 elements of the array should correspond to
         other three atoms of another bencene of other amino acid.
    
    Returns:
        output_pdb (numpy.ndarray): Pi-stacking features ordered as
         1. Pi stacking distance, 2. Dihedral angle between both 
         planes of the aromatic rings, 3. In-plane displacement 
         from the bencene ring of the aromatic side chains,
         4. Principal moments of inertia of the whole molecule and
         its 5. Center of mass.
         
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
    
    return [pi_stacking_distance, dihedral_angle, in_plane_displacement,\
           principal_moments_of_inertia, com]
    #return np.array([pi_stacking_distance, dihedral_angle, in_plane_displacement,\
    #       principal_moments_of_inertia, com], dtype=object)

def save_features(pdb_file, indx_planes, indx_rings, filename="pi_features_results"):
    """
    Function that saves the computed features by the 
     function get_pi_stacking_features() in a csv format.
    """
    
    # Set header for output files
    header = ["Pdb_file"] + ["Pi_stacking_distance"] + ["Dihedral_angle"] + \
             ["in_place_displacement"] + ["principal_moments_of_intertia"] + \
             ["center_of_mass"]

    # Hand-craft extraction of structural features based on pi stacking 
    pi_features = get_pi_stacking_features(pdb_file, indx_planes, indx_rings)
    
    data = [pdb_file] + pi_features
    
    # Save as CSV
    df = pd.DataFrame([data], columns=header)
    # Append to the file if it exists, create a new file if it doesn't
    df.to_csv(f"{filename}.csv", index=False, mode='a', \
              header=not pd.io.common.file_exists(f"{filename}.csv"))
    # Last command ensures that the header is written only once
    
def main():

    import ast

    # Check for the correct number of arguments
    if len(sys.argv) != 5:
        print("Usage: python pi_feat_extraction.py pdb_list.txt \"[1, 2, 3, 4, 5, 6]\"")
        print("  \"[1, 2, 3, 4, 5, 6]\" name_file_pi_features_results")
        print("Where the two arrays are the same position atoms in all the pdb files in the")
        print("  pdb_list.txt")
        sys.exit(1)

    # Parse command-line arguments
    filename = sys.argv[1]
    plane_array = ast.literal_eval(sys.argv[2]) # Converts "[1,2,3,4,5,6]" to [1, 2, 3, 4, 5, 6]
    plane_array = np.array(plane_array) - 1
    ring_array = ast.literal_eval(sys.argv[3])
    ring_array = np.array(ring_array) - 1
    output_file = sys.argv[4]
    
    with open(filename, 'r') as file:
        for line in file:
            current_file = line.strip()
            save_features(current_file,plane_array,ring_array,output_file)

if __name__ == "__main__":
    main()
