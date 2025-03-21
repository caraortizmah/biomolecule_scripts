import sys
from openbabel import openbabel

def xyz_to_pdb(xyz_file, pdb_file):
    # Initialize Open Babel conversion object
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("xyz")
    obConversion.SetOutFormat("pdb")
    
    # Create an Open Babel molecule object
    mol = openbabel.OBMol()
    
    # Read XYZ file
    obConversion.ReadFile(mol, xyz_file)
    
    # Write PDB file
    obConversion.WriteFile(mol, pdb_file)

def main():

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # convert xyz to pdb
    xyz_to_pdb(input_file, output_file)

if __name__ == "__main__":
    main()
