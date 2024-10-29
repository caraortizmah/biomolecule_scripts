import io
import sys

from pymol import editor, cmd

"""
Author:  @caraortizmah

Description: This script does a capping process on a molecule by
completing valence of missing bonds in atoms with Hydrogen atoms.

Requirements: pymol package

"""

def addh(pdb_file):
    """
    addh() caps missing atoms in a pdb file with Hydrogen atoms.
    It retuns a new pdb with the same name and an additional 
    suffix about the Hydrogen capping process.

    Args:
        addh('pdbfile.pdb')
    Return: (new pdb file)
        return 'pdbfile.capH.pdb'
    """
 
    cmd.reinitialize()
   
    try:
         cmd.load(pdb_file)
    except ValueError:
         err.write(f'There is a format error in {pdb_file} file\n')
         return 1
    except:
         err.write(f'There may be an error in {pdb_file} file\n')
         return 1
         
    cmd.h_add()
    newpdbh = pdb_file.split('.pdb')[0]+".capH.pdb"
    cmd.save(newpdbh)
    
    f.write(f'{pdb_file} succesfully capped with Hydrogen\n')

    return 0

if __name__ == "__main__":
    pdbpairaa = sys.argv[1]
    
    # list of pairaa pdb files successfully Hydrogen-capped
    f = open('list_pairaaH.txt', 'a', encoding="utf-8")
    
    #  ...in case of errors...
    err = open('err_list_pairaaH.err', 'a', encoding="utf-8")

    # list of .pdb files that will be opened to run xx()
    # function
    with open(pdbpairaa) as listpdb:
        for line in listpdb:
            line = line.splitlines()[0]

            addh(line)
            
        f.close()
    err.close()

exit()

