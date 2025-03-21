### Python scripts for molecule manipulation

`cappingH.py` caps any molecule with Hydrogen atoms where the Lewis structure show an incomplete
valence. That is valid for organic molecules following the classic and standard valence of atoms
in molecules.

`convert_xyztopdb.py` convert any molecule in xyz file into a pdb file using openbabel.
Run it as:
$ convert_xyztopdb.py input.xyz output.pdb

where you provide the `input.xyz` and make sure that openbabel can be called by python.

You can run the same script for a list of xyz files using `xyztopdb.sh` as
$ ./xyztopdb.sh

the `xyztopdb.sh` script recognizes the `convert_xyztopdb.py` script and run for all the xyz files
placed in the same folder.
