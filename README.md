# CHO_screening
This algorithm is designed to search for carbohydrate binding TRP residues in a given protein PDB file. This was developed by Provart lab from University of Toronto.

Instructions:
download the parameters files and the script cho_screening.pl
prepare: setup the pdb file path in cho_screening.pl
run: perl cho_screening.pl pdb_name
output: if there are binding TRP residues in the given PDB file, those residues will be printed out.


example:
perl cho_screening.pl test

output:
TRP A 254
TRP A 311
TRP B 254
TRP B 311
