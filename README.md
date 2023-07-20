# Madelung potential
A simple script to calculate the Madelung potential (MP)

Created on: 01/03/2020

USAGE: python3 argv.sys[0] 

It will read the CONTCAR automatically. Supply the partial charge to the last column of the CONTCAR file.

NB: NO warranty guaranteed whatsoever, not even implied. The script reads
CONTCAR file and compute the Madelung potential.
It is better to read the final position after the relaxation and associated charges.

This simple script assumes the charges with atomic charges on the atom.
In real cases, it can be approximated with charges obtained from the DFT code such as
Mulliken charge or Bader charge analysis. 
The charges can be appended to the last column in the CONTCAR file and the script 
can be modified to read the charges on each atom and calculate the MP.
