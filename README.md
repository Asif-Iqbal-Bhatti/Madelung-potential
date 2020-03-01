# Madelung potential
Simple script to calculate the Madelung potential (MP)

Author			: Asif Iqbal
Created on	: 01/03/2020

USAGE				: python3 argv.sys[0] 

It will read the CONTCAR automatically.
NB: NO warranty guarranteed whatsoever or even implied. The script reads
CONTCAR file and compute the Madelung potential. 
It is better to read the final position after relaxation and associated charges.

This is a simple script which assumes the charges with atomic charges on the atom.
In real case it can be approximated with charges obtained from DFT code such as
Mulliken charge or Bader charge analysis. 
The charges can be appended to the last line in the CONTCAR file and the script 
can be modified to read the charges on each atom and calculate the MP.
