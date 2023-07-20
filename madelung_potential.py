#!/usr/bin/env python3.9

'''
#=======================================================================
#Simple script to calculate the Madelung potential (MP)
#
#Author      : Asif Iqbal -> @AIB_EM
#Created on  : 14/05/2020
#
#USAGE : python3 argv.sys[0] It will read the CONTCAR automatically.
#NB: NO warranty guarranteed whatsoever or even implied. The script reads
#CONTCAR file and compute the Madelung potential. It is better to read the final
#position after relaxation and associated charges.
#---
#This is a simple script which assumes the charges on the atom with atomic charges.
#In real case it should be approximated with Mulliken charge or Bader charge analysis. 
#The charges can be appended to the 4th column in the CONTCAR file and the script 
#can be modified to read the charges on each atom and calculate the MP.
#=======================================================================
'''

import sys, os
import numpy as np
import warnings

# Constants
e = 1.602176634e-19  # Elementary charge (Coulombs)
d = 1.11265e-10 # 4πϵo  is 1.11265x10-10 C2/(Jm).

contcar_file = 'POSCAR'

if not os.path.exists(contcar_file):
    print('ERROR: CONTCAR does not exist')
    sys.exit(1)

print('Reading CONTCAR ...\n')
with open(contcar_file, 'r') as f:
    lines_contcar = f.readlines()

for line in lines_contcar:
    if line.strip().lower() in ["direct", "d"]:
        warnings.warn('The file should be in Cartesian coordinates. Convert to Cartesian before proceeding.')
        sys.exit(1)
        
#-------------------------------------------------------------------
#                       Defining atomic numbers 
#-------------------------------------------------------------------
# THE CHEMICAL SYMBOL OF ELEMENTS IN THE PERIODIC TABLE, 
# EXTRACTED FROM VESTA CONFIGURATION FILE.
atomic_name = [
    "H" , "D", "He" , "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg", "Al",
    "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co",
    "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I" ,
    "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" ,
    "Np", "Pu", "Am", "XX"]

atomic_number = [
    1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
    38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
    56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91,
    92, 93, 94, 95, 96]

atomic_charge = [
    1, 1, 1, 2, -3, 4, -3, -2, -1, 0, 1, 2, 3, 4, 5, -2, -1, 0, 1, 2,
    3, 4, 2, 2, 2, 2, 2, 2, 2, 1, 2, 3, -4, -3, -2, -1, 0, 1,
    1, 2, 3, 4, 3, 3, 6, 3, 4, 2, 1, 2, 3, 2, 3, 2, -3, -2,
    -1, 0, +1, +2, +3, 3, +3, 3, +3, 3, 3, +3, 3, 3, 3, 3, 4, 5,
    6, +3, 3, 3, 2, 1, 1, 1, 2, 3, 2, 0, 0, 0, 2, 3, 4, 5,
    +3, 0, 0, 0, 0]

atomic_mass = {
'Ru': 102.91, 'Pd': 106.4, 'Pt': 195.09, 'Ni': 58.71, 'Mg': 24.305, 
'Na': 22.99, 'Nb': 92.906, 'Db': 268.0, 'Ne': 20.179, 'Li': 6.941, 'Pb': 207.2, 
'Re': 128.207, 'Tl': 204.37, 'B': 10.81, 'Ra': 226.03, 'Rb': 85.468, 'Ti': 47.9, 
'Rn': 222.0, 'Cd': 112.41, 'Po': 209.0, 'Ta': 180.95, 'Be': 9.0122, 'Fr': 223.0, 
'Te': 127.6, 'Ba': 137.33, 'Os': 190.2, 'La': 0, 'Bh': 272.0, 'Ge': 72.59, 
'Zr': 91.224, 'Tc': 97.0, 'Fe': 55.847, 'Br': 79.904, 'Sr': 87.62, 'Hf': 178.49, 
'Hg': 200.59, 'He': 4.0026, 'C': 12.011, 'Cl': 35.453, 'Rf': 265.0, 'P': 30.974, 
'F': 18.998, 'I': 126.9, 'H': 1.0079, 'Mo': 95.94, 'v': 50.941, 'Ac': 227.0, 
'O': 15.999, 'N': 14.007, 'Kr': 83.8, 'Si': 28.086, 'Sn': 118.69, 'W': 183.84, 
'Y': 88.906, 'Sb': 121.75, 'Bi': 206.98, 'Al': 26.982, 'Sg': 271.0, 'Se': 78.96, 
'Sc': 44.955912, 'Zn': 65.38, 'Co': 58.933, 'Ag': 107.87, 'Mt': 276.0, 
'k': 39.096, 'Ir': 192.22, 'S': 32.06, 'Xe': 131.3, 'Mn': 54.938, 'As': 74.922, 
'Ar': 39.948, 'Au': 196.97, 'At': 210.0, 'Ga': 69.72, 'Hs': 227.0, 'Cs': 132.91, 
'Cr': 51.996, 'Ca': 40.08, 'Cu': 63.546, 'In': 114.82}
#-------------------------------------------------------------------

dict, at_nu = {}, {}
rx, ry, rz = [], [], []
qq, r, charge = [], [], []

#============== READING ATOM TYPES AND # OF ELEMENTS
atoms_types, atom_list = lines_contcar[5].split(), lines_contcar[6].split()
output = f"{'ATOMS TYPE':30s}: {list(zip(atoms_types, atom_list))}"
print(output)

#============== 
# SUMMING ATOM NUMBERS
#============== 
sum_atoms = atom_list
for i in range( len(atoms_types) ):
    dict[i] = atoms_types[i], sum_atoms[i]   
sum_atoms = sum(list(map(int, sum_atoms)))
print(f"{'TOTAL # OF ATOMS':30s}: {sum_atoms}")

#============== 
# FINDING ATOMIC NUMBER FROM A LIST &
# MULTIPLYING THE LIST WITH ATOMIC NUMBERS 
# ACCORDING TO THE POSCAR LIST FORMAT
#============== 

for j in range(len(atoms_types)):
    for k in atomic_name: 
        if dict[j][0] == k:
            element = atomic_name.index(k)
            at_nu[j] = dict[j][0], atomic_number[element]
            #print(dict[j][0],"-->", atomic_number[element])
#print (f'Charge Dictionary -> {at_nu}')

for i in range( len(atoms_types) ):
    charge.extend(at_nu[i][1] for _ in range( int(atom_list[i]) ))
print(f"{'Charge Multiplicity':30s}: {charge}")

#============== INDEXING POSITIONS
for i in lines_contcar:
    if any(s in i for s in ["Cartesian", "C", "cartesian"]):
        lp = lines_contcar.index(i)

coord_q = np.array([list(map(float, lines_contcar[lp+1+i].split())) for i in range(sum_atoms)])
rx, ry, rz, qq = coord_q.T

print(f"{'='*50:<30}")
print(f"CHARGES IN FILE {qq}\n")

temp = 0.0
for i in range(sum_atoms):
    for j in range(sum_atoms):
        if (j > i):  # AVOIDING DOUBLE COUNTING
            r = np.sqrt( (rx[i]-rx[j])**2 + (ry[i]-ry[j])**2 + (rz[i]-rz[j])**2 )
            temp += ( qq[j]/r )
    break
    
print(temp)            
print(temp * e/d)
    
Medulung = temp * (e/d)
print(f"Medulung Potential is: {Medulung:6.3f}" )


