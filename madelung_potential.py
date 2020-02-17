#!/usr/bin/env python3

'''
#-------------------------------------------------------------------
Simple script to calculate the Medulung potential
Author: Asif Iqbal
Created on: 17/02/2020
USAGE: python3 argv.sys[0] It will read the POSCAR automatically.
NB: NO warranty guarranteed whatsoever even implied. This script is written
to read POSCAR file and compute the Madelung potential. 
#-------------------------------------------------------------------
'''

import sys, os
import numpy as np	

e = 1.60217662E10-19 # C
d = 1.11265E10-10 # 4πϵo  is 1.11265x10-10 C2/(J m).
	
if not os.path.exists('POSCAR'):
	print (' ERROR: POSCAR does not exist')
	sys.exit(0)
print('Reading POSCAR ... \n')

#-------------------------------------------------------------------
# 													Defining atomic numbers 
#-------------------------------------------------------------------
# the chemical symbol of elements in the periodic table, extracted from VESTA
# configuration file.
atomic_name = [
    "H" , "D", "He" , "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg", "Al",
    "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co",
    "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I" ,
    "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" ,
    "Np", "Pu", "Am", "XX"
]
# the atomic number of elements in the periodic table, extracted from VESTA
# configuration file.
atomic_number = [
    1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
    38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
    56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91,
    92, 93, 94, 95, 96,
]

atomic_mass = {'Ru': 102.91, 'Pd': 106.4, 'Pt': 195.09, 'Ni': 58.71, 'Mg': 24.305, 
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

x = []; y =[]; z =[]; 
rx=[]; ry =[]; rz=[]; at_nu={}
r=[]; dict={}; temp=0 ; charge=[]

f = open('POSCAR','r')
lines_poscar = f.readlines()
f.close()

atoms_types = lines_poscar[5].split()
atom_list = lines_poscar[6].split()  ### reading 7th lines for reading # of atoms
sum_atoms = atom_list
print("Atom Numbers:")
for i in range( len(atoms_types) ):
	dict[i] = atoms_types[i],sum_atoms[i]
sum_atoms = [int(i) for i in sum_atoms]
sum_atoms = sum(sum_atoms)

#--------------------- Finding atomic number from a list ---------------------

for j in range(len(atoms_types)):
	for k in atomic_name: 
		if dict[j][0] == k:
			element = atomic_name.index(k)
			at_nu[j] = dict[j][0],atomic_number[element]
			print(dict[j][0],"-->",atomic_number[element])
print (at_nu)

#-------------------------Multiplying the list with atomic numbers 
# according to the POSCAR list format

for i in range( len(atoms_types) ):
	for k in range( int(atom_list[i]) ):
		charge.append(at_nu[i][1])
print(charge)

#-------------------------------------------------------------------------

for i in lines_poscar:
	if "Direct" in i:
		lp=lines_poscar.index(i)

print("-"*80)		
for i in range(sum_atoms):
	x, y, z = lines_poscar[lp+1+i].split()
	print (x, y, z)
	x = float(x); y = float(y); z = float(z)
	rx.append(x);	ry.append(y);	rz.append(z)
	##r = [rx, ry, rz]
	
#r = np.matrix(r)	
print("-"*80)
for i in range(sum_atoms):
	for j in range(sum_atoms):
		if (j > i):  # over here avoiding double counting
			
			temp = temp + charge[j] * 1/( np.sqrt( (rx[i]-rx[j])**2 + (ry[i]-ry[j])**2 + (rz[i]-rz[j])**2 ) )
	break
Medulung = temp * (e/d)

print(Medulung) 






