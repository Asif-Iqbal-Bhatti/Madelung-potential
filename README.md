# Madelung Potential Calculator

A Python script for computing the Madelung potential (MP).

Created: 01/03/2020

Usage: python3 argv.sys[0]

The script automatically reads the CONTCAR file and requires the partial charge to be provided in the last column of the CONTCAR file.

Disclaimer: This script comes with absolutely no warranty, neither explicit nor implied. It operates by reading the CONTCAR file and computing the Madelung potential. It is strongly recommended to utilize the final atomic positions post-relaxation, along with associated charges.

This straightforward script assumes atomic charges for each atom. In practical scenarios, these charges can be approximated using charge analysis methods from DFT codes such as Mulliken charge or Bader charge analysis. By appending the charges to the last column in the CONTCAR file, the script can be adapted to read these charges for each atom and calculate the MP accordingly.
