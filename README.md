# Madelung Potential Calculator

A Python script for computing the Madelung potential (MP) and total electrostatic energy of crystal systems from VASP geometry files.

**Created**  :: 14/05/2020  
**Updated**  :: 08/04/2026  
**Author**   :: Asif Iqbal ([@AIB_EM](https://github.com/AIB_EM))

---

## Usage
1. Ensure your structure file is named `CONTCAR` or `POSCAR`.
2. Append your partial charges as a **4th column** in the coordinate section.


## Key Features
Automatic Coordinate Conversion: Automatically detects and converts "Direct" (fractional) coordinates to "Cartesian" for distance calculations.
Periodic Boundary Conditions (PBC): Unlike simple calculators that only sum atoms within one unit cell, this script utilizes a Supercell Expansion (shell summation). It calculates the interaction of each atom with all atoms in a 7×7×7  block of unit cells (controllable via shell_limit) to ensure a more converged Madelung sum.
Scientific Precision: Uses the CODATA 2018 physical constants. It converts results from internal units (e/A˚e/A˚) to Volts (V) for site potentials and electron-Volts (eV) for total energy.
Neutrality Check: Automatically verifies if the system is charge-neutral, as Madelung sums for non-neutral systems are physically divergent.
