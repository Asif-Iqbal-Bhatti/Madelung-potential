#!/usr/bin/env python3

"""
Madelung Potential Calculator for VASP structures.

This script parses a VASP POSCAR/CONTCAR file, extracts atomic positions 
and charges (from the 4th column), and computes the Madelung potential 
at each atomic site using a direct-sum supercell approach.
"""

import numpy as np
import os
import sys

# Physical Constants (CODATA 2018)
EPS0 = 8.8541878128e-12    # Vacuum permittivity (F/m)
E_CHG = 1.602176634e-19    # Elementary charge (C)
A2M = 1e-10                # Angstrom to meter
KE = 1 / (4 * np.pi * EPS0) # Coulomb constant (~8.987e9 N·m²/C²)

def read_vasp(filename):
    """Parses VASP geometry files and extracts lattice and charges.

    Args:
        filename (str): Path to the POSCAR or CONTCAR file.

    Returns:
        tuple: (lattice, coords, charges)
            - lattice (np.array): 3x3 array of lattice vectors.
            - coords (np.array): Nx3 array of Cartesian coordinates in Angstroms.
            - charges (np.array): N-length array of charges from the 4th column.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file format is incorrect.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Structure file '{filename}' not found.")

    with open(filename, 'r') as f:
        lines = f.readlines()

    try:
        # Lattice and Scaling
        scale = float(lines[1].strip())
        lattice = np.array([l.split() for l in lines[2:5]], dtype=float) * scale

        # Species handling
        counts = list(map(int, lines[6].split()))
        total_atoms = sum(counts)

        # Detect start of coordinates (Skip 'Selective Dynamics' if present)
        offset = 8 if "selective" in lines[7].lower() else 7
        coord_type = lines[offset].strip().lower()
        
        coords = []
        charges = []
        for i in range(total_atoms):
            parts = lines[offset + 1 + i].split()
            coords.append(parts[:3])
            # Assign charge from 4th column, default to 0.0 if missing
            charges.append(float(parts[3]) if len(parts) >= 4 else 0.0)

        coords = np.array(coords, dtype=float)
        charges = np.array(charges, dtype=float)

        # Convert Direct (fractional) to Cartesian (Angstrom)
        if coord_type.startswith('d'):
            coords = coords @ lattice

        return lattice, coords, charges
    
    except Exception as e:
        raise ValueError(f"Failed to parse VASP file: {e}")

def calculate_madelung(lattice, coords, charges, shell_limit=3):
    """Computes the electrostatic potential at each atomic site.

    This function performs a direct summation of q/r over periodic images.
    
    Args:
        lattice (np.array): 3x3 lattice vectors.
        coords (np.array): Nx3 atomic positions.
        charges (np.array): N-length array of charges.
        shell_limit (int): Number of unit cells to expand in each direction.
            A value of 3 creates a (2*3+1)^3 = 343 cell supercell.

    Returns:
        np.array: Madelung potential at each site in Volts (V).
    """
    n_atoms = len(coords)
    potentials = np.zeros(n_atoms)

    # 1. Verification: Charge Neutrality (Required by Literature)
    net_charge = np.sum(charges)
    if abs(net_charge) > 1e-3:
        print(f"Warning: System is not neutral (Net: {net_charge:.3f}). "
              "Madelung sums for non-neutral systems are physically divergent.")

    # 2. Generate Periodic Image Offsets
    r_ranges = range(-shell_limit, shell_limit + 1)
    offsets = [i*lattice[0] + j*lattice[1] + k*lattice[2] 
               for i in r_ranges for j in r_ranges for k in r_ranges]
    offsets = np.array(offsets)

    # 3. Double Loop Summation
    for i in range(n_atoms):
        phi_i = 0.0
        pos_i = coords[i]
        
        for image_offset in offsets:
            # Vectorized distance calculation for all atoms in this image box
            pos_j_images = coords + image_offset
            dists = np.linalg.norm(pos_j_images - pos_i, axis=1)
            
            for j in range(n_atoms):
                # Skip self-interaction
                if dists[j] < 1e-5:
                    continue
                phi_i += charges[j] / dists[j]
        
        potentials[i] = phi_i

    # Unit Conversion: (e / Angstrom) -> Volts
    # Conversion factor is approx 14.3996
    volts_conversion = (E_CHG / A2M) * KE
    return potentials * volts_conversion

def main():
    """Main execution block."""
    input_file = 'CONTCAR'
    
    try:
        lattice, coords, charges = read_vasp(input_file)
        print(f"File parsed: {len(coords)} atoms found.")
        
        # Calculate potentials
        potentials = calculate_madelung(lattice, coords, charges, shell_limit=3)
        
        # Calculate total electrostatic energy (Madelung Energy)
        # U = 0.5 * sum(q_i * V_i)
        total_energy_ev = 0.5 * np.sum(charges * potentials)

        # Output results
        print("\n" + "="*50)
        print(f"{'Index':<6} {'Charge (e)':<12} {'Potential (V)':<15}")
        print("-"*50)
        for i, (q, v) in enumerate(zip(charges, potentials)):
            print(f"{i+1:<6} {q:<12.3f} {v:<15.4f}")
        
        print("="*50)
        print(f"Total Madelung Energy: {total_energy_ev:.6f} eV")
        print("="*50)

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
