#!/usr/bin/env python3
"""
Generate a 216 SPC water box (6x6x6 grid) for reference testing.
Box size: 1.8652 nm (density ~1000 kg/m^3 at 300K)
"""

import sys
import os

# SPC water geometry
# O at origin, H1 along +x, H2 in xy-plane
# OH bond = 0.1 nm, HOH angle = 109.47 degrees
import math

oh_bond = 0.1  # nm
angle = 109.47 * math.pi / 180.0

# H positions relative to O at origin
h1 = (oh_bond, 0.0, 0.0)
h2 = (oh_bond * math.cos(angle), oh_bond * math.sin(angle), 0.0)

# Box parameters
n_side = 6  # 6x6x6 = 216 waters
box_length = 1.8652  # nm
spacing = box_length / n_side

out_dir = os.path.dirname(os.path.abspath(__file__))

# Generate configuration
conf_lines = []
conf_lines.append("TITLE")
conf_lines.append("216 SPC water molecules in rectangular box")
conf_lines.append(f"Box: {box_length} x {box_length} x {box_length} nm")
conf_lines.append("END")
conf_lines.append("POSITION")
conf_lines.append("# first 24 chars ignored")

atom_num = 0
mol_num = 0
for ix in range(n_side):
    for iy in range(n_side):
        for iz in range(n_side):
            mol_num += 1
            # O position centered in grid cell
            ox = (ix + 0.5) * spacing
            oy = (iy + 0.5) * spacing
            oz = (iz + 0.5) * spacing
            
            atom_num += 1
            conf_lines.append(f"  {mol_num:3d} SOL   OW   {atom_num:6d}    {ox:.9f}    {oy:.9f}    {oz:.9f}")
            
            atom_num += 1
            hx1 = ox + h1[0]
            hy1 = oy + h1[1]
            hz1 = oz + h1[2]
            conf_lines.append(f"  {mol_num:3d} SOL  HW1   {atom_num:6d}    {hx1:.9f}    {hy1:.9f}    {hz1:.9f}")
            
            atom_num += 1
            hx2 = ox + h2[0]
            hy2 = oy + h2[1]
            hz2 = oz + h2[2]
            conf_lines.append(f"  {mol_num:3d} SOL  HW2   {atom_num:6d}    {hx2:.9f}    {hy2:.9f}    {hz2:.9f}")

conf_lines.append("END")
conf_lines.append("VELOCITY")

atom_num = 0
mol_num = 0
for ix in range(n_side):
    for iy in range(n_side):
        for iz in range(n_side):
            mol_num += 1
            atom_num += 1
            conf_lines.append(f"  {mol_num:3d} SOL   OW   {atom_num:6d}    0.000000000    0.000000000    0.000000000")
            atom_num += 1
            conf_lines.append(f"  {mol_num:3d} SOL  HW1   {atom_num:6d}    0.000000000    0.000000000    0.000000000")
            atom_num += 1
            conf_lines.append(f"  {mol_num:3d} SOL  HW2   {atom_num:6d}    0.000000000    0.000000000    0.000000000")

conf_lines.append("END")
conf_lines.append("GENBOX")
conf_lines.append("# box type (1=rectangular)")
conf_lines.append("1")
conf_lines.append("# box dimensions")
conf_lines.append(f"{box_length} {box_length} {box_length}")
conf_lines.append("90.0 90.0 90.0")
conf_lines.append("0.0 0.0 0.0")
conf_lines.append("0.0 0.0 0.0")
conf_lines.append("END")

with open(os.path.join(out_dir, "water_216_box.conf"), "w") as f:
    f.write("\n".join(conf_lines) + "\n")

# Generate topology
n_atoms = 216 * 3
topo_lines = []
topo_lines.append("TITLE")
topo_lines.append("216 SPC water molecules - topology")
topo_lines.append("END")
topo_lines.append("PHYSICALCONSTANTS")
topo_lines.append("# FPEPSI: 1.0/(4.0*PI*EPS0)")
topo_lines.append("138.9354")
topo_lines.append("# HBAR")
topo_lines.append("0.0635078")
topo_lines.append("# SPDL")
topo_lines.append("299792.458")
topo_lines.append("# BOLTZ")
topo_lines.append("0.00831441")
topo_lines.append("END")
topo_lines.append("TOPVERSION")
topo_lines.append("2.0")
topo_lines.append("END")
topo_lines.append("ATOMTYPENAME")
topo_lines.append("# NRATT: number of van der Waals atom types")
topo_lines.append("2")
topo_lines.append("# TYPE: atom type names")
topo_lines.append("OW")
topo_lines.append("H")
topo_lines.append("END")
topo_lines.append("RESNAME")
topo_lines.append(f"# NRAA2: number of residues in a solute molecule")
topo_lines.append(f"{216}")
topo_lines.append("# AANM: residue names")
for i in range(216):
    topo_lines.append("SOL")
topo_lines.append("END")
topo_lines.append("SOLUTEATOM")
topo_lines.append(f"#   NRP: number of solute atoms")
topo_lines.append(f"   {n_atoms}")
topo_lines.append("# ATNM MRES PANM IAC     MASS       CG  CGC INE")
topo_lines.append("#                                           INE14")

atom_num = 0
for mol in range(1, 217):
    atom_num += 1
    topo_lines.append(f"    {atom_num:4d}  {mol:3d}    OW   1   15.99940   -0.82000  0  2  {atom_num+1:4d}  {atom_num+2:4d}")
    topo_lines.append(f"                                              0")
    atom_num += 1
    topo_lines.append(f"    {atom_num:4d}  {mol:3d}   HW1   2    1.00800    0.41000  0  1  {atom_num+1:4d}")
    topo_lines.append(f"                                              0")
    atom_num += 1
    topo_lines.append(f"    {atom_num:4d}  {mol:3d}   HW2   2    1.00800    0.41000  1  0")
    topo_lines.append(f"                                              0")

topo_lines.append("END")
topo_lines.append("BONDSTRETCHTYPE")
topo_lines.append("#  NBTY: number of covalent bond types")
topo_lines.append("1")
topo_lines.append("#         CB         CHB         B0")
topo_lines.append("     1.87000e+07     3.74000e+05     1.00000e-01")
topo_lines.append("END")
topo_lines.append("BONDH")
topo_lines.append(f"#  NBONH: number of bonds involving H atoms in solute")
topo_lines.append(f"{216 * 2}")

atom_num = 0
for mol in range(216):
    o_idx = atom_num + 1
    h1_idx = atom_num + 2
    h2_idx = atom_num + 3
    topo_lines.append(f"    {o_idx:4d}  {h1_idx:4d}    1")
    topo_lines.append(f"    {o_idx:4d}  {h2_idx:4d}    1")
    atom_num += 3

topo_lines.append("END")
topo_lines.append("BOND")
topo_lines.append("#  NBON: number of bonds NOT involving H atoms in solute")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("BONDANGLEBENDTYPE")
topo_lines.append("#  NTTY: number of bond angle types")
topo_lines.append("1")
topo_lines.append("#         CT         CHT          T0")
topo_lines.append("     6.80000e+02     1.84130e-01     1.09470e+02")
topo_lines.append("END")
topo_lines.append("BONDANGLEH")
topo_lines.append(f"#  NTHEH: number of bond angles involving H atoms in solute")
topo_lines.append(f"{216}")

atom_num = 0
for mol in range(216):
    h1_idx = atom_num + 2
    o_idx = atom_num + 1
    h2_idx = atom_num + 3
    topo_lines.append(f"    {h1_idx:4d}  {o_idx:4d}  {h2_idx:4d}    1")
    atom_num += 3

topo_lines.append("END")
topo_lines.append("BONDANGLE")
topo_lines.append("#  NTHE: number of bond angles NOT involving H atoms in solute")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("IMPDIHEDRALTYPE")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("IMPDIHEDRALH")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("IMPDIHEDRAL")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("TORSDIHEDRALTYPE")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("DIHEDRALH")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("DIHEDRAL")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("CROSSDIHEDRALH")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("CROSSDIHEDRAL")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("LJPARAMETERS")
topo_lines.append("#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2")
topo_lines.append("3")
topo_lines.append("# IAC  JAC           C12            C6          CS12           CS6")
topo_lines.append("    1    1  2.634129e-06  2.617346e-03  2.634129e-06  2.617346e-03")
topo_lines.append("    1    2  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00")
topo_lines.append("    2    2  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00")
topo_lines.append("END")
topo_lines.append("LJEXCEPTIONS")
topo_lines.append("# NEX: number of exceptions")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("SOLVENTATOM")
topo_lines.append("#  NRAM: number of atoms per solvent molecule")
topo_lines.append("0")
topo_lines.append("END")
topo_lines.append("SOLVENTCONSTR")
topo_lines.append("#  NCONS: number of constraints")
topo_lines.append("0")
topo_lines.append("END")

with open(os.path.join(out_dir, "water_216_box.topo"), "w") as f:
    f.write("\n".join(topo_lines) + "\n")

# Generate charge group energy groups string for input file
negr_list = " ".join([str(i*3) for i in range(1, 217)])
print(f"Generated 216 water box: {n_atoms} atoms")
print(f"Box: {box_length} nm")
print(f"NEGR line: 216  {negr_list[:80]}...")
