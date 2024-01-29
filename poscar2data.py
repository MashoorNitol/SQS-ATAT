#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POSCAR to LAMMPS Data Converter

This code converts a VASP POSCAR file to a LAMMPS data file in the atomic format.
During the conversion, the code excludes the mass section in the resulting LAMMPS data file.
The LAMMPS data file is created with the same name as the POSCAR file, but with a different extension.

@author: mash
"""

import numpy as np
import os

class POSCARConverter:
    def __init__(self, poscar_file):
        self.poscar_file = poscar_file
        self.atom_data = {}
        self.positions = []
        self.basevect = None
        self.total_sum = 0
        self.scaling = 0.0
        self.total_atoms = []
        self.atom_types = []
        self.tot_atom_types = 0
        self.num_atoms = []
        self.positions_start = 0
        self.atom_type_index = 1
        self.all_positions = None
        self.coord_type = ''
        self.cart_coords = None
        self.lx = 0.0
        self.xy = 0.0
        self.xz = 0.0
        self.ly = 0.0
        self.yz = 0.0
        self.lz = 0.0

    def convert(self, lammps_file):
        self._read_poscar()
        self._compute_cartesian_coordinates()
        self._write_lammps_data(lammps_file)

    def _read_poscar(self):
        with open(self.poscar_file, 'r') as f:
            lines = f.readlines()
            self.total_atoms = lines[6].strip().split()
            self.total_sum = sum(int(atom) for atom in self.total_atoms)
            self.scaling = float(lines[1].strip())
            self.basevect = self.scaling * np.array([list(map(float, lines[i].split())) for i in range(2, 5)])
            self.atom_types = lines[5].split()
            self.tot_atom_types = len(self.atom_types)
            self.num_atoms = [int(x) for x in self.total_atoms]
            self.positions_start = 8
            for i in range(len(self.atom_types)):
                atom_type = "type " + str(self.atom_type_index)
                self.atom_type_index += 1
                num_atom = self.num_atoms[i]
                atom_positions = []
                for j in range(self.positions_start, self.positions_start + num_atom):
                    position = lines[j].split()
                    atom_positions.append([float(position[0]), float(position[1]), float(position[2])])
                    self.positions.append([float(position[0]), float(position[1]), float(position[2])])
                self.atom_data[atom_type] = np.array(atom_positions)
                self.positions_start += num_atom

        self.all_positions = np.array(self.positions)
        self.coord_type = str(lines[7]).strip()

    def _compute_cartesian_coordinates(self):
        if self.coord_type in ('Direct', 'D', 'direct'):
            self.cart_coords = np.matmul(self.all_positions, self.basevect)
        elif self.coord_type in ('Cartesian', 'C', 'cartesian'):
            self.cart_coords = self.all_positions

        self.lx = np.linalg.norm(self.basevect[0])
        self.xy = np.dot(self.basevect[0], self.basevect[1]) / self.lx
        self.xz = np.dot(self.basevect[0], self.basevect[2]) / self.lx
        self.ly = np.sqrt(np.linalg.norm(self.basevect[1]) ** 2 - self.xy ** 2)
        self.yz = (np.dot(self.basevect[1], self.basevect[2]) - self.xy * self.xz) / self.ly
        self.lz = np.sqrt(np.linalg.norm(self.basevect[2]) ** 2 - self.xz ** 2 - self.yz ** 2)
        
    def _write_lammps_data(self, lammps_file):
        with open(lammps_file, 'w') as f:
            f.write("# LAMMPS data file\n\n")
            f.write('{:d} atoms\n'.format(self.total_sum))
            f.write('{:d} atom types\n\n'.format(self.tot_atom_types))
            f.write('0.0 {:.16f} xlo xhi\n'.format(self.lx))
            f.write('0.0 {:.16f} ylo yhi\n'.format(self.ly))
            f.write('0.0 {:.16f} zlo zhi\n'.format(self.lz))
            f.write('{:.16f} {:.16f} {:.16f} xy xz yz\n\n'.format(self.xy, self.xz, self.yz))
            f.write('Atoms # atomic\n\n')
            index = 1
            for atom_type, positions in self.atom_data.items():
                for position in positions:
                    formatted_coords = ' '.join('{:.16f}'.format(coord) for coord in self.cart_coords[index - 1])
                    line = '{:d} {:s} {}'.format(index, atom_type.split()[1], formatted_coords)
                    f.write(line + "\n")
                    index += 1

def convert_poscar_to_lammps(poscar_file):
    lammps_file = poscar_file.replace('-POSCAR', '.lmp')
    converter = POSCARConverter(poscar_file)
    converter.convert(lammps_file)

files = os.listdir()
for file in files:
    if file.endswith('-POSCAR'):
        convert_poscar_to_lammps(file)

# poscar_file = 'Ti3Al.POSCAR'
# lammps_file = poscar_file.replace('.poscar', '.data')
# 
# converter = POSCARConverter(poscar_file)
# converter.convert(lammps_file)
# 
