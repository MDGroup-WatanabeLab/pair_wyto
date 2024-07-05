import os
import numpy as np
import shutil
import argparse
import random
import itertools
import math
from decimal import Decimal, ROUND_HALF_UP

# The number of processes
mpi_num = 16

# Path to lammps
lmp_path = "lmp_mpi"

# Cut off radial of Si and O
cut_off_Si = 2.8
cut_off_Si_O = 1.8

# Lattice parameter of Si
l_a = 5.431
l_b = 5.431
l_c = 5.431

# Consider the cells around the original cell because of periodic boundary condition
"""
original cell : (x, y) = (0, 0)

  y
  ^
  |  ____________________________
  |  |        |        |        |
  |  |(-1, +1)|(0, +1) |(+1, +1)|
  |  |________|________|________|
  |  |        |        |        |
  |  |(-1, 0) | (0, 0) |(+1, 0) |
  |  |________|________|________|
  |  |        |        |        |
  |  |(-1, -1)|(0, -1) |(+1, -1)|
  |  |________|________|________|
  |
  |---------------------------> x

"""
cell = [0, 1, -1]
cell_ar = list(itertools.product(cell, repeat=2))
cell_ar.remove((0, 0))

# Make lammps input files
# in.dimer.Si : For structural optimization
# in.stable.SiO : For structural optimization
# in.md.SiO : For NVT MD 
def make_in_file():
    with open("in.dimer.Si", "w") as f:
        f.write("units        metal\n")
        f.write("boundary     p p m\n")
        f.write("atom_style   molecular\n")
        f.write("atom_modify  sort 10000 1.0\n")
        f.write("read_data    stable.data\n")
        f.write("pair_style   wyto\n")
        f.write("pair_coeff   * * SiO.wyto Si\n")
        f.write("neighbor     4.0 bin\n")
        f.write("neigh_modify every 1 delay 0 check yes\n")
        f.write("group        fix_lay molecule 255\n")
        f.write("group        mov_lay molecule != 255\n")
        f.write("fix          1 fix_lay setforce 0 0 0\n")
        f.write("min_style    cg\n")
        f.write("minimize     1e-25 1e-25 50000 100000\n")
        f.write("timestep     0.001\n")
        f.write("fix          2 all nvt temp 300 300 0.05\n")
        f.write("velocity     mov_lay create 1000 458127641 mom yes rot yes dist gaussian\n")
        f.write("run          10000\n")
        f.write("min_style    cg\n")
        f.write("minimize     1e-25 1e-25 50000 100000\n")
        f.write("dump         1 all custom 1 stable.final id type x y z mol element\n")
        f.write("dump_modify  1 element Si\n")
        f.write("dump_modify  1 sort id\n")
        f.write("run          0\n")
        f.write("undump       1\n")
    
    with open("in.stable.SiO", "w") as f:
        f.write("units        metal\n")
        f.write("boundary     p p m\n")
        f.write("atom_style   molecular\n")
        f.write("atom_modify  sort 10000 1.0\n")
        f.write("read_data    stable.data\n")
        f.write("pair_style   wyto\n")
        f.write("pair_coeff   * * SiO.wyto Si O\n")
        f.write("neighbor     4.0 bin\n")
        f.write("neigh_modify every 1 delay 0 check yes\n")
        f.write("group        fix_lay molecule 255\n")
        f.write("fix          1 fix_lay setforce 0 0 0\n")
        f.write("min_style    cg\n")
        f.write("minimize     1e-25 1e-25 50000 100000\n")
        f.write("dump         1 all custom 1 stable.final id type x y z mol element\n")
        f.write("dump_modify  1 element Si O\n")
        f.write("dump_modify  1 sort id\n")
        f.write("run          0\n")
        f.write("undump       1\n")        
    
    with open("in.md.SiO", "w") as f:
        f.write("units        metal\n")
        f.write("boundary     p p m\n")
        f.write("atom_style   molecular\n")
        f.write("atom_modify  sort 10000 1.0\n")
        f.write("read_data    md.data\n")
        f.write("pair_style   wyto\n")
        f.write("pair_coeff   * * SiO.wyto Si O\n")
        f.write("neighbor     4.0 bin\n")
        f.write("neigh_modify every 1 delay 0 check yes\n")
        f.write("group        fix_lay molecule 255\n")
        f.write("group        mov_lay molecule != 255\n")
        f.write("fix          1 all nvt temp 1000 1000 0.05\n")
        f.write("fix          2 fix_lay setforce 0 0 0\n")
        f.write("velocity     mov_lay create 1000 458127641 mom yes rot yes dist gaussian\n")
        f.write("timestep     0.0001\n")
        f.write("dump         1 all custom 100 oxi.traj id type x y z mol element\n")
        f.write("dump_modify  1 element Si O\n")
        f.write("dump_modify  1 sort id\n")
        f.write("run          20000\n")
        f.write("undump       1\n")
        f.write("dump         2 all custom 1 md.final id type x y z mol element\n")
        f.write("dump_modify  2 element Si O\n")
        f.write("dump_modify  2 sort id\n") 
        f.write("run          0\n")   
        f.write("undump       2")
    return

# Make Si crystal
def make_config(size):
    pos = []
    # Si atomic position (Si unit cell consists of four layers)
    atom_pos1 = [[0.0, 0.0, 0.0], [1/2, 1/2, 0.0]]
    atom_pos2 = [[1/4, 1/4, 0.0], [3/4, 3/4, 0.0]]
    atom_pos3 = [[0.0, 1/2, 0.0], [1/2, 0.0, 0.0]]
    atom_pos4 = [[1/4, 3/4, 0.0], [3/4, 1/4, 0.0]]
    
    # Positions of Si (Cartesian coordinates)
    # Layers are numbered in order of proximity to the surface
    # The layer farthest from the surface is 255
    for k in range(size[2]*4):
        for j in range(size[1]):
            for i in range(size[0]):
                for l in range(2):
                    if k % 4 == 0:
                        x = (atom_pos4[l][0] + i) * l_a
                        y = (atom_pos4[l][1] + j) * l_b
                        z = (atom_pos4[l][2] + k) * l_c / 4.0
                    elif k % 4 == 1:
                        x = (atom_pos3[l][0] + i) * l_a
                        y = (atom_pos3[l][1] + j) * l_b
                        z = (atom_pos3[l][2] + k) * l_c / 4.0
                    elif k % 4 == 2:
                        x = (atom_pos2[l][0] + i) * l_a
                        y = (atom_pos2[l][1] + j) * l_b
                        z = (atom_pos2[l][2] + k) * l_c / 4.0
                    elif k % 4 == 3:
                        x = (atom_pos1[l][0] + i) * l_a
                        y = (atom_pos1[l][1] + j) * l_b
                        z = (atom_pos1[l][2] + k) * l_c / 4.0
                    if k == 0:
                        pos.append([str(x), str(y), str(z), "255"])
                    else:
                        pos.append([str(x), str(y), str(z), str(size[2] * 4 - k)])
    return pos

# Make data file (lammps input file)
def make_data_file(position, lattice, atom_num, atom_type, f_name):
    # Write data file
    with open(f_name + ".data", "w") as f:
        f.write("data file for lammps input\n\n")
        f.write(" " + str(int(atom_num[0])+int(atom_num[1])) + " atoms\n")
        f.write(" " + str(len(atom_type)) + " atom types\n")
        f.write(" 0.00 " + lattice[0][0] + " xlo xhi\n")
        f.write(" 0.00 " + lattice[1][1] + " ylo yhi\n")
        f.write(" 0.00 " + lattice[2][2] + " zlo zhi\n")
        f.write("\n")
        f.write(" Masses\n\n")
        f.write("1 28.0855\n")
        if len(atom_type) == 2:
            f.write("2 15.9994\n")
        f.write("\n")
        f.write(" Atoms\n\n")
        
        # Write positions
        # Si atom : layer number + 1 + positions
        # O atom : 0 + 2 +  positions
        for i in range(len(position)):
            f.write(str(i+1) + " ")
            if position[i][3] == "255":
                f.write("255 1")
            elif position[i][3] == "0":
                f.write("0 2")
            else:
                f.write(position[i][3] + " 1")            
            for j in range(3):
                f.write(" " + position[i][j])
            f.write("\n")
    return

# Calculate the distance between two atoms
def bond_length(position1, position2):
    a = np.array([float(position1[0]), float(position1[1]), float(position1[2])])
    b = np.array([float(position2[0]), float(position2[1]), float(position2[2])])
    return np.linalg.norm(a - b)
    

# Calculate the position of the midpoint of two atoms
def insert_o(position1, position2, lattice):
    x = (float(position1[0]) + float(position2[0])) * 0.5
    y = (float(position1[1]) + float(position2[1])) * 0.5
    z = (float(position1[2]) + float(position2[2])) * 0.5

    # If the position of the midpoint is outside the cell, move it into the cell
    if x < 0:
        x += float(lattice[0][0])
    elif x > float(lattice[0][0]):
        x -= float(lattice[0][0])
    
    if y < 0:
        y += float(lattice[1][1])
    elif y > float(lattice[1][1]):
        y -= float(lattice[1][1])
    
    if z < 0:
        z += float(lattice[2][2])
    elif z > float(lattice[2][2]):
        z -= float(lattice[2][2])
    
    x = Decimal(str(x)).quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP)
    y = Decimal(str(y)).quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP)
    z = Decimal(str(z)).quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP)

    return [str(x), str(y), str(z)]

# Find all positions where O atom can be inserted
def check_insert(position, position_pbc, layer, lattice):
    pos = []
    for i in range(len(position)):
        for j in range(i+1, len(position)):
            # Original cell
            # Focus on Si atoms in two layer
            # If the distance between two atoms is within the cut_off_si, get the position of the midpoint of two atoms
            if layer <= int(position[i][3]) <= layer+1 and layer <= int(position[j][3]) <= layer+1:
                if not (int(position[i][3]) == layer+1 and int(position[j][3]) == layer+1):
                    if bond_length(position[i], position[j]) < cut_off_Si:
                        tmp_o = insert_o(position[i], position[j], lattice) + ["0"]
                        if tmp_o not in pos:
                            pos.append(tmp_o)

        for j in range(len(position_pbc)):
            # Cells around the original cell (consider periodic boundary condition)
            # Focus on Si atoms in two layer
            # If the distance between two atoms is within the cut_off_si, get the positions of the midpoint of the two atom
            if layer <= int(position[i][3]) <= layer+1 and layer <= int(position_pbc[j][3]) <= layer+1:
                if not (int(position[i][3]) == layer+1 and int(position_pbc[j][3]) == layer+1):
                    if bond_length(position[i], position_pbc[j]) < cut_off_Si:
                        tmp_o = insert_o(position[i], position_pbc[j], lattice) + ["0"]
                        if tmp_o not in pos:
                            pos.append(tmp_o)

    return pos

# Shift the x, y coordinates of an atom by the lattice constant
def pos_adjust(pos, pbc, lattice):
    x = float(pos[0]) + pbc[0]*float(lattice[0][0]) + pbc[1]*float(lattice[1][0])
    y = float(pos[1]) + pbc[0]*float(lattice[0][1]) + pbc[1]*float(lattice[1][1])
    z = pos[2]

    return [str(x), str(y), z]

# Read final file and get lattice vectors and positions of Si and O atoms
def read_final_file(f_name):
    with open(f_name, "r") as f:
        line = f.readlines()
    pos_si = []
    pos_o = []

    # Get lattice vectors
    x = line[5].split()
    y = line[6].split()
    z = line[7].split()
    l_vec = [[str(float(x[1]) - float(x[0])), "0.0", "0.0"], ["0.0", str(float(y[1]) - float(y[0])), "0.0"], ["0.0", "0.0", str(float(z[1]) - float(z[0]))]]    
    
    # Get positions of Si and O atoms
    for i in range(9, len(line)):
        _, atom, *pos = line[i].split()
        if atom == "1":
            pos_si.append(pos)
        else:
            pos_o.append(pos)

    return l_vec, pos_si, pos_o

# Input cell size
parser = argparse.ArgumentParser(description="input two argument")
parser.add_argument("-n", required=True, help="Please input the cell size \"x y z\"")
args = parser.parse_args()
size = args.n.split()

# Non-numeric and integer numbers less than or equal to 0 are not accepted
if len(size) == 3:
    for i in range(len(size)):
        if size[i].isdecimal() and int(size[i]) > 0:
            size[i] = int(size[i])
        else:
            flag = 0
            print("Please input number more than 0")
            exit()
else:
    print("Please input three number")
    exit()

# Lattice vectors
lattice = [[str(l_a * size[0]), "0.0", "0.0"], ["0.0", str(l_b * size[1]), "0.0"], ["0.0", "0.0", str(l_c * size[2])]]

# Positions of Si
position_si = make_config(size)

# Make vacuum region on top of Si crystal
position_o = []
lattice[2][2] = str(float(lattice[2][2]) * 2.0)

# Number of atoms and atomic species
atom_num = [len(position_si), len(position_o)]
atom_type = ["Si"]

# Make lammps input files
make_in_file()

# Make data file
make_data_file(position_si+position_o, lattice, atom_num, atom_type, "stable")

# Structure optimization and NVT MD (make dimer)
# Copy input and output files
os.system("mpirun -np " + str(mpi_num) + " " + lmp_path + " -in in.dimer.Si")
os.mkdir("si_init")
shutil.copy("stable.data", "si_init/stable.data")
shutil.copy("stable.final", "si_init/stable.final")
shutil.copy("log.lammps", "si_init/log.lammps")

# If set insert_rate <= 0.0 or insert_rate >= 1.0, insert all O atoms at once
insert_rate = 0.1
atom_type = ["Si", "O"] 

# Oxidiation process
os.mkdir("oxidation")
for oxi_lay in range(1, int(size[2]*4 / 2) + 1):
    os.mkdir("oxidation/layer" + str(oxi_lay))
    break_flag = 0
    oxidation_count = 0

    # Oxidation of each layer
    while break_flag == 0:
        oxidation_count += 1

        # Read final file after structure optimization and get lattice vectors and positions of Si and O atoms
        lattice, position_si, position_o = read_final_file("stable.final")

        # Get the position of Si in the two layers
        position_si_target = []
        for i in range(len(position_si)):
            if int(position_si[i][3]) == oxi_lay or int(position_si[i][3]) == oxi_lay+1:
                position_si_target.append(position_si[i])

        # Get positions of Si in the cells around the original cell
        position_si_pbc = []   
        for pbc in range(len(cell_ar)):
            for i in range(len(position_si_target)):
                tmp = pos_adjust(position_si_target[i], cell_ar[pbc], lattice)
                position_si_pbc.append(tmp + [position_si_target[i][3]])
            
        # Stores positions of insertable O atoms
        position_o_tmp = check_insert(position_si_target, position_si_pbc, oxi_lay, lattice)
        if (len(position_o_tmp)) == 0:
            break

        # Get positions of O in the cells around the original cell
        position_o_pbc = []   
        for pbc in range(len(cell_ar)):
            for i in range(len(position_o)):
                position_o_pbc.append(pos_adjust(position_o[i], cell_ar[pbc], lattice))

        # Check coordination number of Si
        # if it exceeds 4, remove the position of O from position_o_tmp         
        for i in range(len(position_si_target)):
            count = 0
            re_list_tmp = []
            tmp = []
            for j in range(len(position_o)):
                if bond_length(position_si_target[i], position_o[j]) < cut_off_Si_O:
                    count += 1
            for j in range(len(position_o_pbc)):
                if bond_length(position_si_target[i], position_o_pbc[j]) < cut_off_Si_O:
                    count += 1
            for j in range(len(position_o_tmp)):
                if bond_length(position_si_target[i], position_o_tmp[j]) < cut_off_Si_O:
                    count += 1
                    re_list_tmp.append(position_o_tmp[j])
            
            if count >= 5:
                tmp = random.sample(re_list_tmp, count-4)
                for j in range(len(tmp)):
                    position_o_tmp.remove(tmp[j])

        if (len(position_o_tmp)) == 0:
            break

        # Insert O atoms
        if 0.0 < insert_rate < 1.0:
            insert_num = math.ceil(len(position_o_tmp) * insert_rate)
            try:
                position_o = position_o + random.sample(position_o_tmp, insert_num)
            except ValueError:
                position_o = position_o + position_o_tmp
                break_flag = 1
        else:
            position_o = position_o + position_o_tmp
            break_flag = 1

        # Get the number of Si and O 
        atom_num = [len(position_si), len(position_o)] 

        # Make data file
        make_data_file(position_si+position_o, lattice, atom_num, atom_type, "stable")

        # Structure optimization
        # Copy input and output files
        os.system("mpirun -np " + str(mpi_num) + " " + lmp_path + " -in in.stable.SiO")
        os.mkdir("oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_1st")
        shutil.copy("stable.data", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_1st/stable.data")
        shutil.copy("stable.final", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_1st/stable.final")
        shutil.copy("log.lammps", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_1st/log.lammps")

        # Read final file after structural optimization and get lattice vectors and positions of Si and O atoms
        lattice, position_si, position_o = read_final_file("stable.final")

        # Make data file
        make_data_file(position_si+position_o, lattice, atom_num, atom_type, "md")

        # NVT MD 
        # Copy input and output files
        os.system("mpirun -np " + str(mpi_num) + " " + lmp_path + " -in in.md.SiO")
        os.mkdir("oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_md")
        shutil.copy("stable.data", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_md/stable.data")
        shutil.copy("md.final", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_md/md.final")
        shutil.copy("log.lammps", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_md/log.lammps")
        shutil.copy("oxi.traj", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_md/oxi.traj")

        # Read final file after structural optimization to obtain lattice vectors and positions of Si and O atoms
        lattice, position_si, position_o = read_final_file("md.final")

        # Make data file
        make_data_file(position_si+position_o, lattice, atom_num, atom_type, "stable")

        # Structure optimization
        # Copy input and output files
        os.system("mpirun -np " + str(mpi_num) + " " + lmp_path + " -in in.stable.SiO")
        os.mkdir("oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_2nd")
        shutil.copy("stable.data", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_2nd/stable.data")
        shutil.copy("stable.final", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_2nd/stable.final")
        shutil.copy("log.lammps", "oxidation/layer" + str(oxi_lay) + "/" + str(oxidation_count) + "_stable_2nd/log.lammps")
