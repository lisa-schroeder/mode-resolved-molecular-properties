import re
import numpy as np
import math
from scipy import constants
import mpmath


# method to read in the aoforce file and generate numpy 2d arrays containig the cartesian displacement
def aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce):

    # open file and read all lines in a list called "lines"
    with open(filename_aoforce, 'r') as fp:
        lines = fp.readlines()

        # find out the number of vibrations the molecule has
        no_vibrations = 0
        for line in lines:
            # check if 'mode' present on a current line
            if line.find('   mode   ') != -1:
                mode_list = re.split(r"\s+", line)
                no_vibrations = no_vibrations + len(mode_list) - 3
        # generate empty arrays to save x, y and z coordinates of each vibration
        array_x = np.zeros(shape=(no_vibrations,no_atoms))
        array_y = np.zeros(shape=(no_vibrations,no_atoms))
        array_z = np.zeros(shape=(no_vibrations,no_atoms))
        mass = np.zeros(shape=(no_vibrations))
        freq = np.zeros(shape=(no_vibrations))

        for line in lines:

            # check if string present on a current line
            if line.find('   mode   ') != -1:

                # split line in which mode is mentioned into a list containing each word seperated by a space
                mode_list = re.split(r"\s+", line)
                value_list_freq = re.split(r"\s+", lines[lines.index(line)+2])
                value_list_mass = re.split(r"\s+", lines[lines.index(line)+14+3*no_atoms])

                # save reduced mass and freq or each vibration in the arrays
                for no_column in range(no_columns_aoforce):
                    if len(mode_list) > no_column + 3:
                        if '***' not in value_list_mass[2]:
                            mass[int(mode_list[no_column+2])-1] = float(value_list_mass[no_column+2])
                        freq[int(mode_list[no_column+2])-1] = float(value_list_freq[no_column+2])
                
                
                # save the cartesian displacement of the columns into the 2d array
                for coord in range(no_atoms):

                    # split the lines containing the coords into a list containing the coords
                    value_list_x = re.split(r"\s+", lines[lines.index(line)+13+3*coord])
                    value_list_y = re.split(r"\s+", lines[lines.index(line)+14+3*coord])
                    value_list_z = re.split(r"\s+", lines[lines.index(line)+15+3*coord])
                    
                    # save coordinates of each mode in array, if column is present
                    for no_column in range(no_columns_aoforce):
                        if len(mode_list) > no_column + 3:
                            if value_list_x[0] == "":
                                addval = 1
                            else: 
                                addval = 0
                            array_x[int(mode_list[no_column+2])-1][coord] = float(value_list_x[no_column+3+addval])
                            array_y[int(mode_list[no_column+2])-1][coord] = float(value_list_y[no_column+2])
                            array_z[int(mode_list[no_column+2])-1][coord] = float(value_list_z[no_column+2])

                
    return array_x, array_y, array_z, mass, freq

# method to convert coord file into 3 arrays containing the x, y and z coordinates
def coord_in_arrays(filename_coord, no_atoms):
    # generate empty arrays to save x, y and z coordinates of each vibration
    array_x = np.zeros(shape=(no_atoms))
    array_y = np.zeros(shape=(no_atoms))
    array_z = np.zeros(shape=(no_atoms))
    atoms = []

    # open file and read all lines in a list
    with open(filename_coord, 'r') as fp:
        lines = fp.readlines()
        for line in range(no_atoms):
            value_list = re.split(r"\s+", lines[line+1])
            array_x[line] = value_list[1]
            array_y[line] = value_list[2]
            array_z[line] = value_list[3]
            atoms.append(value_list[4])
     
    # return arrays containing the x, y and z values and sort of atoms of the coord file
    return array_x, array_y, array_z, atoms


# method to get scaling factor for vibration dependent on Temperature according to J. Chem. Phys. 154, 244109 (2021)
def get_scaling_factor(freq, mass, no_vibration, sc_dir, sc_gauss, sc_tmole, T):
    # scaling factor dependent on temperature, scaling factor for EXPECTATION VALUE of distortion
    if T == 0:
        sc_factor = sc_dir*sc_gauss*math.sqrt(1/(2*freq[no_vibration]*mass[no_vibration]))
    else:
        sc_factor = sc_dir*sc_gauss*math.sqrt(1/(2*freq[no_vibration]*mass[no_vibration])*mpmath.coth(1.986e-23*freq[no_vibration]/(2*constants.k*T)))

    sc_energy = sc_tmole*freq[no_vibration]**2*mass[no_vibration]*sc_factor**2

    return sc_factor, sc_energy

# method to add the arrays of the aoforce and coord file dependent on the Temperature at which the vibration is generated
def get_coord_vib(no_vibration, no_atoms, ao_array_x, ao_array_y, ao_array_z, co_array_x, co_array_y, co_array_z, mass, freq, sc_tmole, sc_dir, T, sc_gauss):
    array_x = []
    array_y = []
    array_z = []

    # scaling factor dependent on temperature, scaling factor for EXPECTATION VALUE of distortion
    sc_factor, sc_energy = get_scaling_factor(freq, mass, (no_vibration-1), sc_dir, sc_gauss, sc_tmole, T)

    print("sc_energy = " + str(sc_energy))
    print("freq = " + str(freq[no_vibration-1]))
    print("mass = " + str(mass[no_vibration-1]))
    print("sc_factor = " + str(sc_factor))

    #ZPVE, ZPVE_factor, gauss_width = calculate_ZPVE(freq, (no_vibration-1), mass, sc_tmole)
    #print("ZPVE = " + str(ZPVE))
    #print("ZPVE_factor = " + str(ZPVE_factor))
    #print("gauss_width = " + str(gauss_width))

    # add each elements of one list to other list
    for i in range(no_atoms):
        array_x.append(str(sc_factor*float(ao_array_x[no_vibration-1][i]) + float(co_array_x[i])))
        array_y.append(str(sc_factor*float(ao_array_y[no_vibration-1][i]) + float(co_array_y[i])))
        array_z.append(str(sc_factor*float(ao_array_z[no_vibration-1][i]) + float(co_array_z[i])))

    return array_x, array_y, array_z, sc_energy


# method to add the arrays of the aoforce and coord file dependent on the energy which is pumped into the vibration
def get_coord_vib_energy(no_vibration, no_atoms, sc_energy, ao_array_x, ao_array_y, ao_array_z, co_array_x, co_array_y, co_array_z, mass, freq, sc_tmole, sc_dir):
    array_x = []
    array_y = []
    array_z = []

    # calcuate factor needed for scaling
    sc_factor = sc_dir*math.sqrt(sc_energy/(sc_tmole*freq[no_vibration-1]**2*mass[no_vibration-1]))

    # add each elements of one list to other list
    for i in range(no_atoms):
        array_x.append(str(sc_factor*float(ao_array_x[no_vibration-1][i]) + float(co_array_x[i])))
        array_y.append(str(sc_factor*float(ao_array_y[no_vibration-1][i]) + float(co_array_y[i])))
        array_z.append(str(sc_factor*float(ao_array_z[no_vibration-1][i]) + float(co_array_z[i])))

    return array_x, array_y, array_z


# method to write the final arrays into a coord file
def write_coord_to_file_tmole(array_x, array_y, array_z, no_vibration, atoms, Temp, sc_dir):
    if sc_dir > 0:
        filename = "coord_vib_" + str(no_vibration-6) + "_" + str(Temp)
    if sc_dir < 0:
        filename = "coord_vib_" + str(no_vibration-6) + "_m" + str(Temp)

    file1 = open(filename, "w")
    file1.write("$coord\n")
    for i in range(len(array_x)):
        file1.write("   " + '%.14f'%float(array_x[i]) + "      " + '%.14f'%float(array_y[i]) + "      " + '%.14f'%float(array_z[i]) + "  " + atoms[i] + "\n")
    file1.write("$end")
    file1.close()


# method to write the final arrays into a txt file
def write_coord_to_gaussian(array_x, array_y, array_z, no_vibration, atoms, Temp, sc_dir, sc, name_molecule):
    if sc_dir > 0:
        filename = "vib_" + str(no_vibration) + "_" + str(sc) + ".com"
    if sc_dir < 0:
        filename = "vib_" + str(no_vibration) + "_m" + str(sc) + ".com"

    file1 = open(filename, "w")
    file1.write("%nproc=36\n%mem=144GB\n%chk=" + name_molecule + "_vib_" + str(no_vibration) + "_sc_" + str(sc) + ".chk\n#p sp def2svp b3lyp\n\n" + name_molecule + " neutral S0\n\n0 1\n")
    for i in range(len(array_x)):
        file1.write(atoms[i] + "\t" + '%.14f'%float(array_x[i]) + "\t" + '%.14f'%float(array_y[i]) + "\t" + '%.14f'%float(array_z[i]) + "\n")
    file1.write("\n\n")
    file1.close() 


# method to write the final arrays into a xyz file
def write_coord_to_xyz(array_x, array_y, array_z, no_vibration, atoms, sc_factor, sc_dir):
    if sc_dir > 0:
        filename = "coord_vib_" + str(no_vibration) + "_" + str(sc_factor) + ".xyz"
    if sc_dir < 0:
        filename = "coord_vib_" + str(no_vibration) + "_m" + str(sc_factor) + ".xyz"

    file1 = open(filename, "w")
    file1.write(str(len(array_x)) + "\n\n")
    for i in range(len(array_x)):
        file1.write("\t" + atoms[i] + "\t" +  '%.14f'%float(array_x[i]) + "\t" + '%.14f'%float(array_y[i]) + "\t" + '%.14f'%float(array_z[i]) + "\n")
    file1.close() 


# method to change an array containing atomic numbers to an array containing the symbols of the elements
def get_elements_gaussian(atoms):
    new_atoms = ["0"]*len(atoms)
    for atom in range(len(atoms)):
        if int(atoms[atom]) == 1:
            new_atoms[atom] = "h"
        if int(atoms[atom]) == 6:
            new_atoms[atom] = "c"
        if int(atoms[atom]) == 8:
            new_atoms[atom] = "o"
        if int(atoms[atom]) == 30:
            new_atoms[atom] = "zn"
    return new_atoms

# method to calculate ZPVE in eV, the scaling factor belonging to that energy and the gauss_width
def calculate_ZPVE(freq, no_vibration, mass, sc_tmole):

    ZPVE = 1/2*freq[no_vibration]*0.000124
    ZPVE_factor = math.sqrt(ZPVE/(sc_tmole*freq[no_vibration]**2*mass[no_vibration]))
    gauss_width = ZPVE_factor*np.sqrt(freq[no_vibration]*mass[no_vibration])
    return ZPVE, ZPVE_factor, gauss_width
   

# generate coord file for the vibration for every Temperature defined at the beginning
# generate a file containing the energy which is pumped into each vibration
def generate_coords_tmole(Temp, no_atoms, filename_aoforce, filename_coord, no_columns_aoforce, no_vibration, scale_dir, sc_tmole, sc_gauss, numbering):
    # array to save the energy pumped into each vibration at each temperature
    energy_summary = np.zeros(shape=(len(Temp), no_atoms*3-6))

    ao_array_x, ao_array_y, ao_array_z, mass, freq = aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce)
    print("aoforce read")
    co_array_x, co_array_y, co_array_z, atoms = coord_in_arrays(filename_coord, no_atoms)
    print("coord read")

    # generating the structures for each vibration and temperature
    for no_vib in no_vibration:
        if numbering == "Turbomole":
            actual_vib = no_vib + 6
        if numbering == "Gaussian":
            actual_vib = no_vib 
        for T in range(len(Temp)):
            for sc_dir in scale_dir:
                array_x, array_y, array_z, sc_energy = get_coord_vib(no_vib, no_atoms, ao_array_x, ao_array_y, ao_array_z, co_array_x, co_array_y, co_array_z, mass, freq, sc_tmole, sc_dir, Temp[T], sc_gauss)
                energy_summary[T][no_vib-7] = sc_energy
                write_coord_to_file_tmole(array_x, array_y, array_z, actual_vib, atoms, Temp[T], sc_dir)

    # save file containing the summary of pumped in energies
    file1 = open("energy_pumped_in.txt", "w")

    temp_line = "\t"
    for T in Temp:
        temp_line = temp_line + str(T) + " K\t"
    temp_line = temp_line + "\n"
    file1.write(temp_line)

    for vib in range(no_atoms*3-6):
        if numbering == "Turbomole":
            actual_vib = vib + 7
        if numbering == "Gaussian":
            actual_vib = vib + 1
        file_line = "vib " + str(actual_vib) + "\t"
        for T in range(len(Temp)):
            file_line = file_line + str(T) + " " + str(energy_summary[T][vib]) + "\t"
        file_line = file_line + "\n"
        file1.write(file_line)

    file1.write("$end")
    file1.close() 


# generate coord file for the vibration for every energy defined at the beginning
def generate_coords_energy(no_atoms, filename_aoforce, filename_coord, no_columns_aoforce, no_vibration, scale_dir, sc_tmole, scale_energy, numbering):
    # generate the coords for every vibration and scaling factor
    for no_vib in no_vibration:
        if numbering == "Turbomole":
            actual_vib = no_vib + 6
        if numbering == "Gaussian":
            actual_vib = no_vib 
        for sc_energy in scale_energy:
            for sc_dir in scale_dir:

                ao_array_x, ao_array_y, ao_array_z, mass, freq = aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce)
                co_array_x, co_array_y, co_array_z, atoms = coord_in_arrays(filename_coord, no_atoms)
        
                array_x, array_y, array_z = get_coord_vib_energy(no_vib, no_atoms, sc_energy, ao_array_x, ao_array_y, ao_array_z, co_array_x, co_array_y, co_array_z, mass, freq, sc_tmole, sc_dir)

                write_coord_to_file_tmole(array_x, array_y, array_z, actual_vib, atoms, sc_energy, sc_dir)


# method to creat xyz file from coord file
def coord_to_xyz(filename_coord, no_atoms):
    array_x, array_y, array_z, atoms = coord_in_arrays(filename_coord, no_atoms)
    for i in range(len(array_x)):
        array_x[i] = 0.52918*array_x[i]
        array_y[i] = 0.52918*array_y[i]
        array_z[i] = 0.52918*array_z[i]

    file1 = open((filename_coord + ".xyz"), "w")
    file1.write(str(len(array_x)) + "\n\n")
    for i in range(len(array_x)):
        file1.write("\t" + atoms[i] + "\t" +  '%.14f'%float(array_x[i]) + "\t" + '%.14f'%float(array_y[i]) + "\t" + '%.14f'%float(array_z[i]) + "\n")
    file1.close() 

