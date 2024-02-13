import re
import numpy as np
import math
from scipy import constants
import mpmath


# method to read in the aoforce file and generate numpy 2d arrays containig the cartesian displacement
def dEidR_of_HOMO(filename, no_atoms, no_HOMO):

    # open file and read all lines in a list
    with open(filename, 'r') as fp:
        lines = fp.readlines()

        # generate empty arrays to save x, y and z coordinates of each vibration
        array_x = np.zeros(shape=(no_atoms))
        array_y = np.zeros(shape=(no_atoms))
        array_z = np.zeros(shape=(no_atoms))

        for line in lines:
            # check if 'mode' present on a current line
            if line.find('#') == -1:
                value_list = re.split(r"\s+", line)
            
                # save x coordinates of each atom in array
                if int(value_list[1]) == no_HOMO:
                    if int(value_list[3]) == 1:
                        array_x[int(value_list[2])-1] = float(value_list[4])
                    if int(value_list[3]) == 2:
                        array_y[int(value_list[2])-1] = float(value_list[4])
                    if int(value_list[3]) == 3:
                        array_z[int(value_list[2])-1] = float(value_list[4])

                                      
    return array_x, array_y, array_z


# method to read in the aoforce file and generate numpy 2d arrays containig the cartesian displacement
def aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce):

    # open file and read all lines in a list
    with open(filename_aoforce, 'r') as fp:
        lines = fp.readlines()

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
        array_mass = np.zeros(shape=(no_vibrations))
        array_freq = np.zeros(shape=(no_vibrations))
        array_dipole = np.zeros(shape=(no_vibrations))

        for line in lines:

            # check if string present on a current line
            if line.find('   mode   ') != -1:

                # split line in which mode is mentioned into a list containing each word seperated by a space
                mode_list = re.split(r"\s+", lines[lines.index(line)])
                value_list_freq = re.split(r"\s+", lines[lines.index(line)+2])
                value_list_mass = re.split(r"\s+", lines[lines.index(line)+14+3*no_atoms])
                value_list_dipole = re.split(r"\s+", lines[lines.index(line)+7])

                for no_column in range(no_columns_aoforce):
                    if len(mode_list) > no_column + 3:
                        array_mass[int(mode_list[no_column+2])-1] = float(value_list_mass[no_column+2])
                        if "i" not in value_list_freq[no_column+2]:
                            array_freq[int(mode_list[no_column+2])-1] = float(value_list_freq[no_column+2])
                        array_dipole[int(mode_list[no_column+2])-1] = float(value_list_dipole[no_column+2])
                        
                
                # save values of the columns into the 2d array
                for coord in range(no_atoms):

                    # split the lines containing the coords into a list containing the coords
                    value_list_x = re.split(r"\s+", lines[lines.index(line)+13+3*coord])
                    value_list_y = re.split(r"\s+", lines[lines.index(line)+14+3*coord])
                    value_list_z = re.split(r"\s+", lines[lines.index(line)+15+3*coord])
                    
                    # save x coordinates of each mode in array, if column is present
                    for no_column in range(no_columns_aoforce):
                        if len(mode_list) > no_column + 3:
                            if value_list_x[0] == "":
                                addval = 1
                            else: 
                                addval = 0
                            array_x[int(mode_list[no_column+2])-1][coord] = float(value_list_x[no_column+3+addval])
                            array_y[int(mode_list[no_column+2])-1][coord] = float(value_list_y[no_column+2])
                            array_z[int(mode_list[no_column+2])-1][coord] = float(value_list_z[no_column+2])

                
    return array_x, array_y, array_z, array_mass, array_freq, array_dipole


# method to get Lambda (how changes orbital energy with each vibration, displaced as defined in aoforce without any scaling factor)
def get_Lambda(ao_array_x, ao_array_y, ao_array_z, dE_array_x, dE_array_y, dE_array_z, no_vibrations, no_atoms):

    # dE_dX is array which contains for all vibrations the derivative
    dE_dX = np.zeros(shape=(no_vibrations))

    # displacement X (in ao_array) of each vibration is multiplied with derivative of HOMO per atom in each direction
    for vib in range(no_vibrations):
        for atom in range(no_atoms):
            dE_dX[vib] = dE_dX[vib] + float(ao_array_x[vib][atom]) * dE_array_x[atom]
            dE_dX[vib] = dE_dX[vib] + float(ao_array_y[vib][atom]) * dE_array_y[atom]
            dE_dX[vib] = dE_dX[vib] + float(ao_array_z[vib][atom]) * dE_array_z[atom]

    return dE_dX


# method to calculate the gap between HOMO and LUMO in arrays
def get_gap(HOMO, LUMO):
    gap = np.zeros(shape=(len(HOMO)))
    for vib in range(len(HOMO)):
        gap[vib] = LUMO[vib] - HOMO[vib]
    
    return gap


# writes dE_dX or Delta E to .txt file with the numbers of the vibration as defined in Turbomole (first 6 vibrations are for translation and rotation)
def write_dE_dX_to_file_tmole(filename, dE_dX_HOMO, dE_dX_HOMO_sorted, dE_dX_LUMO, dE_dX_LUMO_sorted, dE_dX_gap, dE_dX_gap_sorted, freq, val, dipole, unit):
    if unit == "eV":
        filename = filename + "_tmole_eV.txt"
        unit_factor = 27.21138
    else:
        filename = filename + "_tmole_Eh.txt"
        unit_factor = 1

    file1 = open(filename, "w")
    file1.write("LUMO\n")
    for vib in range(len(dE_dX_LUMO)):
        if dE_dX_LUMO_sorted[vib]+1 > 0:
            file1.write("vib_" + str(dE_dX_LUMO_sorted[vib]+1) + "\t" + str(val) + " " + str(round(unit_factor*dE_dX_LUMO[dE_dX_LUMO_sorted[vib]], 16)) + "\tfreq " + str(freq[dE_dX_LUMO_sorted[vib]]) + "\tdipole (a.u.) " + str(dipole[dE_dX_gap_sorted[vib]]) + "\n")

    file1.write("\nHOMO\n")
    for vib in range(len(dE_dX_HOMO)):
        if dE_dX_HOMO_sorted[vib]+1 > 0:
            file1.write("vib_" + str(dE_dX_HOMO_sorted[vib]+1) + "\t" + str(val) + " " + str(round(unit_factor*dE_dX_HOMO[dE_dX_HOMO_sorted[vib]], 16)) + "\tfreq " + str(freq[dE_dX_HOMO_sorted[vib]]) + "\tdipole (a.u.) " + str(dipole[dE_dX_gap_sorted[vib]]) + "\n")

    file1.write("\ngap\n")
    for vib in range(len(dE_dX_gap)):
        if dE_dX_gap_sorted[vib]+1 > 0:
            file1.write("vib_" + str(dE_dX_gap_sorted[vib]+1) + "\t" + str(val) + " " + str(round(unit_factor*dE_dX_gap[dE_dX_gap_sorted[vib]], 16)) + "\tfreq " + str(freq[dE_dX_gap_sorted[vib]]) + "\tdipole (a.u.) " + str(dipole[dE_dX_gap_sorted[vib]]) + "\n")

    file1.close() 


# writes dE_dX or Delta E to .txt file with the numbers of the vibration as defined in Gaussian (no vibrations for translation and rotation)
def write_dE_dX_to_file_gaussian(filename, dE_dX_HOMO, dE_dX_HOMO_sorted, dE_dX_LUMO, dE_dX_LUMO_sorted, dE_dX_gap, dE_dX_gap_sorted, freq, val, dipole, unit):
    if unit == "eV":
        filename = filename + "_gaussian_eV.txt"
        unit_factor = 27.21138
    else:
        filename = filename + "_gaussian_Eh.txt"
        unit_factor = 1

    file1 = open(filename, "w")
    file1.write("LUMO\n")
    for vib in range(len(dE_dX_LUMO)):
        if dE_dX_LUMO_sorted[vib]+1-6 > 0:
            file1.write("vib_" + str(dE_dX_LUMO_sorted[vib]+1-6) + "\t" + str(val) + " " + str(round(unit_factor*dE_dX_LUMO[dE_dX_LUMO_sorted[vib]], 16)) + "\tfreq " + str(freq[dE_dX_LUMO_sorted[vib]]) + "\tdipole (a.u.) " + str(dipole[dE_dX_gap_sorted[vib]]) + "\n")

    file1.write("\nHOMO\n")
    for vib in range(len(dE_dX_HOMO)):
        if dE_dX_HOMO_sorted[vib]+1-6 > 0:
            file1.write("vib_" + str(dE_dX_HOMO_sorted[vib]+1-6) + "\t" + str(val) + " " + str(round(unit_factor*dE_dX_HOMO[dE_dX_HOMO_sorted[vib]], 16)) + "\tfreq " + str(freq[dE_dX_HOMO_sorted[vib]]) + "\tdipole (a.u.) " + str(dipole[dE_dX_gap_sorted[vib]]) + "\n")

    file1.write("\ngap\n")
    for vib in range(len(dE_dX_gap)):
        if dE_dX_gap_sorted[vib]+1-6 > 0:
            file1.write("vib_" + str(dE_dX_gap_sorted[vib]+1-6) + "\t" + str(val) + " " + str(round(unit_factor*dE_dX_gap[dE_dX_gap_sorted[vib]], 16)) + "\tfreq " + str(freq[dE_dX_gap_sorted[vib]]) + "\tdipole (a.u.) " + str(dipole[dE_dX_gap_sorted[vib]]) + "\n")

    file1.close() 


# method to calculate dE/dT and Delta E
# dE/dT is derivative of orbital energy at temperature T
# Delta E is total energy of the orbital at temperature T
def get_dE_dsigma(filename_aoforce, filename_dEidR, name_molecule, numbering, no_atoms, no_columns_aoforce, no_HOMO, sc_gauss, no_vibrations, Temp, unit, threshold):
    # filenames to save final data
    filename_dE_dT = name_molecule + "_dE_dT"
    filename_Delta_E = name_molecule + "_Delta_E"
    print("get_dE_dsigma started")
    no_LUMO = no_HOMO + 1

    for T in range(len(Temp)):
        Delta_E_gap_tmp = []
        dE_dT_gap_tmp = []
        for tr in range(threshold):
            print(tr)
            # data needed for calculation is read from aoforce and dEidR file
            ao_array_x, ao_array_y, ao_array_z, mass, freq, dipole = aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce)
            print("aoforce read")
            dE_HOMO_x, dE_HOMO_y, dE_HOMO_z = dEidR_of_HOMO(filename_dEidR, no_atoms, no_HOMO-tr)
            dE_LUMO_x, dE_LUMO_y, dE_LUMO_z = dEidR_of_HOMO(filename_dEidR, no_atoms, no_LUMO+tr)
            print("dEidR read")
        
            # Lambda is derivative of orbital energy with coorinates of the normal modes without any temperature dependent scaling factor
            Lambda_HOMO = get_Lambda(ao_array_x, ao_array_y, ao_array_z, dE_HOMO_x, dE_HOMO_y, dE_HOMO_z, no_vibrations, no_atoms)
            Lambda_LUMO = get_Lambda(ao_array_x, ao_array_y, ao_array_z, dE_LUMO_x, dE_LUMO_y, dE_LUMO_z, no_vibrations, no_atoms)
            #print(Lambda_HOMO)
            print("Lambda calculated")
            
            # 2D arrays will containg the final data for dE/dT and Delta E for each temperature and vibration
            dE_dT_HOMO_tmp = np.zeros(shape=(threshold, no_vibrations))
            dE_dT_LUMO_tmp = np.zeros(shape=(threshold, no_vibrations))
            Delta_E_HOMO_tmp = np.zeros(shape=(threshold, no_vibrations))
            Delta_E_LUMO_tmp = np.zeros(shape=(threshold, no_vibrations))   

            # for every temperature and vibration dE/dT and Delta E will be calculated
            # Delta E is Lambda multiplied with a temperature depenend scaling factor (sigma), ref: J. Chem. Phys. 154, 244109 (2021)
            # dE/dT is derivative of Delta E per temperature
    
            for vib in range(no_vibrations):
                if freq[vib] > 0:
                    dE_dT_HOMO_tmp[tr][vib] = Lambda_HOMO[vib] * get_dsigma_dT(float(freq[vib]), float(mass[vib]), float(Temp[T]))
                    dE_dT_LUMO_tmp[tr][vib] = Lambda_LUMO[vib] * get_dsigma_dT(float(freq[vib]), float(mass[vib]), float(Temp[T]))
                    if T == 0:
                        sc_factor = sc_gauss*math.sqrt(1/(2*freq[vib]*mass[vib]))
                    else:
                        sc_factor = sc_gauss*math.sqrt(1/(2*freq[vib]*mass[vib])*mpmath.coth(1.986e-23*freq[vib]/(2*constants.k*Temp[T])))
                    Delta_E_HOMO_tmp[tr][vib] = Lambda_HOMO[vib] * sc_factor 
                    Delta_E_LUMO_tmp[tr][vib] = Lambda_LUMO[vib] * sc_factor 
            print("Delta E and dE/dT calculated")

            # gap size between HOMO and LUMO will be calculated or Delta E and dE/dT
            dE_dT_gap_tmp.append(get_gap(dE_dT_HOMO_tmp[tr], dE_dT_LUMO_tmp[tr]))
            Delta_E_gap_tmp.append(get_gap(Delta_E_HOMO_tmp[tr], Delta_E_LUMO_tmp[tr]))     
               
        # calculate average between different directions
        dE_dT_gap = np.zeros(shape=(no_vibrations))
        Delta_E_gap = np.zeros(shape=(no_vibrations))
        dE_dT_HOMO = np.zeros(shape=(no_vibrations))
        Delta_E_HOMO = np.zeros(shape=(no_vibrations))
        dE_dT_LUMO = np.zeros(shape=(no_vibrations))
        Delta_E_LUMO = np.zeros(shape=(no_vibrations))
        for tr in range(threshold):
            for vib in range(no_vibrations):
                if freq[vib] > 0:
                    dE_dT_gap[vib] = dE_dT_gap[vib] + dE_dT_gap_tmp[tr][vib] / threshold
                    Delta_E_gap[vib] = Delta_E_gap[vib] + Delta_E_gap_tmp[tr][vib] / threshold
                    dE_dT_HOMO[vib] = dE_dT_HOMO[vib] + dE_dT_HOMO_tmp[tr][vib] / threshold
                    Delta_E_HOMO[vib] = Delta_E_HOMO[vib] + Delta_E_HOMO_tmp[tr][vib] / threshold
                    dE_dT_LUMO[vib] = dE_dT_LUMO[vib] + dE_dT_LUMO_tmp[tr][vib] / threshold
                    Delta_E_LUMO[vib] = Delta_E_LUMO[vib] + Delta_E_LUMO_tmp[tr][vib] / threshold

        # sort dE/dT and Delta E according to their energies
        # sort gap by absolute values
        # vibrations are numbered as in gaussian (vib number of Tmole -6)
        dE_dT_HOMO_sorted = np.flip(np.argsort(np.absolute(dE_dT_HOMO)))
        dE_dT_LUMO_sorted = np.flip(np.argsort(np.absolute(dE_dT_LUMO)))
        Delta_E_HOMO_sorted = np.flip(np.argsort(np.absolute(Delta_E_HOMO)))
        Delta_E_LUMO_sorted = np.flip(np.argsort(np.absolute(Delta_E_LUMO)))
        dE_dT_gap_sorted = np.flip(np.argsort(np.absolute(dE_dT_gap)))
        Delta_E_gap_sorted = np.flip(np.argsort(np.absolute(Delta_E_gap)))
        filename_dE_dT_new = filename_dE_dT + "_" + str(Temp[T]) + "K"
        filename_Delta_E_new = filename_Delta_E + "_" + str(Temp[T]) + "K"

        # write E/dT and Delta E to seperate files for each temperature
        print("write results to file")
        if numbering == 'Gaussian':
            write_dE_dX_to_file_gaussian(filename_dE_dT_new, dE_dT_HOMO, dE_dT_HOMO_sorted, dE_dT_LUMO, dE_dT_LUMO_sorted, dE_dT_gap, dE_dT_gap_sorted, freq, "dE/dT", dipole, unit)
            write_dE_dX_to_file_gaussian(filename_Delta_E_new, Delta_E_HOMO, Delta_E_HOMO_sorted, Delta_E_LUMO, Delta_E_LUMO_sorted, Delta_E_gap, Delta_E_gap_sorted, freq, "Delta_E", dipole, unit)
        if numbering == 'Turbomole':
            write_dE_dX_to_file_tmole(filename_dE_dT_new, dE_dT_HOMO, dE_dT_HOMO_sorted, dE_dT_LUMO, dE_dT_LUMO_sorted, dE_dT_gap, dE_dT_gap_sorted, freq, "dE/dT", dipole, unit)
            write_dE_dX_to_file_tmole(filename_Delta_E_new, Delta_E_HOMO, Delta_E_HOMO_sorted, Delta_E_LUMO, Delta_E_LUMO_sorted, Delta_E_gap, Delta_E_gap_sorted, freq, "Delta_E", dipole, unit)



# method to calculate dsigma/dT
def get_dsigma_dT(freq, mass, Temp):
    k_eV = 8.617333262E-5
    dsigma_dT = freq*0.0001239842573148/(4*k_eV*Temp**2*math.sqrt(2*mass*freq*mpmath.coth(freq*0.0001239842573148/(2*k_eV*Temp)))*mpmath.sinh(freq*0.0001239842573148/(2*k_eV*Temp))**2)
    return dsigma_dT



# method to calculate how much orbital enregy of HOMO, LUMO and gap changes with each vibration
# method sorts those energy and prints it into a .txt file
def read_files_get_dE_dX(filename_aoforce, no_atoms, no_columns_aoforce, filename_dEidR, no_HOMO, no_vibrations, name_molecule, numbering, unit):

    no_LUMO = no_HOMO+1
    
    # generate arrays containing displacement vectors of aoforce
    ao_array_x, ao_array_y, ao_array_z, mass, array_freq, dipole = aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce)

    # get derivative of orbital energy per vibration of HOMO, LUMO and the gap (dE/dX), which is Lambda
    dE_HOMO_x, dE_HOMO_y, dE_HOMO_z = dEidR_of_HOMO(filename_dEidR, no_atoms, no_HOMO)
    dE_LUMO_x, dE_LUMO_y, dE_LUMO_z = dEidR_of_HOMO(filename_dEidR, no_atoms, no_LUMO)
    dE_dX_HOMO = get_Lambda(ao_array_x, ao_array_y, ao_array_z, dE_HOMO_x, dE_HOMO_y, dE_HOMO_z, no_vibrations, no_atoms)
    dE_dX_LUMO = get_Lambda(ao_array_x, ao_array_y, ao_array_z, dE_LUMO_x, dE_LUMO_y, dE_LUMO_z, no_vibrations, no_atoms)
    dE_dX_gap = get_gap(dE_dX_HOMO, dE_dX_LUMO)

    # sort dE/dX energies and write them into a file
    # vibrations are either numbered as in Tmole or in gaussian (vib number of Tmole -6)
    dE_dX_HOMO_sorted = np.argsort(dE_dX_HOMO)
    dE_dX_LUMO_sorted = np.argsort(dE_dX_LUMO)
    dE_dX_gap_sorted = np.argsort(dE_dX_gap)
    filename = name_molecule + "_dE_dX"
    if numbering == 'Gaussian':
        write_dE_dX_to_file_gaussian(filename, dE_dX_HOMO, dE_dX_HOMO_sorted, dE_dX_LUMO, dE_dX_LUMO_sorted, dE_dX_gap, dE_dX_gap_sorted, array_freq, "dE/dX", dipole, unit)
    if numbering == 'Turbomole':
        write_dE_dX_to_file_tmole(filename, dE_dX_HOMO, dE_dX_HOMO_sorted, dE_dX_LUMO, dE_dX_LUMO_sorted, dE_dX_gap, dE_dX_gap_sorted, array_freq, "dE/dX", dipole, unit)


# method to get Delta_E or dE/dT for a vibration dependent on several temperatures
def vib_per_temp(Temp, name_molecule, no_vib, value, unit):

    # arrays which will contain the values for dE/dT or Delta E
    homo = np.zeros(shape=(len(Temp)))
    lumo = np.zeros(shape=(len(Temp)))
    gap = np.zeros(shape=(len(Temp)))

    # get values for every temperature and save them in arrays
    for T in range(len(Temp)):

        # filename in which the values which are evaluated saved and read it
        if value == "dE_dT":
            name_file = name_molecule + "_dE_dT_" + str(Temp[T]) + "K_gaussian_" + str(unit) + ".txt"
        else:
            name_file = name_molecule + "_Delta_E_" + str(Temp[T]) + "K_gaussian_" + str(unit) + ".txt"

        with open(name_file, 'r') as fp:
            lines = fp.readlines()
            
            orb = 1
            for line in lines:

                # check if requested vibration is present in the current line
                tofind = "vib_"  + str(no_vib) + "\t"                
                if line.find(tofind) != -1:
                    value_list = re.split(r"\s+", line)
                    if orb == 1:
                        lumo[T] = value_list[2]
                    if orb == 2:
                        homo[T] = value_list[2]
                    if orb == 3:
                        gap[T] = value_list[2]
                    orb = orb + 1
    return homo, lumo, gap


# method to save Delta E or dE/dT at different temperatures in a seperate file for each vibration
def save_vib_per_temp(name_molecule, no_vib, homo, lumo, gap, Temp, value, unit):

    # filename in which data will be saved
    if value == "dE_dT":
        filename = name_molecule + "_dE_dT_vib_" + str(no_vib) + "_" + str(unit) + ".txt"
    else:
        filename = name_molecule + "_Delta_E_vib_" + str(no_vib) + "_" + str(unit) + ".txt"
    file1 = open(filename, "w")

    for T in range(len(homo)):
        file1.write(str(Temp[T]) + " K\tHOMO " + str(round(homo[T], 16)) + "\tlumo " + str(round(lumo[T], 16)) + "\tgap " + str(round(gap[T], 16)) + "\n")

    file1.close() 
    

# method to calculate Delta E or dE/T for several vibrations in no_vib and at temperatures in Temp
def vibrations_per_temperatures(no_vib, Temp, name_molecule, val, unit):
    for vibr in no_vib:
        for value in val:
            homo, lumo, gap = vib_per_temp(Temp, name_molecule, vibr, value, unit)
            save_vib_per_temp(name_molecule, vibr, homo, lumo, gap, Temp, value, unit)


# method to get most dominating frequencies at each temperatures, either in dE/dT or in Delta E
def dominating_freq(Temp, threshold, no_vibrations, val, unit, name_molecule):

    # get data either for Delta E or dE/dT
    for value in val:
        print(value)

        homo, lumo, gap, freq = read_files(Temp, no_vibrations, val, unit, name_molecule)
        homo_sorted = np.zeros(shape=(len(Temp), no_vibrations-6))
        lumo_sorted = np.zeros(shape=(len(Temp), no_vibrations-6))
        gap_sorted = np.zeros(shape=(len(Temp), no_vibrations-6))

        for T in range(len(Temp)):
            homo_sorted[T] = np.flip(np.argsort(np.absolute(homo[T])))
            lumo_sorted[T] = np.flip(np.argsort(np.absolute(lumo[T])))
            gap_sorted [T]= np.flip(np.argsort(np.absolute(gap[T])))

        # write the results into a txt file
        filename = name_molecule + "_" + value + "_dom_freq.txt"
        file1 = open(filename, "w")
        file1.write("LUMO\n")
        for T in range(len(Temp)):
            lumo_line = "Temp = " + str(Temp[T])
            for thresh in range(threshold):
                lumo_line = lumo_line + "\tvib_" + str(int(lumo_sorted[T][thresh]+1)) + "\t" + str(freq[int(lumo_sorted[T][thresh])]) + "\t" + str(value) + "\t" + str(lumo[T][int(lumo_sorted[T][thresh])])
            lumo_line = lumo_line + "\n"
            file1.write(lumo_line)
 

        file1.write("\nHOMO\n")
        for T in range(len(Temp)):
            homo_line = "Temp = " + str(Temp[T])
            for thresh in range(threshold):
                homo_line = homo_line + "\tvib_" + str(int(homo_sorted[T][thresh]+1)) + "\t" + str(freq[int(homo_sorted[T][thresh])]) + "\t" + str(value) + "\t" + str(homo[T][int(homo_sorted[T][thresh])])
            homo_line = homo_line + "\n"
            file1.write(homo_line)

        file1.write("\ngap\n")
        for T in range(len(Temp)):
            gap_line = "Temp = " + str(Temp[T])
            for thresh in range(threshold):
                gap_line = gap_line + "\tvib_" + str(int(gap_sorted[T][thresh]+1)) + "\t" + str(freq[int(gap_sorted[T][thresh])]) + "\t" + str(value) + "\t" + str(gap[T][int(gap_sorted[T][thresh])])
            gap_line = gap_line + "\n"
            file1.write(gap_line)

        file1.close() 


# method to read in .txt files of certain temperatures
def read_files(Temp, no_vibrations, val, unit, name_molecule):

    homo = np.zeros(shape=(len(Temp), no_vibrations-6))
    lumo = np.zeros(shape=(len(Temp), no_vibrations-6))
    gap = np.zeros(shape=(len(Temp), no_vibrations-6))
    freq = np.zeros(shape=(no_vibrations-6))


    for T in range(len(Temp)):

        if val == "dE_dT":
            name_file = name_molecule + "_dE_dT_" + str(Temp[T]) + "K_gaussian_" + str(unit) + ".txt"
        else:
            name_file = name_molecule + "_Delta_E_" + str(Temp[T]) + "K_gaussian_" + str(unit) +".txt"

        with open(name_file, 'r') as fp:
            lines = fp.readlines()

            for line in lines:
                value_list = re.split(r"\s+", lines[lines.index(line)])
                if lines.index(line) < 1+no_vibrations-6 and lines.index(line) > 0:
                    vib = re.split(r"_", value_list[0])
                    freq[int(vib[1])-1] = value_list[4]
                    lumo[T][int(vib[1])-1] = value_list[2]
                if lines.index(line) < 3+2*(no_vibrations-6) and lines.index(line) > no_vibrations-6+2:
                    vib = re.split(r"_", value_list[0])
                    homo[T][int(vib[1])-1] = value_list[2]
                if lines.index(line) < 5+3*(no_vibrations-6) and lines.index(line) > 2*(no_vibrations-6)+4:
                    vib = re.split(r"_", value_list[0])
                    gap[T][int(vib[1])-1] = value_list[2]

    return homo, lumo, gap, freq


def get_raman(filename_control, no_vibrations):
    # open file and read all lines in a list
    with open(filename_control, 'r') as fp:
        lines = fp.readlines()

        pol = []
        freq = []
        for line in range(len(lines)):
            # check if 'mode' present on a current line
            if lines[line].find('$raman') != -1:
                for i in range(no_vibrations-6):
                    value_list = re.split(r"\s+", lines[line+12+i])
                    print(value_list)
                    if "0.0" not in value_list[2]:
                        pol.append(0.5*(float(value_list[7].replace("D", "E"))+float(value_list[8].replace("D", "E"))))
                        freq.append(value_list[3])
    return freq, pol

