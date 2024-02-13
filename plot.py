
import matplotlib.pyplot as plt
import numpy as np
import math

# plot IR spectrum as from an aoforce file
def plot_IR(filename_aoforce, no_atoms, no_columns_aoforce, name_molecule):

    array_x, array_y, array_z, array_mass, array_freq, array_dipole = aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce)

    fig, ax = plt.subplots()
    ax.stem(array_freq, array_dipole, markerfmt=' ')
    plt.title(name_molecule)
    plt.xlabel("frequency in cm-1")
    plt.ylabel("|dDIP/dQ| in a.u.")

    plt.show()

# plot evib spectrum 
def plot_evib_broadened(Temp, no_vibrations, val, unit, name_molecule, broadening, res, orb):

    homo, lumo, gap, array_freq = read_files(Temp, no_vibrations, val[0], unit, name_molecule)
    if orb == "HOMO":
        data = homo[0]
    elif orb == "LUMO":
        data = lumo[0]
    elif orb == "GAP":
        data = gap[0]
    x_axis = np.arange(0, (np.max(array_freq)+100), np.max(array_freq)/res)
    spectrum = []
    for x in x_axis:
        x_value = 0.0
        for freq in range(len(array_freq)):
            x_value = x_value + 1/(2*np.pi*broadening)*np.exp(-(x-float(array_freq[freq]))**2/broadening**2)*float(np.abs(data[freq]))
            #spectrum[x] = spectrum[x] + 1/(np.pi*broadening)*np.exp((x_axis[x]-array_freq[freq])/broadening**2)*array_dipole[freq]
        spectrum.append(x_value)
        #print(x_value)
    #print(spectrum)

    filename_evib = name_molecule + "_evib_" + orb + "_br" + str(broadening)
    file1 = open(filename_evib, "w")
    for line in range(len(x_axis)):
        file1.write(str(x_axis[line]) + "\t" + str(spectrum[line]) + " \n")
    file1.close() 

    fig, ax = plt.subplots()
    plt.plot(x_axis, spectrum)
    plt.title(name_molecule)
    plt.xlabel("frequency in cm-1")
    plt.ylabel("Delta E in eV")

    plt.show()

def plot_IR_broadened(filename_aoforce, no_atoms, no_columns_aoforce, name_molecule, broadening, res):

    array_x, array_y, array_z, array_mass, array_freq, array_dipole = aoforce_in_arrays(filename_aoforce, no_atoms, no_columns_aoforce)
    x_axis = np.arange(0, 3500, 3500/res)
    print(x_axis)
    spectrum = []
    for x in x_axis:
        x_value = 0.0
        for freq in range(len(array_freq)):
            x_value = x_value + 1/(2*np.pi*broadening)*np.exp(-(x-float(array_freq[freq]))**2/broadening**2)*float(array_dipole[freq])
            #spectrum[x] = spectrum[x] + 1/(np.pi*broadening)*np.exp((x_axis[x]-array_freq[freq])/broadening**2)*array_dipole[freq]
        spectrum.append(x_value)
        #print(x_value)
    print(spectrum)

    filename_IR = name_molecule + "_IR_br" + str(broadening)
    file1 = open(filename_IR, "w")
    for line in range(len(x_axis)):
        file1.write(str(x_axis[line]) + "\t" + str(spectrum[line]) + " \n")
    file1.close() 

    fig, ax = plt.subplots()
    plt.plot(x_axis, spectrum)
    plt.title(name_molecule)
    plt.xlabel("frequency in cm-1")
    plt.ylabel("|dDIP/dQ| in a.u.")

    plt.show()



def plot_Raman_broadened(filename_control, name_molecule, broadening, res, no_vibrations):
    
    array_freq, array_raman = get_raman(filename_control, no_vibrations)
    print(array_freq)

    x_axis = np.arange(0, 3500, 3500/res)

    spectrum = []
    for x in x_axis:
        x_value = 0.0
        for freq in range(len(array_freq)):
            x_value = x_value + 1/(2*np.pi*broadening)*np.exp(-(x-float(array_freq[freq]))**2/broadening**2)*float(array_raman[freq])
            #spectrum[x] = spectrum[x] + 1/(np.pi*broadening)*np.exp((x_axis[x]-array_freq[freq])/broadening**2)*array_dipole[freq]
        spectrum.append(x_value)
        #print(x_value)


    filename_IR = name_molecule + "_Raman_br" + str(broadening)
    file1 = open(filename_IR, "w")
    for line in range(len(x_axis)):
        file1.write(str(x_axis[line]) + "\t" + str(spectrum[line]) + " \n")
    file1.close() 

    fig, ax = plt.subplots()
    plt.plot(x_axis, spectrum)
    plt.title(name_molecule)
    plt.xlabel("frequency in cm-1")
    plt.ylabel("|dDIP/dQ| in a.u.")

    plt.show()



# method to plot dsigma/dT at different temperatures, dependent on a mass and frequency
def plot_dsigma_dT(freq, mass, Temp):

    dsigma_dT = np.zeros(shape=(len(Temp)))

    for T in range(len(Temp)):
        dsigma_dT[T] = get_dsigma_dT(freq, mass, Temp[T])
    
    print("dsigma_dT = " + str(dsigma_dT))
    plt.plot(Temp, dsigma_dT)
    plt.show()
