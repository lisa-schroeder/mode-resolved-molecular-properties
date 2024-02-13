from generate_vibs import generate_coords_tmole, generate_coords_energy, coord_to_xyz
from evib import read_files_get_dE_dX, get_dE_dsigma, vibrations_per_temperatures, dominating_freq, plot_IR, plot_IR_broadened, plot_evib_broadened, plot_Raman_broadened
from plot import plot_IR, plot_IR_broadened, plot_evib_broadened, plot_Raman_broadened, plot_dsigma_dT
import numpy as np

# define necessary variables
name_molecule = "l-p2+"
used_program = "Turbomole"                               # used program can be gaussian or turbomole, numbering of the frequencies changes
#used_program = "Gaussian"
#numbering = 'Turbomole'                                 # numbering of frequencies as in gaussian 
numbering = 'Gaussian'                                   # numbering of frequencies as in turbomole (difference by 6, tmole also counts translation and rotation)
no_atoms = 80                                            # number of atoms the molecule has
no_vibrations = 240                                      # number of vibrations the molecule has in total
#no_vibration = np.arange(7, no_atoms*3+1)               # number of vibration for which properties will be calculated
no_vibration = [35]
filename_aoforce = 'aoforce.out'                         # names of files containing the aoforce and coord data
filename_coord = 'coord'

# generate structures
Temp = [0.01]                                            # temperature at which the displacement will be calculated
scale_dir = [-1, 1]                                      # direction in which atoms will be displaced (1 and -1)

# only for scaling with certain energy
scale_energy = [0]                                       # amount of energy in eV which is pumped into the vibration

# evib
val = ["Delta_E"]
#val = ["Delta_E", "dE_dT"]                              # whether dE/dT or Delta E is calculated for the vibrations at the temperatures
threshold = 2                                            # number of most important vibrations in method dominating_freq
no_HOMO = 211                                            # number of orbital which is HOMO, which will be evaluated
unit = "eV"                                              # unit in eV or Eh (for Delta E and dE/dT, only fits to some orbitals!)
filename_dEidR = 'dEidR_b.dat'                           # names of file containing evib data
filename_control = 'control'                             # name of control file containing raman data
broadening = 8                                           # broadening, neccessary for calculating spectra
res = 1000                                               # resolution of calculated spectrum
orb = "HOMO"                                             # orbital which is used for a evib spectrum



# constants that usually don't need to be changed
no_columns_aoforce = 6                                   # amount of columns in the Turbomole aoforce file
sc_tmole = 5.15003E-07                                   # internal scaling factor of Turbomole
sc_gauss = 10.97213052361718                             # internal scaling for energy per temperature

no_vibration_tmole = []                                  # converting numbering of vibrations either to Tmole or Gaussian
no_vibration_gauss = []
if used_program == "Gaussian":
    for i in range(len(no_vibration)):
        no_vibration_tmole.append(no_vibration[i] + 6)
    no_vibrations_tmole = no_vibrations + 6
    no_vibration_gauss = no_vibration
    no_vibrations_gauss = no_vibrations
if used_program == "Turbomole":
    no_vibration_tmole = no_vibration
    no_vibrations_tmole = no_vibrations
    for i in range(len(no_vibration)):
        no_vibration_gauss.append(no_vibration[i] - 6)
    no_vibrations_gauss = no_vibrations -6




#generate_coords_tmole(Temp, no_atoms, filename_aoforce, filename_coord, no_columns_aoforce, no_vibration_tmole, scale_dir, sc_tmole, sc_gauss, numbering)
#generate_coords_energy(no_atoms, filename_aoforce, filename_coord, no_columns_aoforce, no_vibration_tmole, scale_dir, sc_tmole, scale_energy, numbering)
#coord_to_xyz(filename_coord, no_atoms)
#coord_to_xyz_flipped(filename_coord, no_atoms)

# read_files_get_dE_dX(filename_aoforce, no_atoms, no_columns_aoforce, filename_dEidR, no_HOMO, no_vibrations_tmole, name_molecule, numbering, unit, dir)
get_dE_dsigma(filename_aoforce, filename_dEidR, name_molecule, numbering, no_atoms, no_columns_aoforce, no_HOMO, sc_gauss, no_vibrations_tmole, Temp, unit, threshold)
#vibrations_per_temperatures(no_vibration_tmole, Temp, name_molecule, val, unit)
#dominating_freq(Temp, threshold, no_vibrations_tmole, val, unit, name_molecule)
#plot_IR(filename_aoforce, no_atoms, no_columns_aoforce, name_molecule)
#plot_IR_broadened(filename_aoforce, no_atoms, no_columns_aoforce, name_molecule, broadening, res)
plot_evib_broadened(Temp, no_vibrations, val, unit, name_molecule, broadening, res, orb)
#plot_Raman_broadened(filename_control, name_molecule, broadening, res, no_vibrations)

