import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from itertools import cycle

#Increase font size
plt.rcParams.update({'font.size': 14})
cwd = os.getcwd()

def get_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Plot various quantities using the ledger file''')

    p.add_argument('-plot_G_values', dest='plot_G_values', action='store_true',
                  help='Plot the G_values of grouped (IBK grouping) species in units of [radicals/(100 eV)]')

    p.add_argument('-path', dest='paths', action="append", type=str,
        help='The absolute path to the silo files, can be a list for multiple databases')
    return p.parse_args()

def idx(cIx):
    "Python starts at 0, and I am not dealing with it any longer..."
    return cIx - 1 # Converts cIx to Python index

class Ledger:
    def __init__(self, file, name=""):
        # Lowest x coordinate (start of the plot range)
        self.xmin = 0

        # Guess the default x-axis based on which files are supplied
        if os.path.isfile(file+"streamer_length.txt"):
            self.xaxis = 'length'
        elif os.path.isfile(file+"streamer_Lgap.txt"):
            self.xaxis = 'Lgap'
        else:
            self.xaxis = 'time'

        # Try to load any data that is present
        print_string = ""
        try: self.L = np.loadtxt(file+"streamer_length.txt", dtype=np.dtype('float'))
        except: print_string += (" length,")
        try: self.R = np.loadtxt(file+"streamer_radius.txt", dtype=np.dtype('float'))
        except: print_string += (" radius,")
        try: self.V = np.loadtxt(file+"streamer_volume.txt", dtype=np.dtype('float'))
        except: print_string += (" volume,")
        try: self.Lgap = np.loadtxt(file+"streamer_Lgap.txt", dtype=np.dtype('float'))
        except: print_string += (" Lgap,")
        try: self.energy = np.loadtxt(file+"streamer_energy.txt", dtype=np.dtype('float'))
        except: print_string += (" energy,")
        try: self.Emax = np.loadtxt(file+"streamer_Emax.txt", dtype=np.dtype('float'))
        except: print_string += (" Emax,")
        finally:
            if (not print_string == ""):
                print("Was not able to load:"+print_string+" to data structure "+name)

        # Load the ledger
        file = file + "sim_cs_ledger.txt"

        # Header data is loaded seperately
        self.cIx = np.loadtxt(file, skiprows = 4, max_rows = 1, dtype=np.dtype('i4'))
        self.gas = np.loadtxt(file, skiprows = 5, max_rows = 1, dtype=np.dtype('U'))
        self.coll_type = np.loadtxt(file, skiprows = 6, max_rows = 1, dtype=np.dtype('U'))

        # now load collision data to ledger
        self.ledger = np.loadtxt(file, skiprows = 7, dtype=np.dtype('float'))
        # Extract timestamps
        self.timesteps = self.ledger[:, -1]
        # Remove timestamps from ledger (and first entry which are the null collisions)
        self.ledger = self.ledger[:, 1:-1]
        # Assuming equal timespacing calculate dt
        self.dt = self.timesteps[1] - self.timesteps[0]

        # Finally, always set generic labels (ion/exc/...) to the set
        self.generic_grouping()

        self.name = name

    def generic_grouping(self):
        self.i_N2_ion  = []
        self.i_O2_ion  = []
        self.i_CH4_ion = []

        self.i_N2_exc  = []
        self.i_O2_exc  = []
        self.i_CH4_exc = []

        self.i_N2_attach  = []
        self.i_O2_attach  = []
        self.i_CH4_attach = []

        # Now set the indices according to collision type
        for n in  range(len(self.cIx)):
            if (self.gas[n]=='N2' and self.coll_type[n] =='Ionization'):
                self.i_N2_ion.append(idx(self.cIx[n]))
            elif (self.gas[n]=='O2' and self.coll_type[n] =='Ionization'):
                self.i_O2_ion.append(idx(self.cIx[n]))
            elif (self.gas[n]=='CH4' and self.coll_type[n] =='Ionization'):
                self.i_CH4_ion.append(idx(self.cIx[n]))

            elif (self.gas[n]=='N2' and self.coll_type[n] =='Excitation'):
                self.i_N2_exc.append(idx(self.cIx[n]))
            elif (self.gas[n]=='O2' and self.coll_type[n] =='Excitation'):
                self.i_O2_exc.append(idx(self.cIx[n]))
            elif (self.gas[n]=='CH4' and self.coll_type[n] =='Excitation'):
                self.i_CH4_exc.append(idx(self.cIx[n]))

            elif (self.gas[n]=='N2' and self.coll_type[n] =='Attachment'):
                self.i_N2_attach.append(idx(self.cIx[n]))
            elif (self.gas[n]=='O2' and self.coll_type[n] =='Attachment'):
                self.i_O2_attach.append(idx(self.cIx[n]))
            elif (self.gas[n]=='CH4' and self.coll_type[n] =='Attachment'):
                self.i_CH4_attach.append(idx(self.cIx[n]))
            elif (self.coll_type[n] != 'Elastic'): # Ignore elastic collisions
                print('Collision not recognized')
                print(self.cIx[n])
                print(self.gas[n])
                print(self.coll_type[n])
                print(self.gas[n]=='O2' and self.coll_type[n] =='Atttachment')
                1/0 # I am too lazy to look for the proper clausule

    def set_grouping_all_repo(self):
        # Assuming that by-products are in the lowest energy state
        self.i_O_rad = [17, 28, 30, 31]
        self.i_H2    = [34, 35, 40, 41, 42, 45, 46, 48]
        self.i_H_rad = [33, 36, 40, 43, 47, 48]
        self.i_N_rad = [13, 15, 16]

    def set_grouping_IBK(self):
        # O2:  Itikawa
        # CH4: Song-Bouwman
        # N2:  Kwaguchi
        self.i_O_rad = [36, 47, 49, 50]
        self.i_H2    = [52, 53, 58, 59, 60, 63, 64, 66]
        self.i_H_rad = [51, 54, 58, 61, 65, 66]
        self.i_N_rad = [29, 33, 34]

        self.i_N2_triplets = [0, 1, 2, 3, 4, 5, 6, 7]
        self.i_O2_singlets = [41, 42]

        self.i_N2_vib = [19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
        self.i_O2_vib = [45, 46, 47]
        self.i_CH4_vib = [57, 58]

        AE = self.i_hate_this_IBK()


        #now for activation energies
        self.i_O_rad_AE = np.array([AE[36], AE[47], AE[49], AE[50]])
        self.i_H2_AE    = np.array([AE[52], AE[53], AE[58], AE[59], AE[60], AE[63], AE[64], AE[66]])
        self.i_H_rad_AE = np.array([AE[51], AE[54], AE[58], AE[61], AE[65], AE[66]])
        self.i_N_rad_AE = np.array([AE[29], AE[33], AE[34]])

        self.i_N2_triplets_AE = np.array([AE[0], AE[1], AE[2], AE[3], AE[4], AE[5], AE[6], AE[7]])
        self.i_O2_singlets_AE = np.array([AE[41], AE[42]])

        self.i_N2_vib_AE = [AE[19], AE[20], AE[21], AE[22], AE[23], AE[24], AE[25], AE[26], AE[27], AE[28]]
        self.i_O2_vib_AE = [AE[45], AE[46], AE[47]]
        self.i_CH4_vib_AE = [AE[57], AE[58]]

    def plot_all_repo(self, factor=1, **kwargs):
        # plt.figure()
        self.plot_single_species_weighted(self.i_O_rad, np.array([2, 2, 2, 2])*factor, label="O", color='r', **kwargs)
        self.plot_single_species_weighted(self.i_H2, np.array([1, 1, 1, 2, 1, 1, 2, 1])*factor, label="H$_2$", color='b', **kwargs)
        self.plot_single_species_weighted(self.i_H_rad, np.array([1, 1, 1, 1, 1, 1])*factor, label="H", color='g', **kwargs)
        self.plot_single_species_weighted(self.i_N_rad, np.array([2, 2, 2])*factor, label="N", color='m', **kwargs)

        # plt.legend()
        plt.xlabel(self.xaxis)
        plt.ylabel("particle number")

    def plot_IBK(self, factor=1.0, **kwargs):
        self.plot_single_species_weighted(self.i_N2_triplets, np.array([1, 1, 1, 1, 1, 1, 1, 1])*factor, label="N$_2$ triplets", color='k', **kwargs)
        self.plot_single_species_weighted(self.i_N_rad, np.array([2, 2, 2])*factor, label="N", color='m', **kwargs)
        self.plot_single_species_weighted(self.i_O_rad, np.array([2, 2, 2, 2])*factor, label="O", color='r', **kwargs)
        self.plot_single_species_weighted(self.i_H_rad, np.array([1, 1, 1, 1, 1, 1])*factor, label="H", color='g', **kwargs)
        self.plot_single_species_weighted(self.i_O2_singlets, np.array([1, 1])*factor, label="O$_2$ singlets", color='tab:orange', **kwargs)
        self.plot_single_species_weighted(self.i_H2, np.array([1, 1, 1, 2, 1, 1, 2, 1])*factor, label="H$_2$", color='b', **kwargs)
        # self.plot_single_species_weighted(self.i_N2_vib, np.ones_like(self.i_N2_vib)*factor, label="N$_2$ vib", color='tab:red', **kwargs)
        # self.plot_single_species_weighted(self.i_O2_vib, np.ones_like(self.i_O2_vib)*factor, label="O$_2$ vib", color='tab:blue', **kwargs)
        # self.plot_single_species_weighted(self.i_CH4_vib, np.ones_like(self.i_CH4_vib)*factor, label="CH$_4$ vib", color='tab:brown', **kwargs)
        # plt.legend()
        plt.xlabel(self.xaxis)
        plt.ylabel("particle number")

    def plot_IBK_AE(self, factor=1.0, **kwargs):
        # print(np.array([1, 1, 1, 1, 1, 1, 1, 1])*factor*self.i_N2_triplets_AE)
        self.plot_single_species_weighted(self.i_N2_triplets, np.array([1, 1, 1, 1, 1, 1, 1, 1])*factor*self.i_N2_triplets_AE, label="N$_2$ triplets", color='k', **kwargs)
        self.plot_single_species_weighted(self.i_N_rad, np.array([2, 2, 2])*factor*self.i_N_rad_AE, label="N", color='m', **kwargs)
        self.plot_single_species_weighted(self.i_O_rad, np.array([2, 2, 2, 2])*factor*self.i_O_rad_AE, label="O", color='r', **kwargs)
        self.plot_single_species_weighted(self.i_H_rad, np.array([1, 1, 1, 1, 1, 1])*factor*self.i_H_rad_AE, label="H", color='g', **kwargs)
        self.plot_single_species_weighted(self.i_O2_singlets, np.array([1, 1])*factor*self.i_O2_singlets_AE, label="O$_2$ singlets", color='tab:orange', **kwargs)
        self.plot_single_species_weighted(self.i_H2, np.array([1, 1, 1, 2, 1, 1, 2, 1])*factor*self.i_H2_AE, label="H$_2$", color='b', **kwargs)
        # self.plot_single_species_weighted(self.i_N2_vib, np.ones_like(self.i_N2_vib)*factor*self.i_N2_vib_AE, label="N$_2$ vib", color='tab:red', **kwargs)
        # self.plot_single_species_weighted(self.i_O2_vib, np.ones_like(self.i_O2_vib)*factor*self.i_O2_vib_AE, label="O$_2$ vib", color='tab:blue', **kwargs)
        # self.plot_single_species_weighted(self.i_CH4_vib, np.ones_like(self.i_CH4_vib)*factor*self.i_CH4_vib_AE, label="CH$_4$ vib", color='tab:brown', **kwargs)
        # plt.legend()
        plt.xlabel(self.xaxis)
        plt.ylabel("particle number")

    def plot_G_values(self, **kwargs):
        J_to_eV = 1.0/6.242e+18
        self.plot_IBK(factor=100*J_to_eV, averaged="energy", scale='linear', **kwargs)
        plt.ylabel("G (1/(100 eV))")
        # plt.legend()

    def print_all_G_values(self, **kwargs):
        J_to_eV = 1.0/6.242e+18
        coeff = 100 * J_to_eV / self.energy[-1]
        for n in  range(len(self.cIx)):
            print(str('{:0.3e}'.format(self.ledger[-1, n]*coeff)))

        print("with activation energies")
        thresholds = self.i_hate_this_IBK()
        for n in range(len(thresholds)):
            print('{:0.3f}'.format(thresholds[n]))

    def plot_E_sink(self, **kwargs):

        self.plot_IBK_AE(factor=1.0, averaged=None, scale='linear', **kwargs)
        plt.ylabel("energy spent [eV]")

    def plot_generic_grouped_species(self, print_CH4=True, **kwargs):
        "Plot the grouped CAS, ions/excited states/fragmented fuel/dissociated oxygen"
        self.plot_single_species(self.i_N2_ion, label="N$_2$ ions", linestyle='-', color='tab:blue', **kwargs)
        self.plot_single_species(self.i_O2_ion, label="O$_2$ ions", linestyle='-', color='tab:green', **kwargs)
        if (print_CH4):
            self.plot_single_species(self.i_CH4_ion, label="CH$_4$ ions", linestyle='-', color='tab:brown', **kwargs)

        self.plot_single_species(self.i_N2_exc, label="N$_2$ exc.", linestyle=':', color='tab:blue', **kwargs)
        self.plot_single_species(self.i_O2_exc, label="O$_2$ exc.", linestyle=':', color='tab:green', **kwargs)
        if (print_CH4):
            self.plot_single_species(self.i_CH4_exc, label="CH$_4$ exc.", linestyle=':', color='tab:brown', **kwargs)

        self.plot_single_species(self.i_N2_attach, linestyle='-.', label="N2_attach")
        self.plot_single_species(self.i_O2_attach, linestyle='-.', label="O2_attach")
        self.plot_single_species(self.i_CH4_attach, linestyle='-.', label="CH4_attach")

        plt.legend()
        # plt.xlabel("time")
        plt.ylabel("particle number")
        plt.title("Grouped species")

    def plot_single_species(self, index, **kwargs):
        if self.xaxis == 'length':
            xsteps = self.L * 1e3 # set x-axis to millimeters
        elif self.xaxis == 'Lgap':
            xsteps = self.Lgap * 1e3 # set x-axis to millimeters
        elif self.xaxis == 'time':
            xsteps = self.timesteps * 1e9 # set x-axis to nanoseconds
        elif self.xaxis == 'volume':
            xsteps = self.V * 1e9 #set x-axis to cubic millimeters
        else:
            print("Given xaxis not permitted. Use length, time or volume")
            exit()

        plt.semilogy(xsteps, np.sum(self.ledger[:, index], 1), **kwargs)

    # def plot_single_species_weighted(self, index, label, weight):
    def plot_single_species_weighted(self, index, weight, scale='semilog', averaged=None, **kwargs):
        if self.xaxis == 'length':
            xsteps = self.L * 1e3 # set x-axis to millimeters
        elif self.xaxis == 'Lgap':
            xsteps = self.Lgap * 1e3 # set x-axis to millimeters
        elif self.xaxis == 'time':
            xsteps = self.timesteps * 1e9 # set x-axis to nanoseconds
        elif self.xaxis == 'volume':
            xsteps = self.V * 1e9 #set x-axis to cubic millimeters
        else:
            print("Given xaxis not permitted. Use length, time or volume")
            exit()

        if (averaged == "volume"):
            coeff = self.V # Keep in units of m^(-3)
        elif (averaged == "energy"):
            coeff = self.energy
        elif (averaged is None):
            coeff = np.ones_like(xsteps)
        else:
            print("ERROR: averaging procedure not recognized")
            sys.exit()

        if scale == 'semilog':
            plt.semilogy(xsteps[xsteps>=self.xmin], np.sum(np.tile(weight, (np.shape(self.ledger)[0], 1)) * self.ledger[:, index], 1)[xsteps>=self.xmin]/coeff[xsteps>=self.xmin], **kwargs)
        if scale == 'linear':
            plt.plot(xsteps[xsteps>=self.xmin], np.sum(np.tile(weight, (np.shape(self.ledger)[0], 1)) * self.ledger[:, index], 1)[xsteps>=self.xmin]/coeff[xsteps>=self.xmin], **kwargs)

    def plot_excitations(self, index, label):
        exc_total = np.tile(np.sum(self.ledger[:, index], 1, keepdims=True), (1, np.size(index)))
        exc_relative = np.divide(self.ledger[:, index], exc_total + 1e-10)

        fig, ax = plt.subplots(1)
        ax.stackplot(self.timesteps[1:], np.transpose(exc_relative[1:, :]), labels=np.arange(1, np.size(index)+1).astype(str))
        # Switch order of labels and handles
        handles,labels = ax.get_legend_handles_labels()
        handles = np.flip(handles)
        labels = np.flip(labels)
        # Set x and y margins to zero for better fit
        ax.margins(x = 0)
        ax.margins(y = 0)
        # Move legend left of the plot
        ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xlabel("time")
        plt.ylabel("fraction")
        plt.title("Relative population of excited states of " + label)

    def print_volume_weighted_species(self, timestep):
        # Only plot the final timestep
        for index in range(len(self.cIx)):
            print("cIx = " + str(index+1) + ":  " + str('{:0.2e}'.format(self.ledger[timestep, index]/self.V[timestep])))

    def print_weighted_species(self, timestep):
        # Only plot the final timestep
        for index in range(len(self.cIx)):
            print("cIx = " + str(index+1) + ":  " + str('{:0.2e}'.format(self.ledger[timestep, index])))

    def i_hate_this_IBK(self):
        #This should be in the coll ledger. But it is maybe not worth the effort (anymore)
        activation_energies = np.array([0.988383E-18,
        0.130818e-17,
        0.117808E-17,
        0.176752E-17,
        0.190194E-17,
        0.208041E-17,
        0.205234E-17,
        0.117952E-17,
        0.196347E-17,
        0.134551E-17,
        0.136970E-17,
        0.205936E-17,
        0.200267E-17,
        0.206865E-17,
        0.207232E-17,
        0.209933E-17,
        0.240326E-21,
        0.461427E-19,
        0.437074E-18,
        0.918047E-19,
        0.136986E-18,
        0.181527E-18,
        0.225586E-18,
        0.269005E-18,
        0.311944E-18,
        0.354241E-18,
        0.395898E-18,
        0.142514E-17,
        0.000000E+00,
        0.156272E-17,
        0.267547E-17,
        0.300424E-17,
        0.249635E-17,
        0.389970E-17,
        0.111351E-16,
        0.000000E+00,
        0.000000E+00,
        0.240326E-17,
        0.672914E-18,
        0.980532E-18,
        0.240326E-17,
        0.156533E-18,
        0.260674E-18,
        0.342866E-19,
        0.737001E-19,
        0.111511E-18,
        0.000000E+00,
        0.817110E-18,
        0.193378E-17,
        0.368501E-17,
        0.116959E-16,
        0.000000E+00,
        0.000000E+00,
        0.145798E-17,
        0.120163E-17,
        0.579988E-19,
        0.259553E-19,
        0.000000E+00,
        0.248337E-17,
        0.248337E-17,
        0.357285E-17,
        0.338059E-17,
        0.202355E-17,
        0.259553E-17,
        0.352479E-17,
        0.202355E-17,
        0.355683E-17])
        return 6.242e+18 * activation_energies

# def print_volume_weighted_species_table(fieldz):
#     for index in range(len(fieldz[0].cIx)):
#         print("cIx = " + str(index+1), end="  ")
#         for ii in range(len(fieldz)):
#             print(str('{:0.2e}'.format(fieldz[ii].ledger[-1, index]/fieldz[ii].V[-1])), end=" & ")
#         print() #move to next line

# ============================

args = get_args()
lines = ["-.", "-.", "-.", ":" , ":", ":"]
linecycler = cycle(lines)

# fig = plt.figure(figsize=(11.8/2, 4.8))
plt.figure()
for ii in range(len(args.paths)):
    linestyle = next(linecycler)
    # If no path is given, open the default folder
    try:
        sim = Ledger(file=args.paths[ii])
    except:
        print("WARNING! No path found. Opening the default folder")
        sim = Ledger(file="/home/dennis/codes/afivo-pic/programs/combustion_3d/output/")

    # For now this is hard coded
    sim.set_grouping_IBK()
    if (args.plot_G_values):
        if (ii == 0): # Only show the labels once
            plt.plot([],[], label="N$_2$ triplets", color="k")
            plt.plot([],[], label="N", color="m")
            plt.plot([],[], label="O", color="r")
            plt.plot([],[], label="H", color="g")
            plt.plot([],[], label="O$_2$ singlets", color="tab:orange")
            plt.plot([],[], label="H$_2$", color="b")
            # plt.plot([],[], label="N$_2$ vib", color="tab:red")
            # plt.plot([],[], label="O$_2$ vib", color="tab:blue")
            # plt.plot([],[], label="CH$_4$ vib", color="tab:brown")
            # plt.legend(["N$_2$ triplets", "N", "O", "H",  "O$_2$ singlets", "H$_2$",], loc='center left', bbox_to_anchor=(1, 0.5))
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # sim.plot_E_sink(linestyle=linestyle)
        # sim.plot_IBK()
        sim.plot_G_values(linestyle=linestyle)

plt.xlabel("$L_z$ (mm)")
plt.grid(True, which="both", ls="-", color='0.85')
plt.tight_layout()

sim.print_all_G_values()



plt.show()
