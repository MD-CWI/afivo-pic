import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

#Increase font size
plt.rcParams.update({'font.size': 14})
cwd = os.getcwd()

def get_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Plot various quantities using the ledger file''')
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
        # Remove timestamps from ledger
        self.ledger = self.ledger[:, :-1]
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

    def plot_all_repo(self, factor=1, **kwargs):
        # plt.figure()
        self.plot_single_species_weighted(self.i_O_rad, np.array([2, 2, 2, 2])*factor, label="O", color='r', **kwargs)
        self.plot_single_species_weighted(self.i_H2, np.array([1, 1, 1, 2, 1, 1, 2, 1])*factor, label="H$_2$", color='b', **kwargs)
        self.plot_single_species_weighted(self.i_H_rad, np.array([1, 1, 1, 1, 1, 1])*factor, label="H", color='g', **kwargs)
        self.plot_single_species_weighted(self.i_N_rad, np.array([2, 2, 2])*factor, label="N", color='m', **kwargs)

        # plt.legend()
        plt.xlabel(self.xaxis)
        plt.ylabel("particle number")

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

    # def plot_efficiency_all_repo(self):
    #     self.plot_single_efficiency_weighted(self.i_O_rad, "O rad "+self.name, [2, 2, 2, 2])
    #     self.plot_single_efficiency_weighted(self.i_H2, "H2 "+self.name, [1, 0.5, 1, 2, 1, 1, 2, 1])
    #     self.plot_single_efficiency_weighted(self.i_H_rad, "H rad "+self.name, [1, 1, 1, 1, 1, 1, 1])
    #
    #     plt.legend()
    #     plt.xlabel(self.xaxis)
    #     plt.ylabel("eV$^{-1}$")
    #     plt.title("Efficiency (num. radicals per energy)")

    # def plot_single_efficiency_weighted(self, index, label, weight):
    #     if self.xaxis == 'length':
    #         xsteps = self.L * 1e3 # set x-axis to millimeters
    #     elif self.xaxis == 'time':
    #         xsteps = self.timesteps * 1e9 # set x-axis to nanoseconds
    #     else:
    #         print("Given xaxis not permitted. Use length or time")
    #         exit()
    #
    #     # plt.semilogy(xsteps[xsteps>=self.xmin], np.sum(self.ledger[:, index], 1)[xsteps>=self.xmin]/np.cumsum(self.P[xsteps>=self.xmin]), label=label)
    #     plt.semilogy(xsteps[xsteps>=self.xmin], np.sum(np.tile(weight, (np.shape(self.ledger)[0], 1)) * self.ledger[:, index], 1)[xsteps>=self.xmin]/self.P_integrated[xsteps>=self.xmin]/6.242e+18, label=label)

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
    def plot_single_species_weighted(self, index, weight, scale='semilog', volume_averaged=False, **kwargs):
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

        if volume_averaged:
            factor = self.V # Keep in units of m^(-3)
        else:
            factor = 1.0

        if scale == 'semilog':
            plt.semilogy(xsteps[xsteps>=self.xmin], np.sum(np.tile(weight, (np.shape(self.ledger)[0], 1)) * self.ledger[:, index], 1)[xsteps>=self.xmin]/factor[xsteps>=self.xmin], **kwargs)
        if scale == 'linear':
            plt.plot(xsteps[xsteps>=self.xmin], np.sum(np.tile(weight, (np.shape(self.ledger)[0], 1)) * self.ledger[:, index], 1)[xsteps>=self.xmin]/factor[xsteps>=self.xmin], **kwargs)


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

    def print_volume_weighted_species(self):
        # Only plot the final timestep
        for index in range(len(self.cIx)):
            print("cIx = " + str(index+1) + ":  " + str('{:0.2e}'.format(self.ledger[-1, index]/self.V[-1])))

def print_volume_weighted_species_table(fieldz):
    for index in range(len(fieldz[0].cIx)):
        print("cIx = " + str(index+1), end="  ")
        for ii in range(len(fieldz)):
            print(str('{:0.2e}'.format(fieldz[ii].ledger[-1, index]/fieldz[ii].V[-1])), end=" & ")
        print() #move to next line

# ============================

args = get_args()

# TODO Implement standardized plots as an example

plt.show()
