

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
# from tkinter.filedialog import askopenfilename
from scipy.optimize import curve_fit

#This is for later

# parser = argparse.ArgumentParser()
# parser.add_argument("--select", help= "Open a dialog to manually select a file",
#                                  action="store_true")
# args = parser.parse_args()
#
# if args.select:
#     from tkinter.filedialog import askopenfilename
#     filename = askopenfilename()
#     print(filename)
# else:
#     print("default")

def lin_func(x, a, b):
    return a*x+b

def exp_func(x, a, c):
    return a*np.exp(c*x)

def idx(cIx):
    "Python starts at 0, and I am not dealing with it any longer..."
    return cIx - 1 # Converts cIx to Python index

class Ledger:
    def __init__(self, file, name=""):

        self.get_streamer_length(file+"streamer_length.txt")
        self.get_streamer_radius(file+"streamer_radius.txt")

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

    def get_streamer_length(self, file):
        if os.path.isfile(file):
            self.L = np.loadtxt(file, dtype=np.dtype('float'))
        else:
            print("Cannot find file containing streamer length: " + file)
            print("First define streamer length by running write_streamer_length.py in the tools-directory")
            input("No streamer length was found. Continue anyway?")

    def get_streamer_radius(self, file):
        if os.path.isfile(file):
            self.R = np.loadtxt(file, dtype=np.dtype('float'))
        else:
            print("Cannot find file containing streamer radius: " + file)
            print("First define streamer radius by running write_streamer_length.py in the tools-directory")
            input("No streamer radius was found. Continue anyway?")

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
        self.i_O_rad = [17, 28, 30, 31]
        self.i_H2    = [34, 35, 40, 41, 42, 45, 46, 48]
        self.i_H_rad = [33, 35, 36, 40, 43, 47, 48]

    def plot_all_repo(self):
        # plt.figure()
        self.plot_single_species_weighted(self.i_O_rad, "O rad"+self.name, [2, 2, 2, 2])
        self.plot_single_species_weighted(self.i_H2, "H2"+self.name, [1, 0.5, 1, 2, 1, 1, 2, 1])
        self.plot_single_species_weighted(self.i_H_rad, "H rad"+self.name, [1, 1, 1, 1, 1, 1, 1])

        plt.legend()
        plt.xlabel("length")
        plt.ylabel("no. occurences")
        plt.title("Radicals")

    def set_grouping_biagi_and_lisbon(self):
        self.i_O_rad = [47, 57, 58]
        self.i_H2    = [62, 66]
        self.i_H_rad   = [61, 66, 67, 69]

    def plot_biagi_and_lisbon(self):
        # plt.figure()
        self.plot_single_species_weighted(self.i_O_rad, "BIAGI_IST-L O rad", [2, 2, 2])
        self.plot_single_species_weighted(self.i_H2, "BIAGI_IST-L H2", [1, 0.5])
        self.plot_single_species_weighted(self.i_H_rad, "BIAGI_IST-L H rad", [1, 1, 1, 1])

        plt.legend()
        plt.xlabel("length")
        plt.ylabel("no. occurences")
        plt.title("Radicals")

    def plot_biagi_and_lisbon_seperate(self):
        plt.figure()
        self.plot_single_species_weighted([47], "O- + O", 2.0)
        self.plot_single_species_weighted([57], "O2(EXC A3 DISOC)", 2.0)
        self.plot_single_species_weighted([58], "O2(EXC B3 DISOC)", 2.0)
        plt.title("Amount of oxygen radicals produced")
        plt.legend()

        plt.figure()
        self.plot_single_species([62], "CH2- + H2")
        self.plot_single_species_weighted([66], "CH2 + H2", 0.5)
        plt.title("Amount of molecular hydrogen produced")
        plt.legend()

        plt.figure()
        self.plot_single_species([61], "CH3 + H-")
        self.plot_single_species_weighted([66], "CH2 + H + H", 1.0)
        self.plot_single_species([67], "CH3 + H")
        self.plot_single_species([69], "CH3+ + H")
        plt.title("Amount of hydrogen radicals produced")
        plt.legend()

    def set_grouping_biagi_and_repo(self):
        self.i_O_rad = [47, 57, 58]
        self.i_H2    = [62, 63, 68, 69, 70, 73, 74, 76]
        self.i_H_rad   = [61, 63, 64, 68, 71, 75, 76]

    def plot_biagi_and_repo(self):
        # plt.figure()
        self.plot_single_species_weighted(self.i_O_rad, "BIAGI_REPO O rad", [2, 2, 2])
        self.plot_single_species_weighted(self.i_H2, "BIAGI_REPO H2", [1, 0.5, 1, 2, 1, 1, 2, 1])
        self.plot_single_species_weighted(self.i_H_rad, "BIAGI_REPO H rad", [1, 1, 1, 1, 1, 1, 1])

        plt.legend()
        plt.xlabel("length")
        plt.ylabel("no. occurences")
        plt.title("Radicals")

    def plot_biagi_and_repo_seperate(self):
        plt.figure()
        self.plot_single_species_weighted([47], "O- + O", 2.0)
        self.plot_single_species_weighted([57], "O2(EXC A3 DISOC)", 2.0)
        self.plot_single_species_weighted([58], "O2(EXC B3 DISOC)", 2.0)
        plt.title("Amount of oxygen radicals produced")
        plt.legend()

        plt.figure()
        self.plot_single_species_weighted([62], "CH2- + H2", 1)
        self.plot_single_species_weighted([63], "CH2 + H2", 0.5)
        self.plot_single_species_weighted([68], "CH + H2 + H", 1)
        self.plot_single_species_weighted([69], "C + H2 + H2", 2)
        self.plot_single_species_weighted([70], "CH2 + H2+", 1)
        self.plot_single_species_weighted([73], "CH2+ + H2", 1)
        self.plot_single_species_weighted([74], "C+ + H2 + H2", 2)
        self.plot_single_species_weighted([76], "CH+ + H + H2", 1)
        plt.title("Amount of molecular hydrogen produced")
        plt.legend()

        plt.figure()
        self.plot_single_species_weighted([61], "CH3 + H-", 1.0)
        self.plot_single_species_weighted([63], "CH2 + H + H", 1.0)
        self.plot_single_species_weighted([64], "CH3 + H", 1.0)
        self.plot_single_species_weighted([68], "CH + H2 + H", 1.0)
        self.plot_single_species_weighted([71], "CH3 + H+", 1.0)
        self.plot_single_species_weighted([75], "CH3+ + H", 1.0)
        self.plot_single_species_weighted([76], "CH+ + H + H2", 1.0)
        plt.title("Amount of hydrogen radicals produced")
        plt.legend()


    def plot_generic_grouped_species(self):
        "Plot the grouped CAS, ions/excited states/fragmented fuel/dissociated oxygen"
        plt.figure()
        self.plot_single_species(self.i_N2_ion, "N2_ion")
        self.plot_single_species(self.i_O2_ion, "O2_ion")
        self.plot_single_species(self.i_CH4_ion, "CH4_ion")

        self.plot_single_species(self.i_N2_exc, "N2_exc")
        self.plot_single_species(self.i_O2_exc, "O2_exc")
        self.plot_single_species(self.i_CH4_exc, "CH4_exc")

        self.plot_single_species(self.i_N2_attach, "N2_attach")
        self.plot_single_species(self.i_O2_attach, "O2_attach")
        self.plot_single_species(self.i_CH4_attach, "CH4_attach")

        plt.legend()
        plt.xlabel("time")
        plt.ylabel("no. occurences")
        plt.title("Grouped CAS")

    def plot_single_species(self, index, label, xaxis='length'):
        if xaxis == 'length':
            xsteps = self.L
        elif xaxis == 'time':
            xsteps = self.timesteps
        else:
            print("Given xaxis not permitted. Use length or time")
            1/0

        plt.semilogy(xsteps, np.sum(self.ledger[:, index], 1), label=label)

    def plot_single_species_weighted(self, index, label, weight, xaxis='length'):
        if xaxis == 'length':
            xsteps = self.L
        elif xaxis == 'time':
            xsteps = self.timesteps
        else:
            print("Given xaxis not permitted. Use length or time")
            1/0

        plt.semilogy(xsteps, np.sum(np.tile(weight, (np.shape(self.ledger)[0], 1)) * self.ledger[:, index], 1), label=label)

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

    def plot_single_fit(self, index, label, skip=0, logplot=False):
        # popt, pcov = curve_fit(exp_func, self.timesteps, np.squeeze(np.sum(self.ledger[:, index], 1)), p0=(1e7, 1e5))
        # yy = exp_func(self.timesteps, *popt)
        popt, pcov = curve_fit(lin_func, self.timesteps[skip:], np.log(np.squeeze(np.sum(self.ledger[skip:, index], 1))),  p0=(16, 11))
        y_fit = np.exp(lin_func(self.timesteps, *popt))

        if logplot == True:
            self.plot_single_species(index, label)
            plt.semilogy(self.timesteps, y_fit, label="fit", marker='.')
        else:
            plt.plot(self.timesteps, np.sum(self.ledger[:, index], 1), label=label)
            plt.plot(self.timesteps, y_fit, label="fit")
        plt.title("R = " + "{:.2e}".format(popt[0]) )
        plt.legend()

    def print_single_fit(self, index, skip=0):
        popt, pcov = curve_fit(lin_func, self.timesteps[skip:], np.log(np.squeeze(np.sum(self.ledger[skip:, index], 1))),  p0=(16, 11))
        print("(lambda, A) = " + str(popt))

    def print_volume_weighted_species(self, volume):
        for index in range(len(self.cIx)):
            # print(index)
            # print(self.ledger[-1, index]/volume)
            print("cIx = " + str(index+1) + ":  " + str(self.ledger[-1, index]/volume))

        # print("Total N2 ion")
        # print(self.ledger[-1, self.i_N2_ion])
        # a = np.sum(self.ledger[-1, self.i_N2_ion])
        # print(a)
        # print("Total O2 ion")
        # b = np.sum(self.ledger[-1, self.i_O2_ion])
        # print(b)
        # print("Total CH4 ion")
        # c = np.sum(self.ledger[-1, self.i_CH4_ion])
        #
        # print(c)
        #
        # print("sum pos ion")
        # print(a + b + c)
        #
        #
        # print("Total attach:")
        # d = np.sum(self.ledger[-1, self.i_N2_attach]) + np.sum(self.ledger[-1, self.i_O2_attach]) + np.sum(self.ledger[-1, self.i_CH4_attach])
        # # d = d/volume
        # print(d)
        # print("netto pos ion")
        # print((a + b + c - d))

        # print(self.i_N2_attach)
        # print(self.i_O2_attach)
        # print(self.i_CH4_attach)

# ==============================
fuel = Ledger("/home/dennis/codes/afivo-pic/programs/combustion_3d/output/archive/maybe_sim_thijs/")
# fuel = Ledger("/home/dennis/Documents/drafts/air_methane_streamers/results/bigboy/fuel/")
fuel.print_volume_weighted_species(volume=1)

fuel.set_grouping_all_repo()
fuel.plot_all_repo()

# ========

# fuel = Ledger("/home/dennis/Documents/drafts/air_methane_streamers/results/electrode_fat/fuel/")
# air = Ledger("/home/dennis/Documents/drafts/air_methane_streamers/results/electrode_fat/air/")
#
# plt.figure()
# plt.plot(air.L[-len(air.R):]*1e3, air.R*1e3, 'rd-', label='air')
# plt.plot(fuel.L[-len(fuel.R):]*1e3, fuel.R*1e3, 'ks-', label='fuel')
# plt.xlabel('length (mm)')
# plt.ylabel('radius (mm)')
# plt.legend()


# fuel = Ledger("/home/ddb/results/electrode_fat/fuel/")
# air = Ledger("/home/ddb/results/electrode_fat/air/")

# fuel.set_grouping_all_repo()
# fuel.plot_all_repo()
# plt.figure()
# air.set_grouping_all_repo()
# air.plot_all_repo()
#
# xn_fuel = (fuel.L[:-1] + fuel.L[1:])/2
# xn_air = (air.L[:-1] + air.L[1:])/2
#
# plt.figure()
# plt.plot(xn_fuel/1e-3, np.diff(fuel.L)/np.diff(fuel.timesteps)/1e6, label='fuel')
# plt.plot(xn_air/1e-3, np.diff(air.L)/np.diff(air.timesteps)/1e6, label='air')
# plt.xlabel('length (mm)')
# plt.ylabel('velocity (mm/ns)')
# plt.legend()

# =======
# plt.figure()
# plt.plot(fuel.timesteps/1e-9, fuel.L/1e-3, label='fuel')
# plt.plot(air.timesteps/1e-9, air.L/1e-3, label='air')
# plt.xlabel('time (ns)')
# plt.ylabel('Length (mm)')
# plt.legend()
#


# lisbon.set_grouping_biagi_and_lisbon()
# repo.set_grouping_biagi_and_repo()
#
# full_repo.plot_generic_grouped_species()
# repo.plot_generic_grouped_species()

# plt.figure()
# lisbon.plot_biagi_and_lisbon()
# repo.plot_biagi_and_repo()
#
# # plt.figure()
# # lisbon.plot_biagi_and_lisbon_seperate()
# repo.plot_biagi_and_repo_seperate()
#
#
# # %% Now some plotting
# # ledger.set_grouping_biagi_and_lisbon()
# # ledger.plot_biagi_and_lisbon_seperate()
# # ledger.generic_grouping()
# #
# repo.plot_generic_grouped_species()
# # #
# repo.plot_excitations(repo.i_N2_exc, "N2")
# repo.plot_excitations(repo.i_O2_exc, "O2")
# repo.plot_excitations(repo.i_CH4_exc, "CH4")

plt.show()
