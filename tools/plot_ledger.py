import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from tkinter.filedialog import askopenfilename
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
    def __init__(self, file):

        # Header data is loaded seperately
        self.cIx = np.loadtxt(file, skiprows = 4, max_rows = 1, dtype=np.dtype('i4'))
        self.gas = np.loadtxt(file, skiprows = 5, max_rows = 1, dtype=np.dtype('U'))
        self.coll_type = np.loadtxt(file, skiprows = 6, max_rows = 1, dtype=np.dtype('U'))

        # now load collision data to ledger
        self.ledger = np.loadtxt(file, skiprows = 7)
        # Extract timestamps
        self.timesteps = self.ledger[:, -1]
        # Remove timestamps from ledger
        self.ledger = self.ledger[:, :-1]
        # Assuming equal timespacing calculate dt
        self.dt = self.timesteps[1] - self.timesteps[0]

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

    def set_manual_grouping(self):
        " This grouping of species keeps the dissocation and fragmentation of ..."
        "... O2 and CH4 seperate (respectively)."
        self.i_N2_ion    = [24]
        self.i_O2_ion    = np.arange(39, 41)
        self.i_CH4_ion   = [51]

        self.i_N2_exc    = np.arange(1, 24)
        self.i_O2_exc    = np.arange(26, 39)
        self.i_CH4_exc   = np.arange(45, 51)

        self.i_CH4_frag    = [42, 43, 52]
        self.i_O_diss      = [41]

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

    def plot_grouped_species(self):
        "Plot the grouped CAS, ions/excited states/fragmented fuel/dissociated oxygen"
        plt.figure()
        self.plot_single_species(self.i_N2_ion, "N2_ion")
        self.plot_single_species(self.i_O2_ion, "O2_ion")
        self.plot_single_species(self.i_CH4_ion, "CH4_ion")

        self.plot_single_species(self.i_N2_exc, "N2_exc")
        self.plot_single_species(self.i_O2_exc, "O2_exc")
        self.plot_single_species(self.i_CH4_exc, "CH4_exc")

        self.plot_single_species([self.i_O_diss, self.i_O_diss], "O_diss")
        self.plot_single_species(self.i_CH4_frag, "CH4_frag")

        plt.legend()
        plt.xlabel("time")
        plt.ylabel("no. occurences")
        plt.title("Grouped CAS")

    def plot_grouped_species_derivative(self):
        "Plot the derivative of grouped CAS, ions/excited states/fragmented fuel/dissociated oxygen"
        plt.figure()
        plt.semilogy(self.timesteps, np.gradient(np.sum(self.ledger[:, self.i_N2_ion], 1), self.dt), label="d_N2_ion")
        plt.semilogy(self.timesteps, np.gradient(np.sum(self.ledger[:, self.i_O2_ion], 1), self.dt), label="d_O2_ion")
        plt.semilogy(self.timesteps, np.gradient(np.sum(self.ledger[:, self.i_CH4_ion],1), self.dt), label="d_CH4_ion")

        plt.semilogy(self.timesteps, np.gradient(np.sum(self.ledger[:, self.i_N2_exc], 1), self.dt),
                            label="d_N2_exc")
        plt.semilogy(self.timesteps, np.gradient(np.sum(self.ledger[:, self.i_O2_exc], 1), self.dt),
                            label="d_O2_exc")
        plt.semilogy(self.timesteps, np.gradient(np.sum(self.ledger[:, self.i_CH4_exc], 1), self.dt),
                            label="d_CH4_exc")

        plt.semilogy(self.timesteps, np.gradient(np.sum(2 * self.ledger[:, self.i_O_diss], 1), self.dt),
                            label="d_O_diss") # One occasion produces two dissociated oxygen
        plt.semilogy(self.timesteps, np.gradient(np.sum(self.ledger[:, self.i_CH4_frag], 1), self.dt),
                            label="d_CH4_frag")

        plt.legend()
        plt.xlabel("time")
        plt.ylabel("no. occurences")
        plt.title("Derivative of grouped CAS")

    def plot_single_species(self, index, label):
        plt.semilogy(self.timesteps, np.sum(self.ledger[:, index], 1), label=label)

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


# ==============================

ledger = Ledger("/home/ddb/codes/afivo-pic/programs/combustion_2d/output/sim_cs_ledger.txt")
ledger.set_manual_grouping()

# %% Now some plotting

ledger.plot_grouped_species()

ledger.plot_excitations(ledger.i_N2_exc, "N2")
ledger.plot_excitations(ledger.i_O2_exc, "O2")
ledger.plot_excitations(ledger.i_CH4_exc, "CH4")

plt.show()
