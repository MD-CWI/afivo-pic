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

class Ledger:
    def __init__(self, dt=5e-11):
        self.dt = 2.5e-10 # Output time step

        # Now a set of indices, TODO: should be automatic and supplied by "file"
        self.i_N2_ion    = [24]
        self.i_O2_ion    = np.arange(39, 41)
        self.i_CH4_ion   = [51]

        self.i_N2_exc    = np.arange(1, 24)
        self.i_O2_exc    = np.arange(26, 39)# [26:38]
        self.i_CH4_exc   = np.arange(45, 51) #[45:50]

        self.i_CH4_frag    = [42, 43, 52]
        self.i_O_diss      = [41]


    def get_file(self, file):
        self.ledger = np.genfromtxt(file)
        N = np.size(self.ledger, 0)
        self.timesteps = np.linspace(0, N, N) * self.dt

    def plot_grouped_species(self):
        "Plot the grouped CAS, ions/excited states/fragmented fuel/dissociated oxygen"
        plt.figure()
        self.plot_single_species(self.i_N2_ion, "N2_ion" )
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

    def plot_single_fit(self, index, label, skip=0, logplot=False):
        # popt, pcov = curve_fit(exp_func, self.timesteps, np.squeeze(np.sum(self.ledger[:, index], 1)), p0=(1e7, 1e5))
        # yy = exp_func(self.timesteps, *popt)
        popt, pcov = curve_fit(lin_func, self.timesteps[skip:], np.log(np.squeeze(np.sum(self.ledger[skip:, index], 1))),  p0=(16, 11))
        y_fit = np.exp(lin_func(self.timesteps, *popt))

        if logplot == True:
            self.plot_single_species(index, label)
            plt.semilogy(self.timesteps, y_fit, label=str(index), marker='o')
        else:
            plt.plot(self.timesteps, np.sum(self.ledger[:, index], 1), label=label)
            plt.plot(self.timesteps, y_fit, label="fit")
        plt.title( "{:.2e}".format(popt[1]) )
        plt.legend()

    def print_single_fit(self, index, skip=0):
        popt, pcov = curve_fit(lin_func, self.timesteps[skip:], np.log(np.squeeze(np.sum(self.ledger[skip:, index], 1))),  p0=(16, 11))
        # popt[0] = np.exp(popt[0])
        print("(A, lambda) = " + str(popt))

# ==============================

ledger = Ledger()
ledger.get_file("/home/ddb/codes/afivo-pic/programs/standard_3d/output/sim_werktdit.txt")
# ledger.plot_grouped_species()

plt.figure()
# for ii in range(1, 24):
#     ledger.print_single_fit([ii])
#     ledger.plot_single_fit([ii], str(ii), logplot=True)

ledger.plot_single_fit(ledger.i_N2_exc, "N2_exc", skip=20, logplot=True)
ledger.plot_single_fit(ledger.i_O2_exc, "O2_exc", skip=20, logplot=True)
ledger.plot_single_fit(ledger.i_CH4_exc, "CH4_exc", skip=20, logplot=True)
ledger.plot_single_fit(ledger.i_N2_ion, "N2_ion", skip=20, logplot=True)
ledger.plot_single_fit(ledger.i_O2_ion, "O2_ion", skip=20, logplot=True)
ledger.plot_single_fit(ledger.i_CH4_ion, "CH4_ion", skip=20, logplot=True)
ledger.plot_single_fit(ledger.i_CH4_frag, "CH4_frag", skip=20, logplot=True)
ledger.plot_single_fit(ledger.i_O_diss, "O_diss", skip=20, logplot=True)
# ledger.plot_single_species(ledger.i_N2_ion, "N2_ion")
plt.show()
