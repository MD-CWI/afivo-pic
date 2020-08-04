import numpy as np
import matplotlib.pyplot as plt
import os

# from tkinter.filedialog import askopenfilename
# filename = askopenfilename()
# print(filename)

file = "/home/ddb/results/PIC_3D_airfuel_5e6_double/sim_werktdit.txt"
# file = os.getcwd() + "/output/archive/ledger_phi_25/test_werktdit.txt"
# file = os.getcwd() + "/output/archive/ledger_phi_3/test_werktdit.txt"
# file = os.getcwd() + "/output/archive/ledger_phi_35/test_werktdit.txt"

ledger = np.genfromtxt(file)
N = np.size(ledger, 0)
dt_out = 2.5e-10
timesteps = np.linspace(0, N, N) * dt_out

N2_ion  = ledger[:, 24]
O2_ion  = ledger[:, 39] + ledger[:, 40]
CH4_ion = ledger[:, 51]

N2_exc  = np.sum(ledger[:, 1:23], 1)
O2_exc  = np.sum(ledger[:, 26:38], 1)
CH4_exc = np.sum(ledger[:, 45:50], 1)

CH4_frag = ledger[:, 42] + ledger[:, 43] + ledger[:, 52]
O_diss   = 2 * ledger[:, 41]

d_N2_ion  = np.gradient(ledger[:, 24], dt_out)
d_O2_ion  = np.gradient(ledger[:, 39] + ledger[:, 40], dt_out)
d_CH4_ion = np.gradient(ledger[:, 51], dt_out)

d_N2_exc  = np.gradient(np.sum(ledger[:, 1:23], 1), dt_out)
d_O2_exc  = np.gradient(np.sum(ledger[:, 26:38], 1), dt_out)
d_CH4_exc = np.gradient(np.sum(ledger[:, 45:50], 1), dt_out)

d_CH4_frag = np.gradient(ledger[:, 42] + ledger[:, 43] + ledger[:, 52], dt_out)
d_O_diss   =np.gradient( 2 * ledger[:, 41], dt_out)

N2_exc_total = np.tile(np.sum(ledger[:, 1:23], 1, keepdims=True), (1, 22))
N2_exc_relative = np.divide(ledger[:, 1:23], N2_exc_total+1e-10)

plt.figure()
plt.stackplot(timesteps, np.transpose(N2_exc_relative), labels=np.arange(2,24).astype(str))
plt.legend()
#
# plt.figure()
# plt.stackplot(timesteps, np.transpose(ledger[:, 1:23]))
# plt.legend()
# plt.figure()
# plt.semilogy(timesteps, N2_exc_relative)
#
# plt.figure()
# plt.semilogy(timesteps, N2_ion, label="N2_ion")
# plt.semilogy(timesteps, O2_ion, label="O2_ion")
# plt.semilogy(timesteps, CH4_ion, label="CH4_ion")
# plt.legend()
#
# # plt.figure()
# plt.semilogy(timesteps, N2_exc, label="N2_exc")
# plt.semilogy(timesteps, O2_exc, label="O2_exc")
# plt.semilogy(timesteps, CH4_exc, label="CH4_exc")
# plt.legend()
#
# # plt.figure()
# plt.semilogy(timesteps, CH4_frag, label="CH4_frag")
# plt.semilogy(timesteps, O_diss, label="O_diss")
# plt.xlabel("time")
# plt.ylabel("# occurences")
# plt.legend()
# #
# plt.figure()
# plt.semilogy(timesteps, d_N2_ion, label="d_N2_ion")
# plt.semilogy(timesteps, d_O2_ion, label="d_O2_ion")
# plt.semilogy(timesteps, d_CH4_ion, label="d_CH4_ion")
# plt.legend()
#
# # plt.figure()
# plt.semilogy(timesteps, d_N2_exc, label="d_N2_exc")
# plt.semilogy(timesteps, d_O2_exc, label="d_O2_exc")
# plt.semilogy(timesteps, d_CH4_exc, label="d_CH4_exc")
# plt.legend()
#
# # plt.figure()
# plt.semilogy(timesteps, d_CH4_frag, label="d_CH4_frag")
# plt.semilogy(timesteps, d_O_diss, label="d_O_diss")
# plt.xlabel("time")
# plt.ylabel("# occurences")
# plt.legend()
#
#
plt.show()
