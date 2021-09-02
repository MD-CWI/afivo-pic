#NOTE: this script should be run on Python 2!

import sys
sys.path.append("/opt/visit/3.1.2/linux-x86_64/lib/site-packages") # <- Path your visit folder
import os
import visit
import numpy as np
import matplotlib.pyplot as plt

visit.Launch()

# Clear existing plots
visit.DeleteAllPlots()

# Opening a virtual database representing all wave*.silo files.
# visit.OpenDatabase("/home/dennis/Documents/drafts/air_methane_streamers/results/bigboy/fuel/sim_*.silo database") # <- Path your results
visit.OpenDatabase("/home/dennis/Documents/drafts/air_methane_streamers/results/new_acc_long/sim_*.silo database") # <- Path your results

# Make plots
visit.AddPlot("Pseudocolor", "power_deposition")
visit.DrawPlots()

s = []
t = []
L = []

# plot over time
for states in range(visit.TimeSliderGetNStates()):
    visit.SetTimeSliderState(states)
    visit.Query("Weighted Variable Sum")
    s = np.append(s, [visit.GetQueryOutputValue()])
    visit.Query("Time")
    t = np.append(t, [visit.GetQueryOutputValue()])

# visit.AddPlot("Pseudocolor", "E_v2")
# visit.DrawPlots()
#
# for states in range(visit.TimeSliderGetNStates()):
#     visit.SetTimeSliderState(states)
#     visit.Query("Max")
#     L = np.append(L, [visit.GetQueryOutputObject()['max_coord'][-1]])
# L = (1e-2 - L) - (1e-2 - L[0])

# #####
# # Opening a virtual database representing all wave*.silo files.
# visit.OpenDatabase("/home/ddb/results/electrode_fat/air/sim_*.silo database") # <- Path your results
#
# # Make plots
# visit.AddPlot("Pseudocolor", "power_deposition")
# visit.DrawPlots()
#
# L_air = []
# s_air = []
#
# # plot over time
# for states in range(visit.TimeSliderGetNStates()):
#     visit.SetTimeSliderState(states)
#     visit.Query("Weighted Variable Sum")
#     s_air = np.append(s_air, [visit.GetQueryOutputValue()])
#
# visit.AddPlot("Pseudocolor", "E_v2")
# visit.DrawPlots()
#
# for states in range(visit.TimeSliderGetNStates()):
#     visit.SetTimeSliderState(states)
#     visit.Query("Max")
#     L_air = np.append(L_air, [visit.GetQueryOutputObject()['max_coord'][-1]])
# L_air = (1e-2 - L_air) - (1e-2 - L_air[0])
#
# ###


# plt.figure()
# plt.plot(t, L)
# plt.xlabel('time (s)')
# plt.ylabel('Length (m)')
# plt.title('Streamer length')

plt.figure()
plt.plot(t, s)
plt.xlabel('time (s)')
plt.ylabel('P_dep (J/s)')

plt.figure()
plt.semilogy(t, s)
plt.xlabel('time (s)')
plt.ylabel('P_dep (J/s)')

s[0]=0

s_trapz = np.trapz(s[0:81], t[0:81])

print("Total energy: ")
print(np.sum(s) * (t[1]-t[0]))
print("Traps: "+str(s_trapz))
print(t)
print(s)

print("Final timeframe = "+str(s[-1]))

# plt.figure()
# plt.semilogy(L/1e-3, s, label='fuel')
# plt.semilogy(L_air/1e-3, s_air, label='air')
# plt.xlabel('Length (mm)')
# plt.ylabel('Power deposition (J/s)')
# plt.legend()
plt.show()
