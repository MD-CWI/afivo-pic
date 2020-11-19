#NOTE: this script should be run on Python 2!

import sys
sys.path.append("/home/ddb/applications/visit/3.1.1/linux-x86_64/lib/site-packages") # <- Path your visit folder
import os
import visit
import numpy as np
import matplotlib.pyplot as plt

visit.Launch()

# Clear existing plots
visit.DeleteAllPlots()

# Opening a virtual database representing all wave*.silo files.
visit.OpenDatabase("/home/ddb/codes/afivo-pic/programs/combustion_2d/output/archive/CS_comparison/SUCCES/biagi_repo/sim_*.silo database") # <- Path your results

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

visit.AddPlot("Pseudocolor", "E_v2")
visit.DrawPlots()

for states in range(visit.TimeSliderGetNStates()):
    visit.SetTimeSliderState(states)
    visit.Query("Max")
    L = np.append(L, [visit.GetQueryOutputObject()['max_coord'][1]])

plt.figure()
plt.plot(t, L)
plt.xlabel('time (s)')
plt.ylabel('Length (m)')
plt.title('Streamer length')

plt.figure()
plt.plot(L, s)
plt.xlabel('Length (m)')
plt.ylabel('P_dep')

plt.figure()
plt.semilogy(L, s)
plt.xlabel('Length (m)')
plt.ylabel('P_dep')
plt.show()
