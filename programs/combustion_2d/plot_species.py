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
visit.OpenDatabase("/home/ddb/codes/afivo-pic/programs/dielectric_2d/output/test_*.silo database") # <- Path your results

# Make plots
visit.AddPlot("Pseudocolor", "O_atom")
visit.DrawPlots()

s = []
t = []

# plot over time
for states in range(visit.TimeSliderGetNStates()):
    visit.SetTimeSliderState(states)
    visit.Query("Weighted Variable Sum")
    s = np.append(s, [visit.GetQueryOutputValue()])
    visit.Query("Time")
    t = np.append(t, [visit.GetQueryOutputValue()])


plt.figure()
plt.plot(t, s)
plt.xlabel('time (s)')
plt.ylabel('#Reactions $m^{-2}$ (?) producing $O^-$')
plt.show()
