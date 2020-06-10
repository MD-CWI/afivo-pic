import sys
sys.path.append("/home/ddb/applications/visit/3.1.1/linux-x86_64/lib/site-packages") # <- Path your visit folder
import os
import visit
import numpy as np
import matplotlib.pyplot as plt

# Launch visit
visit.Launch()

# Clear existing plots
visit.DeleteAllPlots()

# Opening a virtual database representing all *.silo files.
visit.OpenDatabase("/home/ddb/codes/afivo-pic/programs/dielectric_2d/output/test_*.silo database") # <- Path your results + " database"

# Make plots
visit.AddPlot("Pseudocolor", "O_atom") # Type of plot || variable to plot
visit.DrawPlots()

s = []
t = []

# loop over time
for states in range(visit.TimeSliderGetNStates()):
    visit.SetTimeSliderState(states)
    visit.Query("Weighted Variable Sum") # Select which type of query to do
    s = np.append(s, [visit.GetQueryOutputValue()]) #save query-result to numpy array
    visit.Query("Time") # Query to get the current time
    t = np.append(t, [visit.GetQueryOutputValue()])

# Make plots using matplotlib as usual
plt.figure()
plt.plot(t, s)
plt.xlabel('time (s)')
plt.ylabel('#Reactions $m^{-2}$ (?) producing $O^-$')
plt.show()
