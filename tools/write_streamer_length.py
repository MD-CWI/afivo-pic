#NOTE: this script should be run on Python 2!

# A post-processing tool to write the last coordinate (=vertical) maximum to a text file

import sys
sys.path.append("/home/ddb/applications/visit/3.1.1/linux-x86_64/lib/site-packages") # <- Path your visit folder
import os
import visit
import numpy as np

output_path = "/home/ddb/results/electrode_fat/fuel/"

visit.Launch()

# Clear existing plots
visit.DeleteAllPlots()

# Opening a virtual database representing all wave*.silo files.
visit.OpenDatabase(output_path + "sim_*.silo database") # <- Path your results

# Make plots
visit.AddPlot("Pseudocolor", "E_v2")
visit.DrawPlots()

L = []

for states in range(visit.TimeSliderGetNStates()):
    visit.SetTimeSliderState(states)
    visit.Query("Max")
    L = np.append(L, [visit.GetQueryOutputObject()['max_coord'][-1]])

L = (1e-2 - L) - (1e-2 - L[0])
np.savetxt(output_path + "streamer_length.txt", L)
