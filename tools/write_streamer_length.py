#NOTE: this script should be run on Python 2!

# A post-processing tool to write the last coordinate (=vertical) maximum to a text file

import sys
sys.path.append("/opt/visit/3.1.2/linux-x86_64/lib/site-packages") # <- Path your visit folder
import os
import visit
import numpy as np
import matplotlib.pyplot as plt

def convert_to_length(z_coord):
    return (z_tip - z_coord)

output_path = "/home/dennis/Documents/drafts/air_methane_streamers/results/electrode_fat/air/"

visit.Launch()

# Clear existing plots
visit.DeleteAllPlots()

# Opening a virtual database representing all wave*.silo files.
visit.OpenDatabase(output_path + "sim_*.silo database") # <- Path your results

# Make plots
visit.AddPlot("Pseudocolor", "E_v2")
visit.DrawPlots()

L = []

# First time we pass through the data to get the length
for states in range(visit.TimeSliderGetNStates()):
    visit.SetTimeSliderState(states)
    visit.Query("Max")
    max_coord = visit.GetQueryOutputObject()['max_coord'][-1]

    if states == 0:
        visit.Query("SpatialExtents")
        x_max = float(visit.GetQueryOutputValue()[1])
        y_max = float(visit.GetQueryOutputValue()[3])
        z_max = float(visit.GetQueryOutputValue()[5])
        z_tip = max_coord # Assume initially, the field is highest at the electrode tip

    current_length = convert_to_length(max_coord)
    L = np.append(L, [current_length])

np.savetxt(output_path + "streamer_length.txt", L)

R = []
start_state = 10
N_slices = 25  # number of slices per mm
slice_size = 1e-3/N_slices
visit.AddOperator("Slice")
for states in range(start_state, visit.TimeSliderGetNStates()): # We start at a later state to ensure inception already occurred
    visit.SetTimeSliderState(states)
    num_steps = max(int((L[states] + 0.25e-3) / slice_size), 10)
    Rx = []
    Ry = []

    for num_slice in range(max(int(num_steps/2), 1), num_steps):
        s = visit.SliceAttributes()
        s.originType = s.Intercept
        s.originIntercept = z_tip - (num_slice * slice_size)
        s.axisType = 2
        s.project2d = 0
        visit.SetOperatorOptions(s)
        visit.DrawPlots()
        # Do the queries here
        visit.Query("Max")
        Rx = np.append(Rx, [[visit.GetQueryOutputObject()['max_coord'][0]]])
        Ry = np.append(Ry, [[visit.GetQueryOutputObject()['max_coord'][1]]])

        #Avoid detecting fields past some threshold (i.e. if the maximum field is on a boundary)
        for ii in range(len(Rx)):
            if np.sqrt((Rx[ii]-0.5*x_max)**2 + (Ry[ii]-0.5*y_max)**2) > 3e-3:
                Rx[ii] = 0.5*x_max
                Ry[ii] = 0.5*y_max

    R = np.append(R, [np.max(np.sqrt((Rx-0.5*x_max)**2 + (Ry-0.5*y_max)**2))])


plt.figure()
plt.plot(L[start_state:]*1e3, R*1e3)
plt.ylabel("Radius (mm)")
plt.xlabel("Length (mm)")

plt.show()

np.savetxt(output_path + "streamer_radius.txt", R)
