#NOTE: this script should be run on Python 2!

# A post-processing tool to write streamer properties to a text file

import sys
sys.path.append("/opt/visit/3.1.2/linux-x86_64/lib/site-packages") # <- Path your visit folder
import os
import visit
import numpy as np
import matplotlib.pyplot as plt
import argparse

def get_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Compute streamer properties and write them to a text file''')
    p.add_argument('-write_all', dest='write_all', action='store_true',
                   help='write parameters to file')
    p.add_argument('-write_length', dest='write_length', action='store_true',
                   help='write the distance of the maximal electric field to the tip of the electrode (only valid for axisymmetric streamers)')
    p.add_argument('-write_Emax', dest='write_Emax', action='store_true',
                   help='write the maximal electric field')
    p.add_argument('-write_radius', dest='write_radius', action='store_true',
                   help='write the electrodynamic radius (only for axisymmetric streamers)')
    p.add_argument('-R_start_state', dest='R_start_state', type=int, default=10,
                   help='Skip the first N time states to ensure R is calculated only after inception')
    p.add_argument('-R_slices_per_mm', dest='R_slices_per_mm', type=float, default=25,
                   help='Number of horizontal slices per mm considered for determining radius')
    p.add_argument('-write_volume', dest='write_volume', action='store_true',
                   help='write the volume of a streamer where the electron density is above a threshold value')
    p.add_argument('-volume_threshold', type=float, default=1.0e19,
                   help='The threshold value of the electron density that defines the streamer volume')
    p.add_argument('-write_Lgap', dest='write_Lgap', action='store_true',
                   help='write the max z coordinate where the electron density is above the volume threshold value')
    p.add_argument('-write_energy_deposition', dest='write_energy_deposition', action='store_true',
                   help='write the time- and space-integrated energy density deposition')
    p.add_argument('-no_win', dest='no_win', action='store_true',
                   help='suppress visit to open windows')
    p.add_argument('-path', dest='path', type=str,
                   help='The absolute path to the silo files')
    p.add_argument('-database_name', dest='database_name', type=str, default="sim_*.silo database",
                  help='The absolute path to the silo files')
    p.add_argument('-no_plots', dest='no_plots', action='store_true',
                   help='Do not show plots of the calculated variables')
    return p.parse_args()

def convert_to_length(z_coord):
    return (z_tip - z_coord)

def initialize_parameters():
    visit.RemoveAllOperators()
    visit.DeleteAllPlots()

    # Make plots
    visit.AddPlot("Pseudocolor", "E")
    visit.DrawPlots()

    visit.SetTimeSliderState(0)
    visit.Query("Max")
    # max_val = visit.GetQueryOutputValue()
    max_coord = visit.GetQueryOutputObject()['max_coord'][-1]

    visit.Query("SpatialExtents")
    x_max = float(visit.GetQueryOutputValue()[1])
    y_max = float(visit.GetQueryOutputValue()[3])
    z_max = float(visit.GetQueryOutputValue()[5])
    z_tip = max_coord # Assume initially, the field is highest at the electrode tip

    return x_max, y_max, z_max, z_tip

def calculate_time_length_and_Emax():
    print("Calculating L and Emax")
    visit.RemoveAllOperators()
    visit.DeleteAllPlots()

    # Make plots
    visit.AddPlot("Pseudocolor", "E")
    visit.DrawPlots()

    t = []
    L = []
    Emax = []

    # First time we pass through the data to get the length
    for states in range(visit.TimeSliderGetNStates()):
        print(str(states+1) + " out of " +  str(visit.TimeSliderGetNStates()))
        visit.SetTimeSliderState(states)
        visit.Query("Max")
        max_val = visit.GetQueryOutputValue()
        max_coord = visit.GetQueryOutputObject()['max_coord'][-1]

        current_length = convert_to_length(max_coord)
        t = np.append(t, [states])
        L = np.append(L, [current_length])
        Emax = np.append(Emax, [max_val])

    np.savetxt(args.path + "/streamer_length.txt", L)
    np.savetxt(args.path + "/streamer_Emax.txt", Emax)
    plot_length(t, L)
    plot_Emax(t, Emax)

    return t, L, Emax

def calculate_R(L):
    print("Calculating R")
    visit.RemoveAllOperators()
    visit.DeleteAllPlots()

    visit.AddPlot("Pseudocolor", "E")
    visit.DrawPlots()

    R = []
    start_state = args.R_start_state
    N_slices = args.R_slices_per_mm  # number of slices per mm
    slice_size = 1e-3/N_slices
    visit.AddOperator("Slice")
    for states in range(start_state, visit.TimeSliderGetNStates()): # We start at a later state to ensure inception already occurred
        print(str(states+1) + " out of " +  str(visit.TimeSliderGetNStates()))
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

    np.savetxt(args.path + "/streamer_radius.txt", R)
    plot_R(t, R)
    return R

def calculate_energy():
    print("Calculating energy deposition")
    visit.RemoveAllOperators()
    visit.DeleteAllPlots()

    visit.AddPlot("Pseudocolor", "energy_deposition")
    visit.DrawPlots()

    E_dep = []
    t = []

    for states in range(visit.TimeSliderGetNStates()):
        print(str(states+1) + " out of " +  str(visit.TimeSliderGetNStates()))
        visit.SetTimeSliderState(states)
        visit.Query("Weighted Variable Sum")

        t = np.append(t, [states])
        E_dep = np.append(E_dep, [visit.GetQueryOutputValue()])

    np.savetxt(args.path + "/streamer_energy.txt", E_dep)
    plot_E_dep(t, E_dep)
    return t, E_dep

def calculate_volume():
    print("Calculating volume")
    visit.RemoveAllOperators()
    visit.DeleteAllPlots()

    visit.AddPlot("Pseudocolor", "electron")
    visit.DefineScalarExpression("I_volume", "if(ge(electron,"+str(args.volume_threshold)+"), I_volume=1, I_volume=0)") # Indicator function defining the volume
    visit.ChangeActivePlotsVar("I_volume")
    visit.DrawPlots()

    V = []
    t = []

    for states in range(visit.TimeSliderGetNStates()):
        print(str(states+1) + " out of " +  str(visit.TimeSliderGetNStates()))
        visit.SetTimeSliderState(states)
        visit.Query("Weighted Variable Sum") # This query derives the volume of the streamer

        t = np.append(t, [states])
        V = np.append(V, [visit.GetQueryOutputValue()])

    np.savetxt(args.path + "/streamer_volume.txt", V)
    plot_volume(t, V)
    return t, V

def calculate_Lgap():
    print("Calculating Lgap")
    visit.RemoveAllOperators()
    visit.DeleteAllPlots()

    visit.AddPlot("Contour", "electron")
    visit.DefineScalarExpression("I_mask", "if(ge(electron,"+str(args.volume_threshold)+"), I_mask=1, I_mask=0)") # Indicator function defining the mask
    visit.ChangeActivePlotsVar("I_mask")
    visit.DrawPlots()

    c = visit.ContourAttributes()
    c.contourNLevels = 1
    c.contourValue = (1.0)
    c.max = 2.0
    visit.SetPlotOptions(c)

    Lgap = []
    t = []

    for states in range(visit.TimeSliderGetNStates()):
        print(str(states+1) + " out of " +  str(visit.TimeSliderGetNStates()))
        visit.SetTimeSliderState(states)
        visit.Query("SpatialExtents") # This query derives the mask of the streamer
        
        t = np.append(t, [states])
        Lgap = np.append(Lgap, [visit.GetQueryOutputValue()[5] - visit.GetQueryOutputValue()[4]])

    np.savetxt(args.path + "/streamer_Lgap.txt", Lgap)
    plot_Lgap(t, Lgap)

    return t, Lgap

def plot_length(t, L):
    plt.figure()
    plt.plot(t, L*1e3)
    plt.xlabel("time (frame)")
    plt.ylabel("length (mm)")
    plt.title("Length")
    plt.tight_layout()

def plot_Emax(t, Emax):
    plt.figure()
    plt.plot(t, Emax)
    plt.xlabel("time (frame)")
    plt.ylabel("electric field (V/m)")
    plt.title("Maximum electric field")
    plt.tight_layout()

def plot_R(t, R):
    plt.figure()
    plt.plot(t[args.R_start_state:], R*1e3)
    plt.ylabel("radius (mm)")
    plt.xlabel("time (mm)")
    plt.title("Radius")
    plt.tight_layout()

def plot_E_dep(t, E_dep):
    plt.figure()
    plt.plot(t, E_dep)
    plt.xlabel("time (frame)")
    plt.ylabel("Energy deposition (J)")
    plt.title("Energy deposition")
    plt.tight_layout()

def plot_volume(t, V):
    plt.figure()
    plt.plot(t, V*1e-9)
    plt.xlabel("time (frame)")
    plt.ylabel("volume (mm$^3$)")
    plt.title("Volume")
    plt.tight_layout()

def plot_Lgap(t, Lgap):
    plt.figure()
    plt.plot(t, Lgap*1e-3)
    plt.xlabel("time (frame)")
    plt.ylabel("Lgap (mm)")
    plt.title("axial length (i.e. Lgap)")
    plt.tight_layout()


args = get_args()
z_tip = 0 # needs to be initialized

if (not args.no_win):
    visit.Launch()
else:
    visit.LaunchNowin()

# Opening a virtual database representing all wave*.silo files.
visit.OpenDatabase(args.path +"/"+ args.database_name) # <- Path your results

[x_max, y_max, z_max, z_tip] = initialize_parameters()

if (args.write_radius or args.write_all):
    [t, L, Emax] = calculate_time_length_and_Emax()
    R = calculate_R(L)
elif (args.write_length or args.write_Emax):
    [t, L, Emax] = calculate_time_length_and_Emax()

if (args.write_energy_deposition or args.write_all): [t, E_dep] = calculate_energy()
if (args.write_volume or args.write_all):                [t, V] = calculate_volume()
if (args.write_Lgap or args.write_all):               [t, Lgap] = calculate_Lgap()

if (not args.no_plots): plt.show()
