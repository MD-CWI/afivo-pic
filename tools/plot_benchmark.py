import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline

file1 = "/home/ddb/cartesius/output/sim_diagnostics"
file2 = "/home/ddb/results/benchmark/sim_diagnostics_C24_1"
file5 = "/home/ddb/results/benchmark/sim_diagnostics_C24_2"

file3 = "/home/ddb/results/benchmark/laptop/test_diagnostics_PHOTOI"
file4 = "/home/ddb/results/benchmark/laptop/test_diagnostics_vanilla"
file000 = "/home/ddb/codes/afivo-pic/programs/combustion_2d/output/sim_diagnostics"

# files = [file1, file2, file5]
files = [file000, file3, file4]

for ii in range(len(files)):
    dat = np.loadtxt(files[ii])

    it            = dat[:, 0]
    n_part        = dat[:, 1]
    wc_time       = dat[:, 2]
    wtime_advance = dat[:, 3] # Watch out: these are CUMULATIVE
    wtime_run     = dat[:, 4] # Watch out: these are CUMULATIVE
    wpercent      = dat[:, 5]

    plt.figure(0)
    plt.plot(it, wc_time)
    plt.xlabel('#iterations')
    plt.ylabel('time')
    plt.title('RACE')

    plt.figure(1)
    plt.semilogy(it, n_part)
    plt.xlabel('#iterations')
    plt.ylabel('#particles')
    plt.title('Particles over time')

    plt.figure(2)
    plt.semilogy(it[1:], np.diff(wc_time))
    plt.xlabel('#iterations')
    plt.ylabel('time')
    plt.title('Calculation time of one iteration')

    plt.figure(3)
    plt.loglog(n_part[1:], np.diff(wc_time))
    plt.xlabel('#particles')
    plt.ylabel('time')
    plt.title('Calculation time of one iteration')

    # plt.figure(4)
    # plt.semilogx(n_part[1:], (np.diff(wc_time) - np.diff(wtime_advance)))
    # plt.xlabel('#particles')
    # plt.ylabel('time')
    # plt.title('NOT push time')

    plt.figure(5)
    plt.loglog(n_part[1:], np.diff(wtime_advance))
    # plt.loglog(n_part[-1], np.diff(wtime_advance)[-1] /  n_part[-1], 'ro')
    plt.xlabel('#particles')
    plt.ylabel('time')
    plt.title('Push time')

    plt.figure(6)
    plt.loglog(n_part[1:], np.diff(wtime_advance) / n_part[1:])
    # plt.loglog(n_part[-1], np.diff(wtime_advance)[-1] /  n_part[-1], 'ro')
    plt.xlabel('#particles')
    plt.ylabel('time')
    plt.title('Push time per particle')

    plt.figure(7)
    plt.semilogx(n_part, wpercent)
    # plt.semilogx(n_part[1:], 100 * np.diff(wtime_advance) / np.diff(wc_time))
    plt.xlabel('#particles')
    plt.title('percentage spent on particles')

    # plt.figure(8)
    # plt.plot(it, wpercent)
    # plt.plot(it[1:], 100 * np.diff(wtime_advance) / np.diff(wc_time))
    # plt.xlabel('#iterations')
    # plt.title('percentage spent on particles')
    #
    # x_new = np.linspace(1, np.size(it[1:]), 75)
    # a_BSpline = make_interp_spline(it[1:], np.diff(wtime_advance))
    # wtime_advance_diff_smooth = a_BSpline(x_new)
    # plt.plot(x_new, wtime_advance_diff_smooth)
    #
    # x_new = np.linspace(1, np.size(it[1:]), 75)
    # a_BSpline = make_interp_spline(it[1:], np.diff(wc_time))
    # wc_time_diff_smooth = a_BSpline(x_new)
    # plt.plot(x_new, wc_time_diff_smooth)

    # plt.figure(9)
    # plt.plot(it[1:], np.diff(wtime_advance))
    # plt.plot(it[1:], np.diff(wc_time))
    # plt.xlabel('#iterations')
    # plt.title('percentage spent on particles')

    # plt.figure(10)
    # plt.plot(x_new, 100 * wtime_advance_diff_smooth / wc_time_diff_smooth)
    # plt.xlabel('#iterations')
    # plt.title('percentage spent on particles')

plt.show()
