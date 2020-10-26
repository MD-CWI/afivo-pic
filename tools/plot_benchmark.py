import numpy as np
import matplotlib.pyplot as plt
import os

file = "/home/ddb/codes/afivo-pic/programs/standard_2d/output/test_diagnostics"

dat = np.loadtxt(file)

it            = dat[:, 0]
n_part        = dat[:, 1]
wc_time       = dat[:, 2]
wtime_advance = dat[:, 3] # Watch out: these are CUMULATIVE
wtime_run     = dat[:, 4] # Watch out: these are CUMULATIVE
wpercent      = dat[:, 5]

plt.figure()
plt.semilogy(it, n_part)
plt.xlabel('#iterations')
plt.ylabel('#particles')
plt.title('Particles over time')

plt.figure()
plt.semilogy(it[1:], np.diff(wc_time))
plt.xlabel('#iterations')
plt.ylabel('time')
plt.title('Calculation time of one iteration')

plt.figure()
plt.loglog(n_part[1:], np.diff(wc_time))
plt.xlabel('#particles')
plt.ylabel('time')
plt.title('Calculation time of one iteration')

plt.figure()
plt.loglog(n_part[1:], np.diff(wtime_advance))
# plt.loglog(n_part[-1], np.diff(wtime_advance)[-1] /  n_part[-1], 'ro')
plt.xlabel('#particles')
plt.ylabel('time')
plt.title('Push time')

plt.figure()
plt.loglog(n_part[1:], np.diff(wtime_advance) / n_part[1:])
# plt.loglog(n_part[-1], np.diff(wtime_advance)[-1] /  n_part[-1], 'ro')
plt.xlabel('#particles')
plt.ylabel('time')
plt.title('Push time per particle')

plt.figure()
plt.semilogx(n_part, wpercent)
# plt.semilogx(n_part[1:], 100 * np.diff(wtime_advance) / np.diff(wc_time))
plt.xlabel('#particles')
plt.title('percentage spent on particles')

plt.figure()
plt.plot(it, wpercent)
plt.xlabel('#iterations')
plt.title('percentage spent on particles')

plt.show()
