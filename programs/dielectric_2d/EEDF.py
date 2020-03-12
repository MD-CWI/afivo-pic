import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os

path = os.getcwd() + '/output/EEDF/'
files = os.listdir(path)

str1 = 'sim_0000'
str2 = '_EEDF.txt'

plt.figure()
for ii in [10, 20, 30, 40, 50, 60, 66]:
	grid = np.genfromtxt(path+str1+str(ii)+str2, skip_header=1, usecols=(0))
	values = np.genfromtxt(path+str1+str(ii)+str2, skip_header=1, usecols=(1))
	plt.step(grid,values, label=str(ii))
	plt.yscale('log')

plt.legend()
plt.xlabel('eV')
plt.ylabel('num elec')
plt.show()



