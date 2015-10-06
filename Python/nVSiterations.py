import matplotlib.pyplot as plt
import numpy as np
import os

n = np.array([100, 150, 200, 250, 300, 350, 400], np.float)
iterations = np.array([13185, 30201, 53770, 84299, 122192, 167112, 219750], np.float)



alpha = np.average(iterations/n**2)
plt.plot(n, iterations, '--o', label='Experimental')
plt.plot(n, alpha*n**2, '--o', label='$\\alpha n^2$')
plt.legend(fontsize=11)
plt.xlabel('$n$ - number of grid points')
plt.ylabel('Iterations needed by jacobi method')
plt.savefig('plot-grid-pointsVSiterations.pdf')
os.system('mv plot-grid-pointsVSiterations.pdf ../Report/')
