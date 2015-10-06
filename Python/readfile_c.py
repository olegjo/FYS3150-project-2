import sys
import numpy as np
import matplotlib.pyplot as plt
import os

def sort(x, A):
	''' sort the vector x ascending, and the matrix A 
	accordingly, so that if element i in x moves to 
	place j, the i-th columb in A also moves to place j. 
	'''
	n = len(x)
	A_new = np.zeros((n,n))
	x_new = np.zeros(n)
	j = 0
	for i in np.argsort(x):
		x_new[j] = x[i]
		A_new[:,j] = A[:,i]
		j += 1
	return x_new, A_new

problem = sys.argv[1]
location = '../results/'

# extracting the eigenvalues
lambdas = []
f = open(location+'results_eigenValues_'+problem+'.txt')
for line in f:
	lambdas.append(float(line))
lambdas = np.array(lambdas)
f.close()

# extracting the eigenvectors
f = open(location+'results_eigenVectors_'+problem+'.txt')
n = len(lambdas)
u = np.zeros((n,n))
i = 0
for line in f:
	numbers = line.split()
	for j in range(n):
		u[i, j] = numbers[j]
	i += 1

# getting the parameters
f = open(location+'parameters_'+problem+'.txt')
lines = []
for line in f:
	lines.append(line)
rho_min = float(lines[0].split()[2])
rho_max = float(lines[1].split()[2])
n = int(lines[2].split()[2])
omega_r = float(lines[3].split()[2])
h = (rho_max - rho_min)/n
rho_i = np.array([rho_min + i*h for i in range(0, n)])

#sort them in ascending order of the eigenvalues
lambdas, u = sort(lambdas, u)
print 'rho_max = %.2f, n = %i, omega_r = %.3f' % (rho_max, n, omega_r)
print lambdas[:3]


plt.plot(rho_i, np.array(u[:, 0])**2*100, label='$\lambda_0=$%.3f \n $\omega_r=$%.2f' % (lambdas[0], omega_r))
plt.legend(fontsize=17)
plt.xlabel("$\\rho$ [-]", fontsize=19)
plt.ylabel("$P(\\rho) = |\psi|^2$ [%]", fontsize=19)
filename = 'plot-'+problem+'.pdf'
plt.savefig(filename)
os.system('mv %s ../Report/' % (filename))







