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
h = (rho_max - rho_min)/n
rho_i = np.array([rho_min + i*h for i in range(0, n)])

#sort them in ascending order of the eigenvalues
lambdas, u = sort(lambdas, u)
print 'rho_max = %.2f, n = %i' % (rho_max, n)
print lambdas[:3]


plt.subplot(3,1,1)
plt.plot(rho_i, np.array(u[:, 0])**2*100, label='$\lambda_0=$%.3f' % lambdas[0])
plt.tick_params(axis='x', labelbottom='off')
plt.legend(fontsize=11)
plt.subplot(3,1,2)
plt.plot(rho_i, np.array(u[:, 1])**2*100, label='$\lambda_1=$%.3f' % lambdas[1])
plt.tick_params(axis='x', labelbottom='off')
plt.ylabel('$P(\\rho)=|\psi|^2$ [%]')
plt.legend(fontsize=11)
plt.subplot(3,1,3)
plt.plot(rho_i, np.array(u[:, 2])**2*100, label='$\lambda_2=$%.3f' % lambdas[2])
plt.legend(fontsize=11)
plt.xlabel("$\\rho$ [-]")
filename = 'plot-'+problem+'.pdf'
plt.savefig('plot-oneElectron.pdf')
os.system('mv plot-oneElectron.pdf ../Report/')




