import sys
import numpy as np
import matplotlib.pyplot as plt

def sort(x, A):
	''' sort the vector x ascending, and the matrix A 
	accordingly. '''
	n = len(x)
	A_new = np.zeros((n,n))
	x_new = np.zeros(n)
	j = 0
	for i in np.argsort(x):
		x_new[j] = x[i]
		A_new[:,j] = A[:,i]
		j += 1
	return x_new, A_new

lambdas = []
f = open(sys.argv[1])
for line in f:
	lambdas.append(float(line))
lambdas = np.array(lambdas)
f.close()

f = open(sys.argv[2])
n = len(lambdas)
u = np.zeros((n,n))
i = 0
for line in f:
	numbers = line.split()
	for j in range(n):
		u[i, j] = numbers[j]
	i += 1

lambdas, u = sort(lambdas, u)
for j in range(len(u[0,:])):
	s = 0
	for i in range(len(u[:,j])):
		s += u[i,j]**2
	print s
plt.plot(np.array(u[:, 0])**2)
plt.plot(np.array(u[:, 1])**2)
plt.plot(np.array(u[:, 2])**2)
plt.show()