import numpy as np
import sys
import time

start = time.time()
n = int(sys.argv[1])
rho_max = float(sys.argv[2])

A = np.zeros([n,n])
V = np.zeros(n)



rho_min = 0
h = (rho_max - rho_min)/float(n)
h_squared = h*h
e_elements = -1./h_squared
d_elements = 2./h_squared

for i in range(n):
	temp = rho_min + (i+1)*h
	V[i] = temp*temp


A[0,0] = d_elements + V[0]
A[0,1] = e_elements
A[-1,-1] = d_elements + V[-1]
A[-1,-2] = e_elements

for i in range(1,n-1):
	A[i,i] = d_elements + V[i]
	A[i,i+1] = e_elements
	A[i,i-1] = e_elements

eigenvals, eigenvecs = np.linalg.eig(A)
print np.sort(eigenvals)[:3]
print 
print 'n =', n
print 'Computation time: %.3f s' % (time.time() - start)
#print A
