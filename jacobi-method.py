import numpy as np

def maxoffdiag(A, k, l):
	'return the max off-diagonal element of A and its position'
	Amax = 0.
	n = np.shape(A)[0]
	for i in np.arange(0, n):
		for j in np.arange(0, n):
			if i != j and abs(A[i,j]) > Amax:
				Amax = A[i,j]
				l, k = i, j
	return k, l, Amax




eps = 1e-8


# eigenvectors for A: (1, 0, -1), (0, 1, 0), (1, 0, 1)
# eigenvalues respec: 1, 2, 3
A = np.array([[2., 0., 1.], [0., 2., 0.], [1., 0., 2.]])



# choosing a tolerance
eps = 1e-8

k,l, max_offdiag = maxoffdiag(A, k=None, l=None)
while max_offdiag**2 > eps:
	tau = (A[l,l] - A[k,k])/(2.*A[k,l])
	if tau > 0:
		tau_sqrt = np.sqrt(tau**2 + 1.)
		t = np.min([-tau + tau_sqrt, -tau - tau_sqrt])
		c = 1./np.sqrt(1 + t**2)
		s = t*c
	else:
		c, s = 1., 0.

	
	# changing the elements with index permutations of k and l
	a_ll = A[l,l]
	a_kk = A[k,k]
	A[k,k] = a_kk*c**2 - 2*A[k,l]*c*s + a_ll*s**2
	A[l,l] = a_ll*c**2 + 2*A[k,l]*c*s + a_kk*s**2
	A[k,l] = 0.
	A[l,k] = 0.
	# changing the rest
	for i in np.arange(0,len(A[:,0])):
		if i != l and i != k:
			A_ik = A[i,k]
			A_il = A[i,l]
			A[i,k] = A_ik*c - A_il*s
			A[k,i] = A[i,k]
			A[i,l] = A_il*c + A_ik*s
			A[l,i] = A[i,l]


	# get ready for next iteration
	k, l, max_offdiag = maxoffdiag(A, k, l)

print A


