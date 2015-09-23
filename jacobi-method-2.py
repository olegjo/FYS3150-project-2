import numpy as np


def max_off_diag(A):
	''' A is a symmetric matrix. Need only check
		the upper (or lower) half.
	'''
	A_max = 0.
	n = np.shape(A)[0]
	for i in np.arange(n):
		for j in np.arange(i+1,n):
			if abs(A[i,j]) > A_max:
				A_max = abs(A[i, j])
				global l
				global k
				l, k = i, j
	return A_max, l, k

eps = 1e-8


A = np.array([[7., 0., 1.], [0., 7., 0.], [1., 0., 7.]])
n = 3
max_off, l, k = max_off_diag(A)
counter = 0
while max_off**2 > eps:
	tau = (A[l,l] - A[k,k])/(2.*A[k,l])
	tau_root = np.sqrt(1 + tau**2)
	# t = min([-tau + tau_root, -tau - tau_root])
	t = -tau - tau_root
	c = 1./np.sqrt(1 + t**2)
	s = t*c

	for i in np.arange(n):
		if i != k and i != l:
			A_ik = A[i,k]
			A[i,k] = A_ik*c - A[i,l]*s
			A[k,i] = A[i,k]
			A[i,l] = A[i,l]*c + A_ik*s
			A[l,i] = A[i,l]

	A_kk = A[k,k]
	A[k,k] = A_kk*c*c - 2*A[k,l]*c*s + A[l,l]*s*s
	A[l,l] = A[l,l]*c*c + 2*A[k,l]*c*s + A_kk*s*s
	A[k,l] = 0
	A[l,k] = 0
	max_off, l, k = max_off_diag(A)
	counter += 1
print A
print counter