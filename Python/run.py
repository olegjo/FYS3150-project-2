import os
import sys
n = sys.argv[1]
rho_max = sys.argv[2]
print 'my function: '
os.system('./jacobi-method %s %s' % (str(n), str(rho_max)))
os.system('python readfile.py outfile_mine.txt')
# print '\nnot mine: '
# os.system('./jacobi %s %s' % (str(n), str(rho_max)))
# os.system('python readfile.py outfile_not_mine.txt')
# print ''
# n = 1000: lowest three eigenvalues:
# [  2.7792062   6.6546918  10.548347 ]
# [Finished in 250.0s]