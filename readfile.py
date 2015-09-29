import sys
import numpy as np

lambdas = []

f = open(sys.argv[1])
for line in f:
	lambdas.append(float(line))

lambdas = np.array(lambdas)
lambdas = np.sort(lambdas)
print lambdas[:3]