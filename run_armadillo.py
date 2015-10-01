import os
import sys
import time

program = sys.argv[1]
os.system('g++ -c -Wall %s -larmadillo -lblas -llapack' % program)
program = program[:-4]
os.system('g++ -o %s %s.o -larmadillo -lblas -llapack' % (program, program))
try:
	os.system('./%s %s' % (program, sys.argv[2]))
except IndexError:
	os.system('./%s' % (program))
