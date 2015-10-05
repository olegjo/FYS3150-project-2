# FYS3150 Project 2
This is a numerical project in the course FYS3150 at Univeristy of Oslo consentrating on the implementation of
Jacobi method/Given's rotation algorithm for finding eigenvalues.

The specific case investigated is the solution of the Schrödinger equation for one and two electrons 
in a harmonic oscillator potential.

## Dependencies
C++11, Python with matplotlib and python

## Benchmarks
A good benckmark calculation is to run 
```
./task_a 220 5
```
from the the build*/task_a/ folder, followed by
```
python readfile.py oneElectron
```
from the Python/ folder.

The first should give output in the terminal 
```
Number of iterations: 65339
```
and the latter results in a plot with the eigenvalues 2.999, 6.999 and 10.998.

## Unit tests
Run 
```
./unit_tests
```
int the build*/unit_tests folder
