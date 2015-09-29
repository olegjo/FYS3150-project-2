#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "jacobi-method.h"
using namespace std;
int main(int argc, char const *argv[])
{
	double ** A;
	double * V;

	int n_step = atoi(argv[1]);
	double rho_min = 0;
	double rho_max = atof(argv[2]);
	double h = (rho_max - rho_min)/n_step;


	A = new double* [n_step]; V = new double [n_step];
	for (int i = 0; i < n_step; i++){
		A[i] = new double [n_step];
		double temp = rho_min + i*h;
		V[i] = temp*temp;
	}


	// creating the matrix A
	// First: first and last rows
	double h_squared = h*h;
	double e_elements = -1.0/h_squared;
	double d_elements = 2.0/h_squared;
	A[0][0] = d_elements + V[0];
	A[0][1] = e_elements;
	A[n_step - 1][n_step - 1] = d_elements + V[n_step - 1];
	A[n_step - 1][n_step - 2] = e_elements;
	// Next: all other rows
	for (int i = 1; i < n_step - 1; i++) {
		A[i][i] = d_elements + V[i];
		A[i][i+1] = e_elements;
		A[i][i-1] = e_elements;
	}


	jacobi_method(A, n_step);


	// write the data to file
	ofstream targetfile;
	targetfile.open("outfile_mine.txt");
	for (int i = 0; i < n_step; i++){
		targetfile << setw(15) << setprecision(8) << A[i][i] << endl;
	}

	targetfile.close();

	delete[] A;
	delete[] V;


	return 0;
}

void jacobi_method(double **A, int n)
{
	int iterations = 0;
	int k, l;
	double eps = 1.0e-8;
	double max_off = max_off_diag(A, &k, &l, n);
	double max_number_iterations = (double) n * (double) n * (double) n; 
	while (max_off*max_off > eps && (double) iterations < max_number_iterations) {
		rotate(A, k, l, n);
		max_off = max_off_diag(A, &k, &l, n);
		iterations++;
	}
	std::cout << "Number of iterations: " << iterations << std::endl;
}

void rotate(double **A, int k, int l, int n)
{
	double tau;
	double s, c, t;
	double A_ik, A_kk;

	if (A[k][l] != 0) {
		tau = (A[l][l] - A[k][k])/(2.0*A[k][l]);
		t = -tau - sqrt(1 + tau*tau);
		c = 1.0/sqrt(1 + t*t);
		s = t*c;
	} else {
		c = 1.0;
		s = 0.0;
	}


	for (int i = 0; i < n; i++) {
		if (i != k && i != l) {
			A_ik = A[i][k];
			A[i][k] = A_ik*c - A[i][l]*s;
			A[k][i] = A[i][k];
			A[i][l] = A[i][l]*c + A_ik*s;
			A[l][i] = A[i][l];
		}
	}

	A_kk = A[k][k];
	A[k][k] = A_kk*c*c - 2*A[k][l]*c*s + A[l][l]*s*s;
	A[l][l] = A[l][l]*c*c + 2*A[k][l]*c*s + A_kk*s*s;
	A[k][l] = 0;
	A[l][k] = 0;
	return;
}


double max_off_diag(double **A, int *k, int *l, int n)
{
	double A_max = 0.0;
	for (int i = 0; i < n; i++){
		for (int j = i+1; j < n; j++){			
			if (fabs(A[i][j]) > A_max){
				A_max = fabs(A[i][j]);
				*l = i; *k = j;				
			}
		}
	}
	return A_max;
}