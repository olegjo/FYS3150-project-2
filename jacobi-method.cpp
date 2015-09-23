#include <iostream>
#include <cmath>
#include "jacobi-method.h"

int main(int argc, char const *argv[])
{
	double ** A;
	int n = 3;

	A = new double* [n];
	for (int i = 0; i < n; i++){
		A[i] = new double [n];
	}
	A[0][0] = 2.0;
	A[0][1] = 0.0;
	A[0][2] = 1.0;
	A[1][0] = 0.0;
	A[1][1] = 2.0;
	A[1][2] = 0.0;
	A[2][0] = 1.0;
	A[2][1] = 0.0;
	A[2][2] = 2.0;
	jacobi_method(A, n);

	std::cout << "lambda_1 = " << A[0][0] << std::endl;

	std::cout << "lambda_2 = " << A[1][1] << std::endl;

	std::cout << "lambda_3 = " << A[2][2] << std::endl;
	return 0;
}

void jacobi_method(double **A, int n)
{
	int iterations = 0;
	int k; int l;
	double eps = 1.0e-8;
	double max_off = max_off_diag(A, &l, &k, n);
	while (max_off*max_off > eps) {
		rotate(A, l, k, n);
		max_off = max_off_diag(A, &l, &k, n);
		iterations++;
	}
	std::cout << "Number of iterations: " << iterations << std::endl;

}

void rotate(double **A, int k, int l, int n)
{
	double tau, tau_root;
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
			A[i][k] = A[i][k]*c - A[i][l]*s;
			A[k][i] = A[i][k];
			A[i][l] = A[i][l]*c + A_ik*s;
			A[l][i] = A[i][l];
		}
	}

	A_kk = A[k][k];
	A[k][k] = A_kk*c*c - 2*A[k][l]*c*s + A[l][l]*s*s;
	A[l][l] = A[l][l]*c*c + 2*A[k][l]*c*s + A[k][k]*s*s;
	A[k][l] = 0;
	A[l][k] = 0;
}


double max_off_diag(double **A, int *l, int *k, int n)
{
	double A_max = 0;
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