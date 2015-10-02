#include "jacobi.h"
#include <iostream>
#include <cmath>

int jacobi_method(double **A, double **R, int n)
{
    int iterations = 0;
    int k, l;
    double eps = 1.0e-8;
    double max_off = max_off_diag(A, &k, &l, n);
    double max_number_iterations = (double) n * (double) n * (double) n;
    while (max_off*max_off > eps && (double) iterations < max_number_iterations) {
        rotate(A, R, k, l, n);
        max_off = max_off_diag(A, &k, &l, n);
        iterations++;
    }
    return iterations;
}

int jacobi_method(double **A, int n)
{
    double **R;
    R = new double* [n];
    for (int i = 0; i < n; i++){
        R[i] = new double [n];
    }
    int returnValue = jacobi_method(A, R, n);
    delete[] R;
    return returnValue;
}




void rotate(double **A, double **R, int k, int l, int n)
{
    double tau;
    double s, c, t;
    double A_ik, A_kk, r_ik, r_il;

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
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
     }

    A_kk = A[k][k];
    A[k][k] = A_kk*c*c - 2*A[k][l]*c*s + A[l][l]*s*s;
    A[l][l] = A[l][l]*c*c + 2*A[k][l]*c*s + A_kk*s*s;
    A[k][l] = 0;
    A[l][k] = 0;

    return;
}

void rotate(double **A, int k, int l, int n)
{
    double **R;
    R = new double* [n];
    for (int i = 0; i < n; i++){
        R[i] = new double [n];
    }
    rotate(A, R, k, l, n);
    delete[] R;
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
