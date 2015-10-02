#include <iostream>
#include <fstream>
#include <iomanip>
#include <../jacobi.h>
#include <../potentials.h>

using namespace std;


int main(int argc, char* argv[])
{
    // setting up the matrix A and R
    double ** A;
    double ** R;
    int n = atoi(argv[1]);
    double rho_min = 0;
    double rho_max = atof(argv[2]);

    A = new double* [n];
    R = new double* [n];
    for (int i = 0; i < n; i++){
        A[i] = new double [n];
        R[i] = new double [n];
    }

    double h = (rho_min - rho_max)/n;
    double h_squared = h*h;
    double d_elements = 2.0/h_squared;
    double e_elements = -1.0/h_squared;
    A[0][0] = d_elements;
    A[0][1] = e_elements;
    A[n-1][n-1] = d_elements;
    A[n-1][n-2] = e_elements;
    for (int i = 1; i < n; i++) {
        A[i][i] = d_elements;
        A[i][i+1] = e_elements;
        A[i][i-1] = e_elements;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i==j) {
                R[i][j] = 1.0;
            }
        else {
                R[i][j] = 0.0;
            }
        }
    }

    // add the potentials to the diagonal
    Potentials::oneElectron(A, n, rho_min, h);

    jacobi_method(A, R, n);


//    // write the data to file
//    ofstream targetfile_eigenValues;
//    targetfile_eigenValues.open("results_eigenValues_oneElectron.txt");
//    for (int i = 0; i < n; i++){
//        targetfile_eigenValues << setw(15) << setprecision(8) << A[i][i] << endl;
//    }
//    targetfile_eigenValues.close();




//    ofstream targetfile_eigenVectors;
//    targetfile_eigenVectors.open("results_eigenVectors_oneElectron.txt");
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            targetfile_eigenVectors << setw(15) << setprecision(8) << R[i][j] << " ";
//        }
//        targetfile_eigenVectors << endl;
//    }
//    targetfile_eigenVectors.close();



    delete[] A;
    delete[] R;

    return 0;
}

