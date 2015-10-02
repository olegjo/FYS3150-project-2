#include "unit_tests.h"
#include "../jacobi.h"
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;


int TEST_JACOBI()
{
    cout << endl << "*****************************" << endl;
    cout << "TESTING FUNCTION: jacobi_method" << endl;
    cout << "CONDITION: for a 2x2 matrix, the number of iterations should be 1" << endl;
    bool allPassed = true;
    double **A;
    A = new double * [2];
    A[0] = new double [2];
    A[1] = new double [2];
    A[0][0] = 3; A[0][1] = 4;
    A[1][0] = 4; A[1][1] = 3;
    int iterations = jacobi_method(A, 2);
    if (iterations == 1) {
        cout << "PASSED" << endl;
    }
    else {
        cout << "FAILED" << endl;
        allPassed = false;
    }

    cout << "NEW CONDITION: Does it make the correct eigenvalues and are the non-diagonal elements zero?" << endl;
    double eps = 1e-8;
    if ((fabs(A[0][0]+1) < eps && fabs(A[1][1]-7) < eps) || (fabs(A[0][0]-7) < eps && fabs(A[1][1]+1) < eps)) {
        cout << "PASSED" << endl;

    }
    else {
        cout << "FAILED" << endl;
        allPassed = false;
    }

    if (allPassed) {
        return 1;
    }
    return 0;

}

int TEST_MAX_OFF_DIAG()
{
    cout << endl << "*****************************" << endl;
    cout << "TESTING FUNCTION: max_off_diag" << endl;
    cout << "CONDITION: checking if it returns the largest element and its correct position" << endl;
    bool allPassed = true;
    double **A;
    A = new double * [3];
    A[0] = new double [3];
    A[1] = new double [3];
    A[2] = new double [3];
    A[0][0] = 3; A[0][1] = 4; A[0][2] = 300;
    A[1][0] = 4; A[1][1] = 3; A[1][2] = -250;
    A[2][0]=300; A[2][1] = -250; A[2][2] = 3;
    int k, l;
    double max_off = max_off_diag(A, &k, &l, 3);
    if (k==2 && l==0 && max_off==300) {
        cout << "PASSED" << endl;
    }
    else {
        cout << "FAILED" << endl;
        allPassed = false;
    }

    cout << endl << "NEW CONDITION: does it also work if the largest element is negative?" << endl;
    A[1][2] = -350;
    max_off = max_off_diag(A, &k, &l, 3);
    if (k==2 && l==1 && max_off==350) {
        cout << "PASSED" << endl;
    }
    else {
        cout << "FAILED" << endl;
        allPassed = false;
    }
    delete[] A;
    if (allPassed) {
        return 1;
    }
    else {
        return 0;
    }
}

int TEST_ROTATE_EIGENVECTOR_ORTHOGONAL()
{
    cout << endl << "*****************************" << endl;
    cout << "TESTING FUNCTION: rotate" << endl;
    cout << "CONDITION: must return orthogonal eigenvectors" << endl;
    bool passedAll = true;
    double **A;
    A = new double * [3];
    A[0] = new double [3];
    A[1] = new double [3];
    A[2] = new double [3];
    A[0][0] = 4; A[0][1] = 2; A[0][2] = 4;
    A[1][0] = 2; A[1][1] = 1; A[1][2] = 2;
    A[2][0] = 4; A[2][1] = 2; A[2][2] = 4;
    double **R;
    R = new double * [3];
    R[0] = new double [3];
    R[1] = new double [3];
    R[2] = new double [3];
    R[0][0] = 1; R[0][1] = 0; R[0][2] = 0;
    R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
    R[2][0] = 0; R[2][1] = 0; R[2][2] = 1;

    int k, l;
    double eps = 1.0e-8;
    int iterations = 0;
    double max_number_iterations = (double) 3 * (double) 3 * (double) 3;
    double max_off = max_off_diag(A, &k, &l, 3);
    while (max_off*max_off > eps && (double) iterations < max_number_iterations) {
        rotate(A, R, k, l, 3);
        max_off = max_off_diag(A, &k, &l, 3);
        double orth12 = R[0][0]*R[0][1] + R[1][0]*R[1][1] + R[2][0]*R[2][1];
        double orth13 = R[0][0]*R[0][2] + R[1][0]*R[1][2] + R[2][0]*R[2][2];
        double orth23 = R[0][1]*R[0][2] + R[1][1]*R[1][2] + R[2][1]*R[2][2];
        double eps = 1e-8;
        if (fabs(orth12) <! eps || fabs(orth13) <! eps || fabs(orth23) <! eps) {
            cout << "FAILED" << endl;
            passedAll = false;
        }
        iterations++;
    }


    if (passedAll) {
        cout << "PASSED" << endl;
        return 1;
    }
    else {
        cout << "FAILED" << endl;
        return 0;
    }
}

void RUN_ALL_TESTS()
{
    int numberOfTests = 3;
    int passed = 0;
    passed += TEST_JACOBI();
    passed += TEST_MAX_OFF_DIAG();
    passed += TEST_ROTATE_EIGENVECTOR_ORTHOGONAL();

    cout << endl << "PASSED " << passed << " of " << numberOfTests << " tests" << endl;
}



