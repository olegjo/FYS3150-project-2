#ifndef JACOBI_H
#define JACOBI_H


class Jacobi
{
public:
    Jacobi();
    void jacobi_method(double **A, int n);
    void jacobi_method(double **A, double **R, int n);

private:
    void rotate(double **A, int k, int l, int n);
    void rotate(double **A, double **R, int k, int l, int n);
    double max_off_diag(double **A, int *k, int *l, int n);
};

#endif // JACOBI_H
