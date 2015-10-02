#include "potentials.h"

namespace Potentials {
    void oneElectron(double **A, int n, double rho_min, double h) {
        for (int i = 0; i < n; i++) {
            double rho_i = rho_min + (i+1)*h;
            A[i][i] += rho_i*rho_i;
        }
    }

    void twoElectrons(double **A, int n, double rho_min, double h, double omega_r) {
        for (int i = 0; i < n; i++) {
            double rho_i = rho_min + (i+1)*h;
            A[i][i] += omega_r*omega_r*rho_i*rho_i + 1.0/rho_i;
        }
    }
}


