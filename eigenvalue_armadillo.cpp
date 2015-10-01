#include "iostream"
#include "armadillo"

using namespace std;
//using namespace arma;

int main(int argc, char const *argv[])
{

	int n = atoi(argv[1]);

	arma::mat A = arma::zeros<arma::mat>(n,n);
	arma::vec V = arma::zeros<arma::vec>(n);

	double rho_min = 0;
	double rho_max = atof(argv[2]);
	double h = (rho_max - rho_min)/n;
	double h_squared = h*h;
	double e_elements = -1.0/h_squared;
	double d_elements = 2.0/h_squared;
	V(0) = rho_min;
	A(0,0) = d_elements + V[0];
	A(0,1) = e_elements;

	for (int i = 1; i < n - 1; i++) {
		V(i) = rho_min + i*h;
		A(i,i) = d_elements + V[i];
		A(i,i+1) = e_elements;
		A(i,i-1) = e_elements;
	}
	V(n-1) = V(n-2) + h;
	A(n-1,n-1) = d_elements + V(n-1);
	A(n-1,n-2) = e_elements;

	arma::vec eigval;
	arma::mat eigvec;

	eig_sys(eigval, eigvec, A, 5);
	//eigval.print();
	return 0;
}
