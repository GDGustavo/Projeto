#include "giv0.1.h"
#include <cmath>
#include <cstdio>
#include <new>
#include <cstdlib>
#include <iostream>
#include <iomanip>

int main(int argc, char** argv )
{
    double* ham = nullptr;
    double* root = nullptr;
    double** vect = nullptr;

    int nn = 9;
    int dim = nn + 1;
    ham = new double[dim * (dim + 1)/2];
    root = new double[dim + 1];
    vect = new double*[dim];
    for(int lin = 0; lin < dim; lin++){
        vect[lin] = new double[dim];
    }

    for(int n = 0; n < dim * (dim+1)/2; n++){
        ham[n] = 0.0;
    }
    double t_tilde = pow(2.0, nn - 0.5);
    for(int lin = 1; lin < dim; lin++){
        int indx = (lin + 1) * (lin+2) / 2 - 2;
        ham[indx] = t_tilde * pow(2.0, -lin + 0.5);
    }
    double cutoff = 5.0;
    root[0] = cutoff;
    int nroot = givens(dim, 0, ham, root, vect, 1);
    std::cout << std::setprecision(4);
    std::cout << "No. of vectors below cutoff == " << cutoff << " is " << nroot << std::endl;
    std::cout << "Full list of eigenvalues follows: " << std::endl;
    for (int n = 0; n < dim; n++) { std::cout << n << ": " <<  '\t' << root[n] << std::endl; }
    std::cout << "Zero-th elements of calculated eigenvectors " << std::endl;
    for (int n = 0; n < dim; n++) { std::cout << n << ": " <<  '\t' << vect[n][0] << std::endl; }
}