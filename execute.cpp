#include "param.h"
#include "iterd.h"
#include "iterd_r.h"
#include "iterN.h"
#include "general.h"
#include "mel.h"
#include "mel_2.h"
#include <new>
#include "proj.h"

int N_max = 3;									// Maximum number of iterations;
static int **dimen_;								// Matrix: How many bases have in (Q,S) sector.  
static int **dimen_out_;								// Matrix: How many bases have in (Q,S) sector.

int main(){

std::cout << "Starting the NRG code..." << std::endl<< std::endl;
std::cout << std::endl;	

read_all_params();    								// Read the parameters to start the code.
double W1 = W_1();								// What is the value of W ? 
double W2 = W_2();
std::cout << std::endl;	

int nq = 2*q_max(0)+1; 								// Total number of possibles q'. 
int ns = ds_max(0) +1; 								// Total number of possibles ds'.

dimen_ = new int*[nq]; 
for(int i=0; i< nq; i++){
	dimen_[i] = new int[ns];
	for(int j=0; j< ns; j++){
		dimen_[i][j] = 0;
	}
}

iter0(W1, dimen_); 								// Initial conditions for N = 0;
iter0_r(W2);
//projection_start_0(0,dimen_);

//  Matrix |Q' dS' r'|f_N^+|Q dS r|: N = 0;
//  Eigenvalues: |Q dS r|H_N|Q dS r|: N = 0;
//  EigenVectors: |Q dS r||Q dS p|: N = 0;
//  Dim(0); 

std::cout << "Done! Iteraction 0 complete!" << std::endl<< std::endl;

// 							ITERATIVE SOLUTION 
for(int ni = 1; ni <= N_max; ni ++){   
	nq = 2*q_max(ni)+1; 							// Total number of possibles q'. 
	ns = ds_max(ni) +1; 							// Total number of possibles ds'.
	dimen_out_ = new int*[nq]; 						// Matrix to save how many bases have for ni.
	for(int i=0; i< nq; i++){
		dimen_out_[i] = new int[ns];
		for(int j=0; j< ns; j++){
			dimen_out_[i][j] = 0;
		}
	}

	iterN(ni, dimen_, dimen_out_);						// Solving the hamiltonian for N=ni
	nq = 2*q_max(ni-1)+1; 
	for (int q = 0; q < nq; q++){
		delete[] dimen_[q];
	}
	delete[] dimen_;
	dimen_ = NULL;

	dimen_ = dimen_out_;							// New imput number of bases
}

	nq = 2*q_max(N_max)+1; 
	eigen_delete(N_max,dimen_);  // Delete all elements saved in eigen matrix
	mel_nw_delete(N_max,dimen_); // Delete all elements saved in mel_nw_
	mel_ne_delete(N_max,dimen_); // Delete all elements saved in mel_ne_

	for (int q = 0; q < nq; q++){
		delete[] dimen_out_[q];
	}
	delete[] dimen_out_;
	dimen_out_ = NULL;


return 0;
}
