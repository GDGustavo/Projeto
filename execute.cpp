#include "param.h"
#include "iterd.h"
#include "iterd_r.h"
#include "iterN.h"
#include "iterN_r.h"
#include "mel.h"
#include "mel_2.h"
#include <new>
#include "proj.h"

int N_max = 4;								// Maximum number of iterations;
static int **dimen_p_;							// Matrix: How many bases have in (Q,S)_{N-1} sector.  
static int **dimen_;							// Matrix: How many bases have in (Q,S)_N sector.

int main(){

std::cout << "Starting the NRG code..." << std::endl<< std::endl;
std::cout << std::endl;	

read_all_params();    								// Read the parameters to start the code.
double W1 = W_1();								// What is the value of W ? 
double W2 = W_1();
std::cout << std::endl;	

int nq = 2*(2)+1; 								// Total number of possibles q in N=0. 
int ns = 2+1; 									// Total number of possibles dsin N=0.

dimen_p_ = new int*[nq];
for(int i=0; i< nq; i++){
	dimen_p_[i] = new int[ns];
	for(int j=0; j< ns; j++){
		dimen_p_[i][j] = 0;
	}
}

iter0_l(W1, dimen_p_); 								// Initial conditions for "Left Side" N = 0;
//iter0_r(W2);									// Initial condictions for "Right Side" N = 0;

//  Matrix |Q' dS' r'|f_N^+|Q dS r|: N = 0;
//  Eigenvalues: |Q dS r|H_N|Q dS r|: N = 0;
//  EigenVectors: |Q dS r||Q dS p|: N = 0;
//  Dim(0); 

std::cout << "Done! Iteraction 0 complete!" << std::endl<< std::endl<< std::endl;

// 						ITERATIVE SOLUTION N > 0
for(int ni = 1; ni <= N_max; ni ++){   
	nq = 2*(ni+2)+1; 								// Total number of possibles q. 
	ns = ni + 2 +1; 								// Total number of possibles ds.

	dimen_ = new int*[nq]; 							// Matrix to save how many bases have for N=ni.
	for(int i=0; i< nq; i++){
		dimen_[i] = new int[ns];
		for(int j=0; j< ns; j++){
			dimen_[i][j] = 0;
		}
	}

	iterN_L(ni, dimen_p_, dimen_);						// Solve the 'Left Side' hamiltonian for N=ni
	//iterN_R(ni, dimen_p_, dimen_);						// Solve the 'Right Side' hamiltonian for N=ni

	nq = 2*(ni-1+2)+1; 							// Delete the dimen_p_ matrix. 
	for (int q = 0; q < nq; q++){
		delete[] dimen_p_[q];
	}
	delete[] dimen_p_;
	dimen_p_ = NULL;

	dimen_p_ = dimen_;							// New imput: number of bases

	std::cout << "Done! Iteraction " << ni <<  " complete!" << std::endl<< std::endl<< std::endl<< std::endl;
}

	nq = 2*(N_max+2)+1;
	eigen_delete(N_max,dimen_p_); 						// Delete all elements saved in eigen matrix
	mel_nw_delete(N_max,dimen_p_); 						// Delete all elements saved in mel_nw_
	mel_ne_delete(N_max,dimen_p_); 						// Delete all elements saved in mel_ne_
	//projection_delete(N_max,dimen_p_);						// Delete all elements saved in projection_

	for (int q = 0; q < nq; q++){
		delete[] dimen_[q];
	}
	delete[] dimen_;

return 0;
}



// Até a iteração N = 2, todas as autoenergias foram corrigidas. 
// 05/01/2022 até o N =  4, 2:40 min. 
// Projeções conferidas até N= 3. Mesmos valores de W1 e W2. 
//  


