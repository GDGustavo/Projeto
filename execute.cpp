#include "param.h"
#include "iterd.h"
#include "iterd_r.h"
#include "iterN.h"
#include "iterN_r.h"
#include "mel.h"
#include "mel_2.h"
#include "proj.h"
#include <time.h> 

int N_max = 5;								// Maximum number of iterations;
static int **dimen_L_p;							// Number of basis on the Left (Q,S)_{N-1} sector.
static int **dimen_L;							// Number of basis on the Left (Q,S)_N sector.
static int **dimen_R_p;							//  //   	        Right  //  	.   
static int **dimen_R;							//  //                    Right  //  	.

int main(){
time_t beg_time = time(NULL);


std::cout << "Starting the NRG code..." << std::endl<< std::endl;
std::cout << std::endl;	

read_all_params();    								// Read the parameters to start the code.
double W1 = W_1();								// What is the value of W ? 
double W2 = W_1();
std::cout << std::endl;	

int nq = 2*(2)+1; 								// Total number of possibles q in N=0. 
int ns = 2+1; 									// Total number of possibles dsin N=0.

dimen_L_p = new int*[nq];
for(int i=0; i< nq; i++){
	dimen_L_p[i] = new int[ns];
	for(int j=0; j< ns; j++){
		dimen_L_p[i][j] = 0;
	}
}



iter0_l(W1, dimen_L_p); 								// Initial condictions for "Left Side" N = 0;
iter0_r(W2);									// Initial condictions for "Right Side" N = 0;

dimen_R_p = new int*[nq];
for(int i=0; i< nq; i++){
	dimen_R_p[i] = new int[ns];
	for(int j=0; j< ns; j++){
		dimen_R_p[i][j] = dimen_L_p[i][j];
	}
}

//  Matrix |Q' dS' r'|f_N^+|Q dS r|: N = 0;
//  Eigenvalues: |Q dS r|H_N|Q dS r|: N = 0;
//  EigenVectors: |Q dS r||Q dS p|: N = 0;
//  Dim(0); 

std::cout << "Done! Iteraction 0 complete!" << std::endl<< std::endl<< std::endl;

// 						ITERATIVE SOLUTION N > 0
for(int ni = 1; ni <= N_max; ni ++){   
	nq = 2*(ni+2)+1; 								// Total number of possibles q. 
	ns = ni + 2 +1; 								// Total number of possibles ds.

	dimen_L = new int*[nq]; 							// Matrix to save the basis in the Left sector.
	for(int i=0; i< nq; i++){
		dimen_L[i] = new int[ns];
		for(int j=0; j< ns; j++){
			dimen_L[i][j] = 0;
		}
	}

	dimen_R = new int*[nq]; 							// Matrix to save the basis in the right sector.
	for(int i=0; i< nq; i++){
		dimen_R[i] = new int[ns];
		for(int j=0; j< ns; j++){
			dimen_R[i][j] = 0;
		}
	}

	iterN_L(ni, dimen_L_p, dimen_L);						// Solve the 'Left Side' hamiltonian for N=ni
//	iterN_R(ni, dimen_R_p, dimen_R);						// Solve the 'Right Side' hamiltonian for N=ni

	nq = 2*(ni-1+2)+1; 							// Delete the dimen_L_p matrix. 
	for (int q = 0; q < nq; q++){
		delete[] dimen_L_p[q];
	}
	delete[] dimen_L_p;
	dimen_L_p = NULL;
	dimen_L_p = dimen_L;							// Left: New imput of number of basis
										 
	for (int q = 0; q < nq; q++){						// Delete the dimen_R_p matrix.
		delete[] dimen_R_p[q];
	}
	delete[] dimen_R_p;
	dimen_R_p = NULL;
	dimen_R_p = dimen_R;							// Right: New imput of number of basis

	std::cout << "Done! Iteraction " << ni <<  " complete!" << std::endl<< std::endl<< std::endl<< std::endl;
}

time_t end_time = time(NULL);

float time_ = difftime(end_time, beg_time) /60;    					// Converte in seconds to minutes

std::cout << "The code is executed in " << time_  <<  " minutes." << std::endl;

return 0;
}


// Até a iteração N = 3, todas as autoenergias foram corrigidas. 
// 18/01/2022 até o N = 5, tempo de execução:  4.40 min. (Cada setor) 
// Projeções conferidas até N= 4. Mesmos valores de W1 e W2.
// Energia de corte ultravioleta E_uv inserido no Bra "Left".


