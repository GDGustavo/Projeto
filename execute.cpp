#include "param.h"
#include "iterd.h"
#include "iterd_r.h"
#include "iterN.h"
#include "iterN_r.h"
#include "mel.h" 
#include "mel_2.h"
#include "proj.h"
#include "param.h"
#include <time.h> 

int N_max = fN_max();								// Maximum number of iterations;
static int **dimen_p;								// Number of basis on the (Q,S)_{N-1} sector.
static int **dimen;								// Number of basis on the (Q,S)_N sector.
int main(){
time_t beg_time = time(NULL);

std::cout << "Starting the NRG code..." << std::endl;
std::cout << std::endl;	

read_all_params();    								// Read the parameters to start the code.
double W1 = W_1();								// What is the value of W ? 
double W2 = W_2();
std::cout << std::endl;	

int nq = 2*(2)+1; 								// Total number of possibles q in N=0. 
int ns = 2+1; 									// Total number of possibles dsin N=0.

dimen_p = new int*[nq];
for(int i=0; i< nq; i++){
	dimen_p[i] = new int[ns];
	for(int j=0; j< ns; j++){
		dimen_p[i][j] = 0;
	}
}

iter0_l(W1, dimen_p); 								// Initial condictions for "Left Side" N = 0;
iter0_r(W2);									// Initial condictions for "Right Side" N = 0;

//  Matrix |Q' dS' r'|f_N^+|Q dS r|: N = 0;
//  Eigenvalues: |Q dS r|H_N|Q dS r|: N = 0;
//  EigenVectors: |Q dS r||Q dS p|: N = 0;
//  Dim(0); 

std::cout << "Done! Iteraction 0 complete!" << std::endl<< std::endl;

// 						ITERATIVE SOLUTION N > 0
for(int ni = 1; ni <= N_max; ni ++){   
	nq = 2*(ni+2)+1; 								// Total number of possibles q. 
	ns = ni + 2 +1; 								// Total number of possibles ds.

	dimen = new int*[nq]; 							// Matrix to save the basis in the Left sector.
	for(int i=0; i< nq; i++){
		dimen[i] = new int[ns];
		for(int j=0; j< ns; j++){
			dimen[i][j] = 0;
		}
	}

	iterN_L(ni, dimen_p, dimen);						// Solve the 'Left Side' hamiltonian for N=ni
	iterN_R(ni, dimen_p, dimen);						// Solve the 'Right Side' hamiltonian for N=ni

	nq = 2*(ni-1+2)+1; 							// Delete the dimen_p matrix. 
	for (int q = 0; q < nq; q++){
		delete[] dimen_p[q];
	}
	delete[] dimen_p;
	dimen_p = NULL;
	dimen_p = dimen;								// Left: New imput of number of basis
										 

	std::cout << "Done! Iteraction " << ni <<  " complete!" << std::endl<< std::endl<< std::endl;

	float time_ = difftime(time(NULL), beg_time);

	//std::cout << "The code takes " << time_  <<  " seconds util this point." << std::endl << std::endl;
}
//*/

time_t end_time = time(NULL);

float time_ = difftime(end_time, beg_time) /60;    					// Converte in seconds to minutes

std::cout << "The code is executed in " << time_  <<  " minutes." << std::endl;

return 0;
}


// Até a iteração N = 3, todas as autoenergias foram corrigidas. 
// 09/02/2022 até o N = 30, tempo de execução: 1.0 min. (Cada setor)
// 09/02/2022 até o N = 30, tempo de execução: 40 min (Ambos os setores + projeção).
// Projeções conferidas até N= 30. Mesmos valores de W1 e W2.
// Energias e projeções conferidas em um caso particular (anotado no caderno) até N = 2. 
// Energia de corte ultravioleta E_uv inserido no Bra "Left".


