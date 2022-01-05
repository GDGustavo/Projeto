#include "mel.h"
// This function will salve somes bigs matrixs for the NRG code. 
// 

static double ***mel_ne_; 								//|Q+1 dS+1 r'|f_(N)^+|Q S r|
static double ***mel_nw_;								//|Q+1 dS-1 r'|f_(N)^+|Q S r|
static double ***eigen_erg_;							//     |Q dS r'|H_{N}|Q s r|
static double ***eigen_vect_;							//       |Q dS r'|Q dS p'| 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 			       |Q+1 dS+1 r'|f_(N)^+|Q S r| and |Q+1 dS-1 r'|f_(N)^+|Q S r|
//This function starts the both mel_ne and mel_nw matrix.
void mel_start(int N){

int nq = 2*(N + 2) +1;
int ns = 1*(N + 2) +1; 

// Here i'm creating the pointer to save the sectors (Q,dS)
mel_ne_ = new double**[nq];							
mel_nw_ = new double**[nq];	
for(int q = 0; q < nq; q++){					
	mel_ne_[q] = new double*[ns];
	mel_nw_[q] = new double*[ns];
}										

}

// MEL_NE FUNCTIONS; 


// Here i'm allocating memory for the (Q,dS) sector. I need delete and re-alloc memory every iteration. 

void mel_ne_alloc_memory(int q, int ds, long nk){

mel_ne_[q][ds] = new double[nk];							

}

// Functions to write and read in the matrix. 

void mel_ne_write(int q, int ds, long k, double value) {
	
mel_ne_[q][ds][k] = value;

}

double mel_ne_read(int q, int ds, long k) {

return mel_ne_[q][ds][k];
 
}

// Delete all elements from the matrix.

void mel_ne_delete(int N, int **dimen){

int nq = 2*(N+2)+1; 								
int ns =  (N+2) +1; 									
for(int q = 0; q < nq-1; q++){
	for(int ds=0; ds < ns -1; ds++){
		int dim  = dimen[q][ds];
		int dim2 = dimen[q+1][ds+1]; 
		if ((dim>0)&&(dim2>0)){
			delete[] mel_ne_[q][ds];
		}
	}
	delete[] mel_ne_[q];
}
delete[] mel_ne_;
mel_ne_ = NULL;
}


// MEL_NW FUNCTIONS; 


// Here i'm allocating memory for the (Q,dS) sector. I need delete and re-alloc memory every iteration. 

void mel_nw_alloc_memory(int q, int ds, long nk){

mel_nw_[q][ds] = new double[nk];							

}


// Functions to write and read in the matrix. 

void mel_nw_write(int q, int ds, long k, double value) {
	
mel_nw_[q][ds][k] = value;

}

double mel_nw_read(int q, int ds, long k) {

return mel_nw_[q][ds][k];

}

// Delete all elements from the matrix.

void mel_nw_delete(int N, int **dimen){

int nq = 2*(N+2)+1; 								
int ns = 1*(N+2)+1; 									
for(int q = 0; q < nq-1; q++){
	for(int ds=1; ds < ns; ds++){
		int dim  = dimen[q][ds];
		int dim2 = dimen[q+1][ds-1]; 
		if ((dim>0)&&(dim2>0)){
			delete[] mel_nw_[q][ds];
		}
	}
	delete[] mel_nw_[q];
}
delete[] mel_nw_;
mel_nw_ = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 						EIGENVALUES AND EIGENVECTORS FUNCTIONS

void eigen_start(int N){

int nq = 2*(N + 2) +1;
int ns =  (N + 2) + 1; 

// Here i'm creating the pointer to save the sectors (Q,dS)
eigen_erg_  = new double**[nq];							
eigen_vect_ = new double**[nq];	
for(int q = 0; q < nq; q++){						
	eigen_erg_[q]  = new double*[ns];
	eigen_vect_[q] = new double*[ns];
}										
}


void eigen_delete(int N, int **dimen){

int nq = 2*(N+2)+1; 								
int ns = 1*(N+2)+1;								

for(int q = 0; q < nq; q++){
	for(int ds=0; ds < ns; ds++){
		int dim  = dimen[q][ds];
		if (dim>0){
			delete[] eigen_vect_[q][ds];
			delete[] eigen_erg_[q][ds];
		}
	}
	delete[] eigen_vect_[q];
	delete[] eigen_erg_[q];
}

delete[] eigen_vect_;
delete[] eigen_erg_;

eigen_vect_ = NULL;
eigen_erg_  = NULL;
}

void eigen_erg_alloc_memory(int q, int ds, long nk){

eigen_erg_[q][ds] = new double[nk];							

}

void eigen_vect_alloc_memory(int q, int ds, long nk){

eigen_vect_[q][ds] = new double[nk];						

}


// Functions to write and read in the matrix of the EIGENVALUES. 

void eigen_erg_write(int q, int ds, long k, double value) {
	
eigen_erg_[q][ds][k] = value;

}

double eigen_erg_read(int q, int ds, long k) {

return eigen_erg_[q][ds][k];

}


// Functions to write and read in the matrix of the EIGENVECTORS. 

void eigen_vect_write(int q, int ds, long k, double value) {
	
eigen_vect_[q][ds][k] = value;

}

double eigen_vect_read(int q, int ds, long k) {

return eigen_vect_[q][ds][k];
}
