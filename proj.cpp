#include "proj.h"
static double*** projection_;
static double*** projection_p_;

void projection_start(int N){
int nq = 2*(N+2) + 1;
int ns = 1*(N+2) + 1;

// Here i'm creating the pointer to save the sectors (Q,dS)
projection_ = new double**[nq];								
for(int q = 0; q < nq; q++){					
	projection_[q] = new double*[ns];
}						

}

// Here i'm allocating memory for the (Q,dS) sector. I need delete and re-alloc memory every iteration. 
void projection_alloc_memory(int q, int ds, long nk){

projection_[q][ds] = new double[nk];							

}

// Functions to write and read in the matrix. 

void projection_write(int q, int ds, long k, double value) {
	
projection_[q][ds][k] = value;

}

double projection_read(int q, int ds, long k) {

return projection_[q][ds][k];

}

// Delete all elements from the matrix.

void projection_delete(int N, int **dimen){

int nq = 2*(N+2) +1; 								
int ns = 1*(N+2) +1; 									
for(int q = 0; q < nq; q++){
	for(int ds=0; ds < ns; ds++){
		int dim  = dimen[q][ds];
		if (dim>0){
			delete[] projection_[q][ds];
		}
	}
	delete[] projection_[q];
}
delete[] projection_;
projection_ = NULL;
}

void save_projection(int N, int **dimen_){

int nq = 2*(N + 2) +1;						// Past nq
int ns = 1*(N + 2) +1; 						// Past ns

projection_p_ = new double**[nq];								
for(int q = 0; q < nq; q++){					
	projection_p_[q] = new double*[ns];
	for (int ds = 0; ds < ns; ds++){
		int dim = dimen_[q][ds];
		if (dim > 0){
			projection_p_[q][ds] = new double[(long) dim*dim];
			for(long k = 0; k < (long) dim*dim; k++){
				projection_p_[q][ds][k] = projection_read(q,ds,k);
			}
		}	
	}
}		

}


double save_projection_read(int q, int ds, long k){

return projection_p_[q][ds][k];
}

void save_projection_delete(int N, int **dimen_){

int nq = 2*(N +2) +1; 								
int ns = 1*(N +2) +1;
 									
for(int q = 0; q < nq; q++){
	for(int ds=0; ds < ns; ds++){
		int dim  = dimen_[q][ds];
		if (dim>0){
			delete[] projection_p_[q][ds];
		}
	}
	delete[] projection_p_[q];
}
delete[] projection_p_;
projection_p_ = NULL;
}



