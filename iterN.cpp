#include "iterN.h"

// For interactions more or equal than 1, we are supposedly able to find the bases and solve the hamiltonian, using the informations: 
// DIMEN_(N-1)[Q;dS] - How many states have the last iteraction for each (Q;dS) sector. 
// The matrix |Q' S' r'|f_(N-1)^+|Q S r| for the (N-1) iteraction. Since this matrix is importante to find the Hamiltonian H_N;
// The eigenvalues and eigenvectors for the (N-1) iteraction.

static int **NS_;
static int **NE_;
static int **NN_;
static int **NW_;

// Populating the primitive base; 
void Populate_N(int N, int **dimen_){
// DIMEN = DIM(N-1); 

int nq_past = 2*(N-1+2)+1;
int ns_past = 1*(N-1+2)+1;

// In the mathematical bases: 
// Operations Sulth |(Q-1); dS |_N <= 	|Q ; dS|_(N-1)
// Operations East  | Q ;(dS+1)|_N <= 	|Q ; dS|_(N-1)
// Operations North |(Q+1); dS |_N <= 	|Q ; dS|_(N-1)
// Operations West  | Q ;(dS-1)|_N <= 	|Q ; dS|_(N-1)

// But in the computaional bases: QcN = QN + [N+2] => (Q + [(N-1)+2])+1, then:
// Operations Sulth |Q  ; dS   |_N <= 	|Q ; dS|_(N-1)
// Operations East  |Q+1;(dS+1)|_N <= 	|Q ; dS|_(N-1)
// Operations North |Q+2; dS   |_N <= 	|Q ; dS|_(N-1)
// Operations West  |Q+1;(dS-1)|_N <= 	|Q ; dS|_(N-1)

for(int q = 0; q < nq_past; q++){
	for (int ds = 0; ds < ns_past; ds++){
		int dim = dimen_[q][ds];
		NS_[q][ds] += dim;
		NE_[q+1][ds+1] += dim;
		NN_[q+2][ds] += dim; 
		if(ds > 0){
			NW_[q+1][ds-1] += dim;
		}	
	}
}

// Done! Bases are already ready to use.
}

// Creating a new function to find the gender of  states {(q, ds, p)} in the primitive Bases;
// gen_ = 0; gen. South;
// gen_ = 1; gen. East;
// gen_ = 2; gen. North;
// gen_ = 3; gen. West;
// gen_ =-1; error;

int genreN(int q, int ds, int p){

int N1 = NS_[q][ds];     // NS

int gen_ = -1;

if (p <= N1){ gen_ = 0;}
else {
	int N2 = N1 + NE_[q][ds];// NS + NE
	if ((N1<p)&&(p<= N2)){
		gen_ =1;
	}
	else {
		int N3 = N2 + NN_[q][ds];// NS + NE + NN
		if ((N2<p)&&(p<= N3)){
			gen_ =2;
		}
		else {
			int N4 = N3 + NW_[q][ds];// NS + NE + NN + NW
			if ((N3<p)&&(p<= N4)){
				gen_ =3;
			}
		}
	}
}

return gen_;}

// To Find the charge, spin and r of the Father state with genre = {S,E,N,W}. 
/* 
int find_charge_spin_r(int q, int ds, int genre, int p, int w){

int result[3];
int r = -1;

result[0] = -1;  
result[1] = -1;
result[2] = -1;

if (genre == 0){ 						// Sulth
	result[0] = q;
	result[1] = ds;
	result[2] = p; 
}
else{	
	if(genre == 1){					// East
		result[0] = q - 1;
		result[1] = ds- 1;
		result[2] = p - NS_[q][ds];
	}
	else{
		if(genre == 2){				// North
			result[0] = q - 2;
			result[1] = ds;
			result[2] = p - NS_[q][ds] - NE_[q][ds];
		}
		else{
			if(genre == 3){			// west
				result[0] = q - 1;
				result[1] = ds +1;
				result[2] = p - NS_[q][ds] - NE_[q][ds] - NN_[q][ds];
			}
		}		
	}
}

if (w==0){r = result[0];}
if (w==1){r = result[1];}
if (w==2){r = result[2];}

return r;
// q'    = result[0]; To access the fathers charge. 
// ds'   = result[1]; To access the fathers spin.
// r'    = result[2]; To access the fathers inherent number. 
}  
//*/

int find_charge(int q, int genre){
int r = -1;
if (genre == 0){ 						// Sulth
	r = q;
}
else{	
	if(genre == 1){					// East
		r = q - 1;
	}
	else{
		if(genre == 2){				// North
			r = q - 2;
		}
		else{
			if(genre == 3){			// west
				r = q - 1;
			}
		}		
	}
}
return r; // q' = result[0]; To access the fathers charge.  
}  

int find_spin(int ds, int genre){

int r = -1;

if (genre == 0){ 						// Sulth
	r = ds;
}
else{	
	if(genre == 1){					// East
		r = ds- 1;
	}
	else{
		if(genre == 2){				// North
			r = ds;
		}
		else{
			if(genre == 3){			// west
				r = ds +1;
			}
		}		
	}
}

return r; // ds'   = result[1]; To access the fathers spin.
}  

int find_father(int q, int ds, int genre, int p){

int r = -1;

if (genre == 0){ 						// Sulth
	r = p; 
}
else{	
	if(genre == 1){					// East
		r = p - NS_[q][ds];
	}
	else{
		if(genre == 2){				// North
			r = p - NS_[q][ds] - NE_[q][ds];
		}
		else{
			if(genre == 3){			// west
				r = p - NS_[q][ds] - NE_[q][ds] - NN_[q][ds];
			}
		}		
	}
}
return r; // r'    = result[2]; To access the fathers inherent number. 
}  

// THE FUNCTION TO SOLVE THE ITERACTION N. 

void iterN_L(int N, int **dimen_p_, int **dimen_){

std::cout << "Diagonalization process for 'Left Side' N=" << N << " is starting now..."  << std::endl<< std::endl;


int nq = 2*(N+2)+1; // charge: {-(N+2),...,-1, 0, 1,..., (N+2)} => {0,1,..,2(N+2)}
int ns = 1*(N+2)+1;  // double spin: {0,1,..(N+2)} => {0,1,..,(N+2)}

// Creating the Basis and Making all matrix elements zero.
NS_ = new int*[nq];
NE_ = new int*[nq];
NN_ = new int*[nq];
NW_ = new int*[nq]; 
for(int i=0; i< nq; i++){
	NS_[i] = new int[ns];
	NE_[i] = new int[ns];
	NN_[i] = new int[ns];
	NW_[i] = new int[ns];
	for(int j=0; j< ns; j++){
			NS_[i][j] = 0;
			NE_[i][j] = 0;
			NN_[i][j] = 0;
			NW_[i][j] = 0;
	}
}

Populate_N(N, dimen_p_);
double lamb_ = lamb();								    	// Reading Lambda from parameters.txt.
double E_uv_ = (double) E_uv();
double aux = (1-pow(lamb_, (float) -N-2))/sqrt((1-pow(lamb_, (float) -2*(N-1)-1))*(1-pow(lamb_,(float) -2*(N-1)-3)));	
double D_N = (1-pow(lamb_, (float) -1))*(pow(lamb_, (float) -(N-1)/2))/log(lamb_);		//
//double t_N = aux*D_N;									// Calculating the coupling t_(N-1).
double t_N = aux; 									// Att, Scaled Hamiltonian.

std::cout << "t_{N-1}/D_N= " << '\t' << t_N <<std::endl;
std::cout << "D_N = "<< '\t' << '\t' << D_N <<std::endl;
std::cout << "t_{N-1}= " << '\t' <<  t_N*D_N <<std::endl;
std::cout << "E_uv_= " << '\t' <<  E_uv_ <<std::endl;


// H_N[p',p] = H_{N-1}[p',p] + t_{N-1}*M_N[p',p] + t_{N-1}*M_N[p,p']; 			Non scaled Hamiltonian
// H_N[p',p] = sqrt(\Lambda)*H_{N-1}[p',p] + t'_{N-1}*M_N[p',p] + t'_{N-1}*M_N[p,p']; 	Scaled Hamiltonian

// Creating here the 2D pointer to hamiltonian sector (q,ds). 
double ***HN_ = new double**[nq];
for(int q = 0; q < nq; q++){					
	HN_[q] = new double*[ns];
}

for(int q = 0; q < nq; q++){
	for (int ds = 0; ds < ns; ds++){
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		if (dim > 0){
			HN_[q][ds] = new double[(long) (dim+1)*dim/2];			// Allocating memory to sector (q,ds).
			for(int p2=0; p2 < dim; p2++){
				for(int p1=0; p1<= p2; p1++){
					int g1 = genreN(q,ds,p1+1);			// g1(p1)
					int g2 = genreN(q,ds,p2+1);			// g2(p2)
					int q1 = find_charge(q,g1);			// |q1' ds1' r1|
					int ds1= find_spin(ds,g1);				// 
					int r1 = find_father(q,ds,g1,p1+1);			//
					int q2 = find_charge(q,g2);			// |q2' ds2' r2|
					int ds2= find_spin(ds,g2);				// 
					int r2 = find_father(q,ds,g2,p2+1);			//
					int dim1 = dimen_p_[q1][ds1];				
					int dim2 = dimen_p_[q2][ds2];
					long k_mel_1 = (r2-1)*dim1+(r1-1);
					long k_mel_2 = (r1-1)*dim2+(r2-1); 
					long k = (p2*(p2+1)/2) + p1; 			// Memory address starting in k=0;
					HN_[q][ds][k] = 0; 
					if ((g1==g2)&&(q1==q2)&&(r1==r2)&&(ds1==ds2)){			// Diagonal terms.
					//HN_[q][ds][k]= eigen_erg_read(q1,ds1,(long) r1-1);
					HN_[q][ds][k]= sqrt(lamb_)*eigen_erg_read(q1,ds1,(long) r1-1);  	// Att, scaled H
					}
					if ((g2==0)&&(g1==1)){						// SE terms.
					HN_[q][ds][k] = t_N*mel_ne_read(q1,ds1,k_mel_1);
					}
					if ((g2==1)&&(g1==0)){						// ES terms. 
					HN_[q][ds][k] = t_N*mel_ne_read(q2,ds2,k_mel_2);
					}					
					if ((g2==1)&&(g1==2)){						// EN terms. 
					HN_[q][ds][k]=(sqrt(ds)/sqrt(ds+1))*t_N*mel_nw_read(q1,ds1,k_mel_1);
					}		
					if ((g2==2)&&(g1==1)){						// NE terms. 
					HN_[q][ds][k]=(sqrt(ds)/sqrt(ds+1))*t_N*mel_nw_read(q2,ds2,k_mel_2);
					}								
					if ((g2==0)&&(g1==3)){						// SW terms. 
					HN_[q][ds][k]= t_N*mel_nw_read(q1,ds1,k_mel_1);
					}	
					if ((g2==3)&&(g1==0)){						// WS terms. 
					HN_[q][ds][k]=t_N*mel_nw_read(q2,ds2,k_mel_2);
					}				
					if ((g2==3)&&(g1==2)){						// WN terms. 
					HN_[q][ds][k]=(-sqrt(ds+2)/sqrt(ds+1))*t_N*mel_ne_read(q1,ds1,k_mel_1);
					}	
					if ((g2==2)&&(g1==3)){						// NW terms. 
					HN_[q][ds][k]=(-sqrt(ds+2)/sqrt(ds+1))*t_N*mel_ne_read(q2,ds2,k_mel_2);
					}
				} // end for p1		
			} // end for p2			
		} // end if dim> 0	
	} // end for ds
} // end for q

// Here It's cleaning the memory, since we don't need keep the information anymore. 

eigen_delete(N-1,dimen_p_);  // Delete all elements saved in eigen matrix
mel_nw_delete(N-1,dimen_p_); // Delete all elements saved in mel_nw_
mel_ne_delete(N-1,dimen_p_); // Delete all elements saved in mel_ne_

eigen_start(N);

// Ok. The Hamiltonian is ready to diagonalization process. 
for (int q=0; q < nq; q++) {
	for (int ds=0; ds < ns; ds++) {
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		dimen_[q][ds] = dim;
		if(dim > 0) {
			std::cout <<'\t'<<'\t'<<'\t' << "["<< (q - N - 2) << ";" << ds;
			std::cout << "] Sector;"<<std::endl<<"dim ="<< dim << std::endl;
			double *eigen_values = new double[dim];                        	// Matrix to find the E-Energies
			double **eigen_vectors = new double*[dim];		               	// Matrix to find the E-Vectors
			long Nel_max = (1+dim)*dim/2;	                      		// Maximum number of elements of LTM
			double* Hamiltonian  = new double[Nel_max];                     	// H[q,ds] to use in diagonalization
			double* Hamiltonian_ = new double[Nel_max];                     	// H[q,ds] to use in diagonalization
			for (long k=0; k<Nel_max; k++) {					// Save the values to diagolize
				Hamiltonian[k] = HN_[q][ds][k];
				Hamiltonian_[k]= HN_[q][ds][k]; 
			}
			delete[] HN_[q][ds];						// Delete the HN in the sector (q,ds)
			for (int k=0; k<dim; k++) {					// Make all elements zero
				eigen_values[k] = 0;
				eigen_vectors[k] = new double[dim];
				for(int h=0; h< dim; h++){
					eigen_vectors[k][h] = 0;
				}
			}
			eigen_values[0] = E_uv_;							// Cut-off Energy
			int ret = givens(dim, 0 , Hamiltonian_ , eigen_values, eigen_vectors, 0);	// Find the Cut-off dim
			delete[] Hamiltonian_;

			for (int k=0; k<dim; k++) {						// Make all elements zero
				eigen_values[k] = 0;
				for(int h=0; h< dim; h++){
					eigen_vectors[k][h] = 0;
				}
			}

			givens(dim, abs(ret) , Hamiltonian , eigen_values, eigen_vectors, 1);		// Solve the H with cut-off
			delete[] Hamiltonian;
			Hamiltonian = NULL;

			int dim_c = abs(ret);							// Dim bellow  the cut-off
			dimen_[q][ds] = dim_c;							// Save the Dim[q.ds]

			eigen_erg_alloc_memory(q,ds,(long) dim_c);				// Alloc memory to save E-Energies
			eigen_vect_alloc_memory(q,ds,(long) dim_c*dim);			// Alloc memory to save E-Vectors

			std::cout << "Ultraviolet Cut-off Energy (non scaled): " << D_N*E_uv_ << std::endl;
			std::cout << "Number of Basis bellow the ultraviolet cut-off Energy: " << ret << std::endl;	
			std::cout << "Eigen values (Non scaled): " << std::endl;			
			for (long k=0; k<dim_c; k++) {
				eigen_erg_write(q,ds,k,eigen_values[k]);
				std::cout << D_N*eigen_values[k] << ";";
			}
			delete[] eigen_values;
			eigen_values = NULL;
			//std::cout << std::endl << "Eigen vectors matrix:" << std::endl;
			for (int i=0; i<dim_c; i++) {						// Line 0, 1... dim_cut_off
				for (int j=0; j<dim; j++) { 					// Collum 0, 1 ... dim
			//		std::cout << eigen_vectors[i][j] << ";";
					eigen_vect_write(q,ds,(long) i*dim+j,eigen_vectors[i][j]); 
				}
			//	std::cout << std::endl;
				delete[] eigen_vectors[i];
			}
			delete[] eigen_vectors;
			eigen_vectors= NULL;
			std::cout << std::endl<< std::endl;
		} //end if dim > 0
	} //end for ds
	delete[] HN_[q];
} // end for q
delete[] HN_;
HN_ = NULL;

std::cout << "Diagonalization process already finished for 'Left Side' N = "<< N<<"!" << std::endl<< std::endl<< std::endl<< std::endl;


mel_start(N);


for(int q=0;q<(nq-1);q++){
	//|Q+1 dS+1 r2|(f_N^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS+1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS+1 p2|(f_N^+)|Q dS p1|
	for(int ds=0;ds<(ns-1);ds++){
		int dim = dimen_[q][ds];
		int dim2 = dimen_[q+1][ds+1]; 
		if (dim*dim2>0){
			int N11 = NS_[q][ds];		// Sup limit for g1 = 0
			int N12 = N11 + NE_[q][ds]; 	// Sup limit for g1 = 1
			int N13 = N12 + NN_[q][ds];	// Sup Limit for g1 = 2
			int N14 = N13 + NW_[q][ds];	// Sup Limit for g1 = 3 (Dim without the dim Cut-off) 
			int N21 = NS_[q+1][ds+1];		// Sup Limit for g2 = 0
			int N22 = N21 + NE_[q+1][ds+1];	//
			int N23 = N22 + NN_[q+1][ds+1];	//
			int N24 = N23 + NW_[q+1][ds+1];	// Sup Limit for g2 = 3 (Dim without the dim Cut-off) 
			mel_ne_alloc_memory(q,ds,(long) dim*dim2);						// 
			for(long k=0;k<dim*dim2;k++){
				double sum = 0;
				for(int p1=0; p1<N11; p1++){ 				// g1(p1) = 0;
					for(int p2=N21; p2<N22; p2++){			// g2(p2) = 1;
						//int l1 = find_father(q,ds,0,p1+1);	// Father 1? 
						//int l2 = find_father(q+1,ds+1,1,p2+1);	// Father 2?
						int l1 = p1;
						int l2 = p2 - N21;
						if (l2==l1){				// delta(l1,l2)
							int r2 = k/dim;      		// line
							int r1 = k - r2*dim; 		// Collum
							double aux2=eigen_vect_read(q+1,ds+1,(long)r2*N24 +p2);
							double aux1=eigen_vect_read(q,ds,(long)r1*N14 +p1);
							sum = sum + aux2*aux1;
						}
					}  						
				}
				for(int p1=N13; p1< N14; p1++){ 				// g1(p1) = 3;
					for(int p2=N22; p2<N23; p2++){			// g2(p2) = 2;
						//int l1 = find_father(q,ds,3,p1+1);	// Father 1? 
						//int l2 = find_father(q+1,ds+1,2,p2+1);	// Father 2?
						int l1 = p1 - N13;
						int l2 = p2 - N22;
						if (l2==l1){				// delta(l1,l2)
							int r2 = k/dim;      		// line
							int r1 = k - r2*dim; 		// Collum
							double aux2=eigen_vect_read(q+1,ds+1,(long)r2*N24 +p2);
							double aux1=eigen_vect_read(q,ds,(long)r1*N14 +p1);
							sum = sum + sqrt(ds+1)/sqrt(ds+2)*aux2*aux1;
						}
					}  					
				}
				mel_ne_write(q,ds,k,sum);
			} //end for k
		} // end if dim*dim2
	} // end for ds
	//|Q+1 dS-1 r2|(f_N^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS-1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS-1 p2|(f_N^+)|Q dS p1|
	for(int ds=1;ds<ns;ds++){
		int dim = dimen_[q][ds];
		int dim2 = dimen_[q+1][ds-1];
		if (dim*dim2>0){
			int N11 = NS_[q][ds];		// Sup limit for g1 = 0
			int N12 = N11 + NE_[q][ds]; 	// Sup limit for g1 = 1
			int N21 = NS_[q+1][ds-1];		// Sup Limit for g2 = 0
			int N13 = N12 + NN_[q][ds];	// Sup Limit for g1 = 2
			int N14 = N13 + NW_[q][ds];	// Sup Limit for g1 = 3 (Dim without the dim Cut-off)
			int N22 = N21 + NE_[q+1][ds-1];	//
			int N23 = N22 + NN_[q+1][ds-1];	//
			int N24 = N23 + NW_[q+1][ds-1];	// Sup Limit for g2 = 3 (Without the dim Cut-off) 
			mel_nw_alloc_memory(q,ds,(long) dim*dim2);						// 
			for(long k=0;k<dim*dim2;k++){
				double sum = 0;
				for(int p1=0; p1<N11; p1++){ 				// g1(p1) = 0;
					for(int p2=N23; p2< N24; p2++){			// g2(p2) = 3;
						int l1 = p1;
						int l2 = p2 - N23;
						if (l2==l1){				// delta(l1,l2)
							int r2 = k/dim;      		// line
							int r1 = k - r2*dim; 		// Collum
							double aux2=eigen_vect_read(q+1,ds-1,(long)r2*N24 +p2);
							double aux1=eigen_vect_read(q,ds,(long)r1*N14 +p1);
							sum = sum + aux2*aux1;
						}
					}  						
				}
				for(int p1=N11; p1< N12; p1++){ 				// g1(p1) = 1;
					for(int p2=N22; p2<N23; p2++){			// g2(p2) = 2;
						int l1 = p1 - N11;
						int l2 = p2 - N22;
						if (l2==l1){				// delta(l1,l2)
							int r2 = k/dim;      		// line
							int r1 = k - r2*dim; 		// Collum
							double aux2=eigen_vect_read(q+1,ds-1,(long)r2*N24 +p2);
							double aux1=eigen_vect_read(q,ds,(long)r1*N14 +p1);
							sum = sum -sqrt(ds+1)/sqrt(ds)*aux2*aux1;
						}
					}  						
				}
				mel_nw_write(q,ds,k,sum);
			} // end for k 
		} // end if dim*dim2
	} //end for ds

} //end for q

for (int q = 0; q < nq; q++){
	delete[] NS_[q];
	delete[] NE_[q];
	delete[] NN_[q];
	delete[] NW_[q];
}
delete[] NS_;
delete[] NE_;
delete[] NN_;
delete[] NW_;
NS_ = NULL;
NE_ = NULL;
NN_ = NULL;
NW_ = NULL;

}
