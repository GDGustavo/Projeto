#include "iterN_r.h"

static int **NS_;
static int **NE_;
static int **NN_;
static int **NW_;

// Populating the primitive base; 
void Populate_N_R(int N, int **dimen_){

int nq_past = 2*(N-1+2)+1;
int ns_past = 1*(N-1+2)+1;

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

//Creating a new function to find the gender of  states {(q, ds, p)} in the primitive Bases;
// gen_ = 0; gen. South;
// gen_ = 1; gen. East;
// gen_ = 2; gen. North;
// gen_ = 3; gen. West;
// gen_ =-1; error;

int genreN_R(int q, int ds, int p){

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


int find_charge_R(int q, int genre){
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

int find_spin_R(int ds, int genre){

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

int find_father_R(int q, int ds, int genre, int p){

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

void iterN_R(int N, int **dimen_p_, int **dimen_out_){

std::cout << "Diagonalization process for 'Right Side' N= "<< N << " is starting now..."  << std::endl<< std::endl;


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

Populate_N_R(N, dimen_p_);
double lamb_ = lamb();								    // It's reading Lambda from parameters.txt.
double aux1 = (1-pow(lamb_,-N-2))/sqrt((1-pow(lamb_,-2*(N-1)-1))*(1-pow(lamb_,-2*(N-1)-3)));
double t_N = aux1*(1-pow(lamb_,-1))*(pow(lamb_,-(N-1)/2))/log(lamb_); 			    // It's calculating the coupling t_(N-1).
std::cout << "t_{N-1}= "<< t_N <<std::endl;

// H_N[p',p] = H_(N-1)[p',p] + t_(N-1)*M_N[p',p] + t_(N-1)*M_N[p,p'];

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
					int g1 = genreN_R(q,ds,p1+1);				// g1(p1)
					int g2 = genreN_R(q,ds,p2+1);				// g2(p2)
					int q1 = find_charge_R(q,g1);			// |q1' ds1' r1|
					int ds1= find_spin_R(ds,g1);			// 
					int r1 = find_father_R(q,ds,g1,p1+1);		//
					int q2 = find_charge_R(q,g2);			// |q2' ds2' r2|
					int ds2= find_spin_R(ds,g2);			// 
					int r2 = find_father_R(q,ds,g2,p2+1);		//
					int dim1 = dimen_p_[q1][ds1];				
					int dim2 = dimen_p_[q2][ds2];
					long k_mel_1 = (r2-1)*dim1+(r1-1);
					long k_mel_2 = (r1-1)*dim2+(r2-1); 
					long k = (p2*(p2+1)/2) + p1; 			// Memory address starting in k=0;
					HN_[q][ds][k] = 0; 
					if ((g1==g2)&&(q1==q2)&&(r1==r2)&&(ds1==ds2)){			// Diagonal terms.
					HN_[q][ds][k]= eigen2_erg_read(q1,ds1,(long) r1-1); 
					}
					if ((g2==0)&&(g1==1)){						// SE terms.
					HN_[q][ds][k] = t_N*mel2_ne_read(q1,ds1,k_mel_1);
					}
					if ((g2==1)&&(g1==0)){						// ES terms. 
					HN_[q][ds][k] = t_N*mel2_ne_read(q2,ds2,k_mel_2);
					}					
					if ((g2==1)&&(g1==2)){						// EN terms. 
					HN_[q][ds][k]=(sqrt(ds)/sqrt(ds+1))*t_N*mel2_nw_read(q1,ds1,k_mel_1);
					}		
					if ((g2==2)&&(g1==1)){						// NE terms. 
					HN_[q][ds][k]=(sqrt(ds)/sqrt(ds+1))*t_N*mel2_nw_read(q2,ds2,k_mel_2);
					}								
					if ((g2==0)&&(g1==3)){						// SW terms. 
					HN_[q][ds][k]= t_N*mel2_nw_read(q1,ds1,k_mel_1);
					}	
					if ((g2==3)&&(g1==0)){						// WS terms. 
					HN_[q][ds][k]=t_N*mel2_nw_read(q2,ds2,k_mel_2);
					}				
					if ((g2==3)&&(g1==2)){						// WN terms. 
					HN_[q][ds][k]=(-sqrt(ds+2)/sqrt(ds+1))*t_N*mel2_ne_read(q1,ds1,k_mel_1);
					}	
					if ((g2==2)&&(g1==3)){						// NW terms. 
					HN_[q][ds][k]=(-sqrt(ds+2)/sqrt(ds+1))*t_N*mel2_ne_read(q2,ds2,k_mel_2);
					}
				} // end for p1		
			} // end for p2			
		} // end if dim> 0	
	} // end for ds
} // end for q

// Here It's cleaning the memory, since we don't need keep the information anymore. 

eigen2_delete(N-1,dimen_p_);  // Delete all elements saved in eigen matrix
mel2_nw_delete(N-1,dimen_p_); // Delete all elements saved in mel_nw_
mel2_ne_delete(N-1,dimen_p_); // Delete all elements saved in mel_ne_

eigen2_start(N);

// Ok. The Hamiltonian is ready to diagonalization process. 
for (int q=0; q < nq; q++) {
	for (int ds=0; ds < ns; ds++) {
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		if(dim > 0) {
			std::cout <<'\t'<<'\t'<<'\t' << "["<< (q - N - 2) << ";" << ds;
			std::cout << "] Sector;"<<std::endl<<"dim ="<< dim << std::endl;
			eigen2_erg_alloc_memory(q,ds,(long) dim);
			eigen2_vect_alloc_memory(q,ds,(long) dim*dim);
			double *eigen_values = new double[dim];                        	// Changing the size
			double **eigen_vectors = new double*[dim];		               	// Changing the size
			long Nel_max = (1+dim)*dim/2;	                      		// Maximum number of elements in LTM
			double* Hamiltonian = new double[Nel_max];                     	// Change the size
			for (long k=0; k<Nel_max; k++) {					// Saving the values to diagolize
				Hamiltonian[k]= HN_[q][ds][k];
			}
			delete[] HN_[q][ds];						// Deteling the HN in the sector (q,ds)
			for (int k=0; k<dim; k++) {					// Making all elements zero
				eigen_values[k] = 0;
				eigen_vectors[k] = new double[dim];
				for(int h=0; h< dim; h++){
					eigen_vectors[k][h] = 0;
				}
			}
			givens(dim, dim , Hamiltonian , eigen_values, eigen_vectors, 0);	// Solving the H
			delete[] Hamiltonian;
			Hamiltonian = NULL;
			std::cout << "Eigen values: " << std::endl;
			for (long k=0; k<dim; k++) {
				std::cout << eigen_values[k] << ";";
				eigen2_erg_write(q,ds,k,eigen_values[k]); 
			}
			delete[] eigen_values;
			eigen_values = NULL;
			//std::cout << std::endl; << "Eigen vectors matrix:" << std::endl;
			for (int i=0; i<dim; i++) {
				for (int j=0; j<dim; j++) { 
			//		std::cout << eigen_vectors[i][j] <<"  " <<'\t' ;
					eigen2_vect_write(q,ds,(long) i*dim+j,eigen_vectors[i][j]); 
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

std::cout << "Diagonalization process already finished for 'Right Side' N = "<< N<<"!" << std::endl<< std::endl<< std::endl<< std::endl;


mel2_start(N);


for(int q=0;q<(nq-1);q++){
	//|Q+1 dS+1 r2|(f0^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS+1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS+1 p2|(f0^+)|Q dS p1|
	for(int ds=0;ds<(ns-1);ds++){
		int dim = dimen_out_[q][ds];
		int dim2 = dimen_out_[q+1][ds+1]; 
		if (dim*dim2>0){
			mel2_ne_alloc_memory(q,ds,(long) dim*dim2);						// 
			for(long k=0;k<dim*dim2;k++){
				double sum = 0;
					for(int p1=0; p1<dim; p1++){
						for(int p2=0; p2<dim2; p2++){
							int g1 = genreN_R(q,ds,p1+1);	 	// What is the genre p1? 
							int g2 = genreN_R(q+1,ds+1,p2+1);		// What is the genre p2?
							if ((g1==0)&&(g2==1)){
								int l1 = find_father_R(q,ds,g1,p1+1);		//Father 1?
								int l2 = find_father_R(q+1,ds+1,g2,p2+1);		//Father 2?
								if (l2==l1){				// delta(l1,l2)
									int r2 = k/dim;      		// line
									int r1 = k - r2*dim; 		// Collum
									double aux2=eigen2_vect_read(q+1,ds+1,(long)r2*dim2 +p2);
									double aux1=eigen2_vect_read(q,ds,(long)r1*dim +p1);
									sum = sum + aux2*aux1;
								}
							} 
							else{
								if ((g1==3)&&(g2==2)){
									int l1 = find_father_R(q,ds,g1,p1+1);	//Father 1? 
									int l2 = find_father_R(q+1,ds+1,g2,p2+1);	//Father 2?
									if (l2==l1){			//delta(l1,l2)
									int r2 = k/dim;      			// line
									int r1 = k - r2*dim; 			// Collum
									double aux2=eigen2_vect_read(q+1,ds+1,(long)r2*dim2 +p2);
									double aux1=eigen2_vect_read(q,ds,(long)r1*dim +p1);
									sum = sum + sqrt(ds+1)/sqrt(ds+2)*aux2*aux1;
									}
								
								}
							}
						}
					}
				mel2_ne_write(q,ds,k,sum);
			}// end for k  
		}// end if dim*dim2
	}// end for ds
	//|Q+1 dS-1 r2|(f0^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS-1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS-1 p2|(f0^+)|Q dS p1|
	for(int ds=1;ds<ns;ds++){
		int dim = dimen_out_[q][ds];
		int dim2 = dimen_out_[q+1][ds-1];
		if (dim*dim2>0){
			mel2_nw_alloc_memory(q,ds,(long) dim*dim2);						// 
			for(long k=0;k<dim*dim2;k++){
				double sum = 0; 
					for(int p1=0; p1<dim; p1++){
						for(int p2=0; p2<dim2; p2++){
							int g1 = genreN_R(q,ds,p1+1); 		// What is the genre p1? 
							int g2 = genreN_R(q+1,ds-1,p2+1);		// What is the genre p2?
							if ((g1==0)&&(g2==3)){
								int l1 = find_father_R(q,ds,g1,p1+1);		// Father 1?
								int l2 = find_father_R(q+1,ds-1,g2,p2+1);		// Father 2?
								if (l2==l1){
									int r2 = k/dim;      // line
									int r1 = k - r2*dim; // Collum
									double aux2=eigen2_vect_read(q+1,ds-1,(long)r2*dim2 +p2);
									double aux1=eigen2_vect_read(q,ds,(long) r1*dim +p1);
									sum = sum + aux2*aux1;
								}	
							}
							else { 
								if ((g1==1)&&(g2==2)){
									int l1 = find_father_R(q,ds,g1,p1+1);	// Father 1?
									int l2 = find_father_R(q+1,ds-1,g2,p2+1);	// Father 2?
									if (l2==l1){
									int r2 = k/dim;      // line
									int r1 = k - r2*dim; // Collum
									double aux2=eigen2_vect_read(q+1,ds-1,(long)r2*dim2 +p2);
									double aux1=eigen2_vect_read(q,ds,(long) r1*dim +p1);
									sum = sum -sqrt(ds+1)/sqrt(ds)*aux2*aux1;
									}	
								}
							}
						}
					}
				mel2_nw_write(q,ds,k,sum);
			}// end for k  
		}// end if dim*dim2 
 	}// end for ds
}// end for q


// Projetcions Beteween the "Right" and "Left" sectors. 

save_projection(N-1, dimen_p_);         // Save the past projection matrix 
projection_delete(N-1, dimen_p_);       // Delete the projection matrix from the iteraction N-1
projection_start(N);                  // Starting the adress in the sector (q, ds) to save the matrix projection

//std::cout <<"Bases Projection for 'Right Side' N = " << N <<";" << std::endl << std::endl;

for(int q=0; q<nq; q++){
	for(int ds = 0; ds < ns; ds++){
		int dim = dimen_out_[q][ds];
		if(dim> 0){
			projection_alloc_memory(q,ds, (long) dim*dim);
			for(long k = 0; k < dim*dim; k++){
				double sum = 0;
				for(int p_L = 0; p_L < dim; p_L ++){
					for(int p_R = 0; p_R < dim; p_R ++){
						int g_L = genreN_R(q, ds, p_L+1);
						int g_R = genreN_R(q, ds, p_R+1);
						if (g_L==g_R){
							int r_L = k/dim;					//Line
							int r_R = k - r_L*dim;				//Colum
							int q_L = find_charge_R(q,g_L); 
							int ds_L= find_spin_R(ds,g_L);
							int dm_L= dimen_p_[q_L][ds_L];
							int l_L = find_father_R(q,ds,g_L,p_L+1);
							int l_R = find_father_R(q,ds,g_R,p_R+1); 
							double aux_1 = eigen_vect_read(q,ds, (long) r_L*dim +p_L);
							double aux_2 = eigen2_vect_read(q,ds,(long) r_R*dim + p_R);
							double aux_3 = save_projection_read(q_L,ds_L,(long)(l_L-1)*dm_L +l_R-1); 
							sum = sum + aux_1*aux_2*aux_3;
						}
					}
				}
				if(abs(sum) < 0.0000000001){
					sum = 0;
				}
				projection_write(q,ds,k,sum);
				//std::cout <<"Proj["<<q-N-2<<";"<< ds<<"]("<<r_L + 1 << ";" <<r_R + 1 << ") = " <<'\t'; 
				//std::cout << projection_read(q,ds,k) << std::endl;
			} // end for k
		} // end if dim > 0
	} // end for ds
} // end for q

save_projection_delete(N-1, dimen_p_);  // Delete the saved matrix

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