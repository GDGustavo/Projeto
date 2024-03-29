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

void iterN_R(int N, int **dimen_p_, int **dimen_){

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
double lamb_ = lamb();								    	// Reading Lambda from parameters.txt.
double aux1 = (1-pow(lamb_,(float) -N-2))/sqrt((1-pow(lamb_, (float) -2*(N-1)-1))*(1-pow(lamb_, (float) -2*(N-1)-3)));
double D_N  = (1-pow(lamb_, (float) -1))*(pow(lamb_, (float) -(N-1)/2))/log(lamb_); 
double t_N  = aux1;								    	// Calculating the coupling t_(N-1).
double E_f  = 1e20;									// Fundamental Energy

// H_N[p',p] = H_(N-1)[p',p] + t_(N-1)*M_N[p',p] + t_(N-1)*M_N[p,p'];

// Creating here the 2D pointer to hamiltonian sector (q,ds). 
double ***HN_ = new double**[nq];
for(int q = 0; q < nq; q++){					
	HN_[q] = new double*[ns];
}

#pragma omp parallel for
for(int q = 0; q < nq; q++){
	for (int ds = 0; ds < ns; ds++){
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		if (dim > 0){
			HN_[q][ds] = new double[(long) (dim+1)*dim/2];			// Allocating memory to sector (q,ds).
			for(int p2=0; p2 < dim; p2++){
				for(int p1=0; p1<= p2; p1++){
					int g1 = genreN_R(q,ds,p1+1);			// g1(p1)
					int g2 = genreN_R(q,ds,p2+1);			// g2(p2)
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
					HN_[q][ds][k]= sqrt(lamb_)*eigen2_erg_read(q1,ds1,(long) r1-1); 
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
			double *eigen_values = new double[dim];                        	//
			double **eigen_vectors = new double*[dim];		               	//
			for (int k=0; k<dim; k++) {					//
				eigen_vectors[k] = new double[dim];
			}
			long Nel_max = (1+dim)*dim/2;	                      		// Maximum number of elements in LTM
			double* Hamiltonian = new double[Nel_max];                     	// Change the size
			for (long k=0; k<Nel_max; k++) {					// Saving the values to diagolize
				Hamiltonian[k]= HN_[q][ds][k];
			}
			delete[] HN_[q][ds];						// Deteling the HN in the sector (q,ds)
			
			int ret = dimen_[q][ds];						// Number of states bellow cut-off.
			givens(dim, ret , Hamiltonian , eigen_values, eigen_vectors, 0);	// Solving the H
			delete[] Hamiltonian;
			Hamiltonian = NULL;

			eigen2_erg_alloc_memory(q,ds,(long) ret);
			eigen2_vect_alloc_memory(q,ds,(long) ret*dim);

			for (long k=0; k< ret ; k++) {
				eigen2_erg_write(q,ds,k,eigen_values[k]);
				if (E_f > eigen_values[k]){
					E_f = eigen_values[k];
				}
			}
			delete[] eigen_values;
			eigen_values = NULL;
			for (int i=0; i< dim; i++) {
				int signal = 1;
				if (eigen_vectors[i][0] < 0){
					signal = -1;
				}
				else{
					if((eigen_vectors[i][0] == 0)&&(eigen_vectors[i][1]<0)){
						signal = -1;						
					}
				}
				for (int j=0; j<dim; j++) {
					if(i < ret){  
						eigen2_vect_write(q,ds,(long) i*dim+j,(double) signal*eigen_vectors[i][j]); 
					}				
				}
				delete[] eigen_vectors[i];
			}
			delete[] eigen_vectors;
			eigen_vectors= NULL;
		} //end if dim > 0
	} //end for ds
	delete[] HN_[q];
} // end for q
delete[] HN_;
HN_ = NULL;

// Subtracting the Fundamental energy from energy spectre.
std::cout << "Fundamental Energy (Non Escaled) for N = "<< N <<":" << '\t' << -D_N*(eigen2_erg_read(0,0,0)-E_f) << std::endl<< std::endl;
for (int q=0; q < nq; q++) {
	for (int ds=0; ds < ns; ds++) {
		int dim = dimen_[q][ds];
		if(dim > 0) {
			for (long k=0; k<dim; k++) {
				double Energy = eigen2_erg_read(q,ds,k) - E_f;
				eigen2_erg_write(q,ds,k,Energy);	        	
			}

		}// end if dim>0
	} //end for ds
}//end for q

// Printing the eigen energies and eigen vectors
for (int q=0; q < nq; q++) {
	for (int ds=0; ds < ns; ds++) {
		int dim = dimen_[q][ds];
		if((dim > 0)&&(abs(q-N-2) + ds <= 1)&&(ds==0)&&(q-N-2>=0)) {
			std::cout << "["<< (q - N - 2) << ";" << ds;
			std::cout << "] Sector"<< '\t' <<"dim ="<< dim << std::endl;
			std::cout << "Eigen values (Escaled): ";
	
			for (long k=0; k<dim; k++) {
				double Energy = eigen2_erg_read(q,ds,k);
				std::cout << "E[" << k+1 << "]=";        	
				std::cout << Energy << ";";         
			}
			if (N == fN_max()){
				std::ofstream file2;
				file2.open("Energy2.txt");
				//Saving the Energies 
				for (long k=0; k<dim; k++) {
					double Energy = eigen2_erg_read(q,ds,k);        
					file2 << std::setprecision(12) << Energy;
					file2 << std::endl;
				}
  				//Closing the file
  				file2.close();
			}

				/*  Printing The Vectors.
				std::cout << std::endl << "Eigen vectors matrix:" << std::endl;
				for (int i=0; i<dim; i++) {
					std::cout << '\t';
					for (int j=0; j<dim_t; j++) {
						long k =i*dim +j;
						double vector = eigen2_vect_read(q,ds,k);
						std::cout << vector << "; ";
					}
					std::cout << std::endl;
				}
				// */
				std::cout << std::endl<< std::endl;
		}// end if dim>0
	} //end for ds
}//end for q

if(N == fN_max()){

eigen2_delete(N,dimen_);  // Delete all elements saved in eigen matrix

}

if(N < fN_max()){

mel2_start(N);

for(int q=0;q<(nq-1);q++){
	//|Q+1 dS+1 r2|(f_N^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS+1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS+1 p2|(f_N^+)|Q dS p1|
	for(int ds=0;ds<(ns-1);ds++){
		int dim = dimen_[q][ds];
		int dim2 = dimen_[q+1][ds+1]; 
		if (dim*dim2>0){
			int N11 = NS_[q][ds];		// Sup limit for g1 = 0
			int N12 = N11 + NE_[q][ds]; 	// Sup limit for g1 = 1
			int N13 = N12 + NN_[q][ds];	// Sup Limit for g1 = 2
			int N14 = N13 + NW_[q][ds];	// Sup Limit for g1 = 3
			int N21 = NS_[q+1][ds+1];		// Sup Limit for g2 = 0
			int N22 = N21 + NE_[q+1][ds+1];	//
			int N23 = N22 + NN_[q+1][ds+1];	//
			int N24 = N23 + NW_[q+1][ds+1];	//
			mel2_ne_alloc_memory(q,ds,(long) dim*dim2);					//
			#pragma omp parallel for
			for(long k=0;k<dim*dim2;k++){
				double sum = 0;
				int r2 = k/dim;      						// line
				int r1 = k - r2*dim; 						// Collum
				for(int p1=0; p1<N11; p1++){ 			// g1(p1) = 0;	// g2(p2) = 1;
					int p2 = p1 + N21;
					if (p2<N22){						// delta(l1,l2)
						double aux2=eigen2_vect_read(q+1,ds+1,(long)r2*N24 +p2);
						double aux1=eigen2_vect_read(q,ds,(long)r1*N14 +p1);
						sum += aux2*aux1;
					}						
				}
				for(int p1=N13; p1< N14; p1++){ 			// g1(p1) = 3;	// g2(p2) = 2;	
					int p2 = (p1 - N13) + N22;					// delta(l1,l2)
					if (p2<N23){				
						double aux2=eigen2_vect_read(q+1,ds+1,(long)r2*N24 +p2);
						double aux1=eigen2_vect_read(q,ds,(long)r1*N14 +p1);
						sum += sqrt(ds+1)/sqrt(ds+2)*aux2*aux1;
					}					
				}
				mel2_ne_write(q,ds,k,sum);
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
			int N13 = N12 + NN_[q][ds];	// Sup Limit for g1 = 2
			int N14 = N13 + NW_[q][ds];	// Sup Limit for g1 = 3
			int N21 = NS_[q+1][ds-1];		// Sup Limit for g2 = 0
			int N22 = N21 + NE_[q+1][ds-1];	//
			int N23 = N22 + NN_[q+1][ds-1];	//
			int N24 = N23 + NW_[q+1][ds-1];	//
			mel2_nw_alloc_memory(q,ds,(long) dim*dim2);					// 
			#pragma omp parallel for
			for(long k=0;k<dim*dim2;k++){
				double sum = 0;
				int r2 = k/dim;      						// line
				int r1 = k - r2*dim; 						// Collum
				for(int p1=0; p1<N11; p1++){ 			// g1(p1) = 0;	// g2(p2) = 3;
					int p2 = p1 + N23; 					// delta(l1,l2)	
					if (p2<N24){				
						double aux2=eigen2_vect_read(q+1,ds-1,(long)r2*N24 +p2);
						double aux1=eigen2_vect_read(q,ds,(long)r1*N14 +p1);
						sum += aux2*aux1;
					}					
				}
				for(int p1=N11; p1< N12; p1++){ 			// g1(p1) = 1;	// g2(p2) = 2;
					int p2 = (p1 - N11) + N22;					// delta(l1,l2)
					if (p2<N23){				
						double aux2=eigen2_vect_read(q+1,ds-1,(long)r2*N24 +p2);
						double aux1=eigen2_vect_read(q,ds,(long)r1*N14 +p1);
						sum -= sqrt(ds+1)/sqrt(ds)*aux2*aux1;
					}  						
				}
				mel2_nw_write(q,ds,k,sum);
			} // end for k 
		} // end if dim*dim2
	} //end for ds
} //end for q

} // end if N < N_max

// The importants sectors for the projections are such that the values are between : 
int q_0 = 1;
if( fN_max()%2 == 0){
q_0 = 0;
}
int ds_0 = 0;

int soma_max = q_0 + ds_0 + (fN_max() - N); 
int soma_min = q_0 + ds_0 - (fN_max() - N);
int q_max = q_0 + (fN_max() - N); 
int q_min = q_0 - (fN_max() - N); 
int ds_max = ds_0 + (fN_max() - N);
int ds_min = ds_0 - (fN_max() - N);

// Projetcions Beteween the "Right" and "Left" sectors.


// Save the eigenvector in matrix here to acess more fast
double*** vect_L_ = new double**[nq];								
double*** vect_R_ = new double**[nq];
for(int q = 0; q < nq; q++){					
	vect_L_[q] = new double*[ns];
	vect_R_[q] = new double*[ns];
	for (int ds = 0; ds < ns; ds++){
		int dim = dimen_[q][ds];
		int dim2= NS_[q][ds]+NE_[q][ds]+NN_[q][ds]+NW_[q][ds];		
		if (dim > 0){
			vect_L_[q][ds] = new double[(long) dim*dim2];
			vect_R_[q][ds] = new double[(long) dim*dim2];
			for(long k = 0; k < (long) dim*dim2; k++){
				vect_L_[q][ds][k] = eigen_vect_read(q,ds,k);
				vect_R_[q][ds][k] = eigen2_vect_read(q,ds,k);
			}
		}	
	}
}		
eigen_vect_delete(N,dimen_); // Delete all elements saved in eigen matrix
eigen2_vect_delete(N,dimen_); // Delete all elements saved in eigen matrix

 // Save the past projection matrix 
double*** projection_past_ = new double**[nq-2];								
for(int q = 0; q < nq-2; q++){					
	projection_past_[q] = new double*[ns-1];
	for (int ds = 0; ds < ns-1; ds++){
		int dim = dimen_p_[q][ds];		
		if (dim > 0){
			projection_past_[q][ds] = new double[(long) dim*dim];
			//#pragma omp parallel for
			for(long k = 0; k < (long) dim*dim; k++){
				projection_past_[q][ds][k] = projection_read(q,ds,k);
			}
		}	
	}
}		
projection_delete(N-1, dimen_p_);       // Delete the projection matrix from the iteraction N-1
projection_start(N);                    // Starting the adress in the sector (q, ds) to save the matrix projection

for(int q=0; q<nq; q++){
	for(int ds = 0; ds < ns; ds++){
		int dim = dimen_[q][ds];
		if((N==fN_max())&&(((q-N-2)!=(q_0))||(ds != ds_0))){// we could only calculate in the ground state charge-spim 
			dim = 0;
		}
		int selec = 1;
		int charge = q - N - 2; 
		int soma = charge + ds; 
		if((soma <= soma_max)&&(soma >= soma_min)){
			if((charge <= q_max)&&(charge >= q_min)){
				if((ds <= ds_max)&&(ds >= ds_min)){
					selec = 0;
				}
			}
		}
		if(dim> 0){
			int N1 = NS_[q][ds];							// Sup limit for g1 = 0
			int N2 = N1 + NE_[q][ds]; 							// Sup limit for g1 = 1
			int N3 = N2 + NN_[q][ds];							// Sup Limit for g1 = 2
			int N4 = N3 + NW_[q][ds];							// Sup Limit for g1 = 3
			if (selec == 0){
				std::cout << "Starting Proj. ["<< (q - N - 2) << ";" << ds;
				std::cout << "] Sector, " <<"dim =="<< dim << std::endl;
				projection_alloc_memory(q,ds, (long) dim*dim);
				#pragma omp parallel for
				for(long k = 0; k < dim*dim; k++){
					double sum = 0;
					int r_L = k/dim;							// Line
					int r_R = k - r_L*dim;						// Colum
					#pragma omp parallel for
					for(int p_L = 0; p_L < N1 ; p_L ++){				// g1 = 0 
						for(int p_R = 0; p_R < N1; p_R ++){				// g2 = 0
							sum += vect_L_[q][ds][(long)r_L*N4 +p_L]*vect_R_[q][ds][(long) r_R*N4 +p_R]*projection_past_[q][ds][(long)(p_L)*dimen_p_[q][ds] +p_R];
						}
					}
					#pragma omp parallel for
					for(int p_L = N1; p_L < N2; p_L ++){				// g1 = 1
						for(int p_R = N1; p_R < N2; p_R ++){			// g2 = 1
							sum += vect_L_[q][ds][(long)r_L*N4 +p_L]*vect_R_[q][ds][(long) r_R*N4 +p_R]*projection_past_[q-1][ds-1][(long)(p_L-N1)*dimen_p_[q-1][ds-1]+p_R-N1]; 
						}
					}
					#pragma omp parallel for
					for(int p_L = N2; p_L < N3; p_L ++){				// g1 = 2
						for(int p_R = N2; p_R < N3; p_R ++){			// g2 = 2
							sum += vect_L_[q][ds][(long)r_L*N4 +p_L]*vect_R_[q][ds][(long) r_R*N4 +p_R]*projection_past_[q-2][ds][(long) (p_L-N2)*dimen_p_[q-2][ds] +p_R-N2];
						}
					}
					#pragma omp parallel for
					for(int p_L = N3; p_L < N4; p_L ++){				// g1 = 3
						for(int p_R = N3; p_R < N4; p_R ++){			// g2 = 3
							sum += vect_L_[q][ds][(long)r_L*N4 +p_L]*vect_R_[q][ds][(long) r_R*N4 +p_R]*projection_past_[q-1][ds+1][(long)(p_L-N3)*dimen_p_[q-1][ds+1]+p_R-N3]; 
						}
					}
				//matrix_to_save[k] = sum;
					projection_write(q,ds,k,sum);
				} // end for k
			} // end if selec == 0
			else{
				projection_alloc_memory(q,ds, (long) 1);
			} 
		} // end if dim > 0
	} // end for ds
} // end for q

// Delete the past projection matrix
for(int q = 0; q < nq-2; q++){
	for(int ds=0; ds < ns-1; ds++){
		int dim  = dimen_p_[q][ds];
		if (dim>0){
			delete[] projection_past_[q][ds];
		}
	}
	delete[] projection_past_[q];
}
delete[] projection_past_;
projection_past_ = NULL;

// Delete the eigenvalues matrix
for(int q = 0; q < nq; q++){
	for(int ds=0; ds < ns; ds++){
		int dim  = dimen_[q][ds];
		if (dim>0){
			delete[] vect_R_[q][ds];
			delete[] vect_L_[q][ds];
		}
	}
	delete[] vect_R_[q];
	delete[] vect_L_[q];
}
delete[] vect_R_;
delete[] vect_L_;
vect_R_ = NULL;
vect_L_ = NULL;

// Printing the projections in the last iteraction
if (N == fN_max()){

std::ofstream file3;
file3.open("Projections.txt");

for(int q=0; q<nq; q++){
	for(int ds = 0; ds < ns; ds++){
		int dim = dimen_[q][ds];
		double proj2 = 0;
		double proj_max = 0;
		int r_R_p_max = 0;
		if((dim> 0)&&(abs(q-N-2) + ds <= 1)&&(ds==0)&&(q-N-2>=0)){
			for(int r_R = 0; r_R < 0*dim+1; r_R++){						// Line /left side
				for(int r_L = 0; r_L < dim; r_L++){					// Colum/ Right side
					long k = r_L*dim + r_R;						
					double sum = projection_read(q,ds,k);
					std::cout <<"Proj["<<q-N-2<<";"<< ds<<"]("<<r_L + 1 << ";" <<r_R + 1 << ") = "; 
					std::cout << sum << "; ";
					file3 << std::setprecision(16) << sum;
					file3 << std::endl;						  
					  proj2 = proj2 + sum*sum;
					  if(abs(sum)>proj_max){ 
						proj_max = abs(sum);
						r_R_p_max = r_L; 
					  }
				}
				file3 << std::endl;
			} // end for k
		 	std::cout << std::endl;
			std::cout <<"sum_n of Proj(0,n)²   = " << proj2 << std::endl;
			std::cout <<"Maximum value of Proj = " << proj_max << "; r =" << (r_R_p_max +1);
			std::cout << std::endl<< std::endl;
		} // end if dim > 0
	} // end for ds
} // end for q

file3.close();

} // end if N==N_max


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
