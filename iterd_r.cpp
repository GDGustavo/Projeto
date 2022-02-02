#include "iterd_r.h"

static int **NS_;
static int **NE_;
static int **NN_;
static int **NW_;
static int **dimen_;

//Creating a new function to find the gender of  states {(q, ds, p)} in the primitive Bases;
// gen_ = 0; gen. South;
// gen_ = 1; gen. East;
// gen_ = 2; gen. North;
// gen_ = 3; gen. West;
// gen_ =-1; error;

int genre_R(int q, int ds, int p){

int N1 = 0; 
int N2 = 0; 
int N3 = 0; 
int N4 = 0; 

N1 = NS_[q][ds];     // NS
N2 = N1 + NE_[q][ds];// NS + NE
N3 = N2 + NN_[q][ds];// NS + NE + NN
N4 = N3 + NW_[q][ds];// NS + NE + NN + NW

int gen_ = -1;

if (p <= N1){ gen_ = 0;}
else {
	if ((N1<p)&&(p<= N2)){
		gen_ =1;
	}
	else {
		if ((N2<p)&&(p<= N3)){
			gen_ =2;
		}
		else {
			if ((N3<p)&&(p<= N4)){
				gen_ =3;
			}
			else {
				if (N4 < p){
					gen_ = -1;
				}
			}
		}
	}
}

return gen_;}

void iter0_r(double W){

int nq = 2*(2+0)+1; // charge: {-2,-1, 0, 1, 2} => {0,1,2,3,4}
int ns = (0+2)+1;  // double spin: {0, 1 2} => {0,1,2}
double lamb_ = lamb();
double D_0   = (1-pow(lamb_, (float) -1))*(pow(lamb_,(float) 1/2))/log(lamb_);


// Creating the Basis and Making all matrix elements zero.
NS_ = new int*[nq];
NE_ = new int*[nq];
NN_ = new int*[nq];
NW_ = new int*[nq];
dimen_ = new int*[nq];
for(int i=0; i< nq; i++){
	NS_[i] = new int[ns];
	NE_[i] = new int[ns];
	NN_[i] = new int[ns];
	NW_[i] = new int[ns];
	dimen_[i] = new int[ns];
	for(int j=0; j< ns; j++){
			NS_[i][j] = 0;
			NE_[i][j] = 0;
			NN_[i][j] = 0;
			NW_[i][j] = 0;
			dimen_[i][j] = 0;
	}
}

 NS_[0][0] = 1; // vacuum   |-2 0 1|
 NS_[1][1] = 1; // Cd       |-1 1 1|
 NS_[2][0] = 1; // CdCd     | 0 0 1|
 NE_[1][1] = 1; // f0       |-1 1 2|
 NE_[2][2] = 1; // (f0Cd)_T | 0 2 1|
 NE_[3][1] = 1; // f0CdCd   |+1 1 1|
 NN_[2][0] = 1; // f0f0     | 0 0 2|
 NN_[3][1] = 1; // f0f0Cd   |+1 1 2|
 NN_[4][0] = 1; // fofoCdCd |+2 0 1|
 NW_[2][0] = 1; // (f0Cd)_S | 0 0 3|

int dim_max = 3;  // This Number can be obteined by the ocupation basis. 
int p_max = (dim_max+1)*dim_max/2;

// Building the Hamiltonian for N=0; (Number of elements = p_max²)

double ***H0_ = new double**[nq];
for(int q = 0; q < nq; q++){
	H0_[q] = new double*[ns];
	for(int ds=0; ds < ns; ds++){
		H0_[q][ds] = new double[p_max];
		for(int k=0; k < p_max; k++){
			H0_[q][ds][k] = 0;
		}
	}	
}

H0_[0][0][0] = 0;

H0_[1][1][0] = E_d();
H0_[1][1][1] = sqrt(2)*V_0();
H0_[1][1][2] = 2*W;

H0_[2][0][0] = 2*E_d() + U();
H0_[2][0][1] = 0;
H0_[2][0][2] = 4*W;
H0_[2][0][3] = -2*V_0();
H0_[2][0][4] = -2*V_0();
H0_[2][0][5] = 2*W + E_d();

H0_[2][2][0] = 2*W + E_d();

H0_[3][1][0] = 2*W + 2*E_d()+U();
H0_[3][1][1] = -sqrt(2)*V_0();
H0_[3][1][2] = 4*W + E_d();

H0_[4][0][0] = 4*W + 2*E_d() + U();

eigen2_start(0); 						

// Finding the Bases that diagonalize the Hamiltonian and the eigenvalues for each Base.
// Solving the Hamiltonian for N= 0 using the functions givens;
std::cout << "Diagonalization process for 'Right Side' N = 0 is starting now...." << std::endl<< std::endl;
for (int q=0; q < nq; q++) {
	for (int ds=0; ds < ns; ds++) {
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		dimen_[q][ds] = dim;
		if(dim > 0) {
			std::cout <<'\t'<< "N=0" <<'\t'<<'\t' << "["<< (q - 2) << ";" << ds;
			std::cout << "] Sector;"<<std::endl<<"dim ="<< dim << std::endl;
			eigen2_erg_alloc_memory(q,ds,(long) dim);
			eigen2_vect_alloc_memory(q,ds,(long) dim*dim);
			double *eigen_values = new double[dim];	            		// Changing the size
			double **eigen_vectors = new double*[dim];	  			// Changing the size
			int Nel_max = (1+dim)*dim/2;                                		// Maximum number of elements in LTM
			double *Hamiltonian = new double[Nel_max];             		// Change the size
			for (int k=0; k<Nel_max; k++) {					// Saving the values to diagolize
				Hamiltonian[k]= H0_[q][ds][k]/D_0;
			}
			delete[] H0_[q][ds];
			for (int k=0; k<dim; k++) {					// Making all elements zero
				eigen_values[k] = 0;
				eigen_vectors[k] = new double[dim];
				for(int h=0; h< dim; h++){
					eigen_vectors[k][h] = 0;
				}
			}
			givens(dim, dim , Hamiltonian , eigen_values, eigen_vectors, 0); 		// Solving the H
			delete[] Hamiltonian;
			Hamiltonian = NULL;
			std::cout << "Eigen values: " << std::endl << '\t';
			for (long k=0; k<dim; k++) {
				eigen2_erg_write(q,ds,k,eigen_values[k]);	        			// Saving the eigenvalues 
				std::cout << D_0*eigen_values[k] << "; ";          
			}
			delete[] eigen_values;
			eigen_values = NULL; 
			//std::cout << std::endl << "Eigen vectors matrix:" << std::endl;
			for (int i=0; i<dim; i++) {
				std::cout << '\t';
				for (int j=0; j<dim; j++) { 
					eigen2_vect_write(q,ds,(long)(i*dim +j),eigen_vectors[i][j]);	// Saving the eigenvectors
				//	std::cout << eigen_vectors[i][j]<< ";" <<'\t';
				}
				//std::cout << std::endl;
				delete[] eigen_vectors[i];	
			}
			std::cout << std::endl;
			delete[] eigen_vectors;
			eigen_vectors = NULL;
		}// end if dim>0
	} //end for ds
	delete[] H0_[q]; 
}//end for q
delete[] H0_;
H0_ = NULL;

mel2_start(0);									// Starting the Matrix |Q' S' r'|f_N^+|Q S r|

std::cout << "Diagonalization process already finished for 'Right' N = 0!" << std::endl<< std::endl<< std::endl<< std::endl;

//|Q+1 dS+1 r2|(f0^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS+1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS+1 p2|(f0^+)|Q dS p1|
for(int q=0;q<(nq-1);q++){
	for(int ds=0;ds<(ns-1);ds++){
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		int dim2 = NS_[q+1][ds+1] + NE_[q+1][ds+1] + NN_[q+1][ds+1] + NW_[q+1][ds+1]; 
		if ((dim>0)&&(dim2>0)){
			mel2_ne_alloc_memory(q,ds,(long) dim*dim2);						// 
			for(long k=0;k<dim*dim2;k++){
				double sum = 0;
				double term = 0;
				int r2 = k/dim;      // line
				int r1 = k - r2*dim; // Collum
					for(int p1=0; p1<dim; p1++){
						for(int p2=0; p2<dim2; p2++){
							term = 0; 
							int g1 = genre_R(q,ds,p1+1);    	 	// What is the genre p1? 
							int g2 = genre_R(q+1,ds+1,p2+1);		// What is the genre p2?
							if ((g1==0)&&(g2==1)){term = 1;}
							if ((g1==3)&&(g2==2)){term = sqrt(ds+1)/sqrt(ds+2);}
							double aux2 = eigen2_vect_read(q+1,ds+1, (long) r2*dim2 +p2);
							double aux1 = eigen2_vect_read(q,ds, (long) r1*dim +p1);
							sum = sum + term*aux2*aux1;
						}
					}
				mel2_ne_write(q,ds,k,sum);
			}
		}
	}
}

//|Q+1 dS-1 r2|(f0^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS-1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS-1 p2|(f0^+)|Q dS p1|
for(int q=0;q<(nq-1);q++){
	for(int ds=1;ds<ns;ds++){
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		int dim2 = NS_[q+1][ds-1] + NE_[q+1][ds-1] + NN_[q+1][ds-1] + NW_[q+1][ds-1];
		if ((dim>0)&&(dim2>0)){
			mel2_nw_alloc_memory(q,ds, (long) dim*dim2);						// 
			for(long k=0;k<dim*dim2;k++){
				double sum = 0; 
				double term = 0;
				int r2 = k/dim;      // line
				int r1 = k - r2*dim; // Collum
					for(int p1=0; p1<dim; p1++){
						for(int p2=0; p2<dim2; p2++){
							term = 0; 
							int g1 = genre_R(q,ds,p1+1);     		// What is the genre p1? 
							int g2 = genre_R(q+1,ds-1,p2+1);		// What is the genre p2?
							if ((g1==0)&&(g2==3)){ term = 1;}
							if ((g1==1)&&(g2==2)){term = -sqrt(ds+1)/sqrt(ds);}
							double aux2 = eigen2_vect_read(q+1,ds-1,(long) r2*dim2 +p2);
							double aux1 = eigen2_vect_read(q,ds,(long) r1*dim +p1);
							sum = sum + term*aux2*aux1;	
						}
					}
				mel2_nw_write(q,ds,k,sum);
			}
		}
	}
}

// Projetcions Beteween the "Right" and "Left" sectors. 

std::cout <<"Bases Projection for 'Right Side' N = 0;" << std::endl << std::endl;

projection_start(0);

for(int q=0; q<nq; q++){
	for(int ds = 0; ds < ns; ds++){
		int dim = dimen_[q][ds];
		if(dim> 0){
			projection_alloc_memory(q,ds, (long) dim*dim);
			for(long k = 0; k < dim*dim; k++){
				double sum = 0;
				int r_l = k/dim;						// Line
				int r_r = k - r_l*dim;					// Colum
				for(int p = 0; p < dim; p++){
					double aux_l = eigen_vect_read(q,ds, (long) r_l*dim +p);
					double aux_r = eigen2_vect_read(q,ds,(long) r_r*dim +p);
					sum = sum + aux_l*aux_r;
				}
				if(abs(sum) < 0.0000000001){
					sum = 0;
				}
				projection_write(q,ds,k,sum);
				std::cout <<"Proj["<<q-2<<";"<< ds<<"]("<<r_l + 1 << ";" <<r_r + 1 << ") = " <<'\t'; 
				std::cout << projection_read(q,ds,k) << std::endl;
			} // end for k
		} // end if dim > 0
	} // end for ds
} // end for q


for (int q = 0; q < nq; q++){
	delete[] NS_[q];
	delete[] NE_[q];
	delete[] NN_[q];
	delete[] NW_[q];
	delete[] dimen_[q];
}
delete[] NS_;
delete[] NE_;
delete[] NN_;
delete[] NW_;
delete[] dimen_;
NS_ = NULL;
NE_ = NULL;
NN_ = NULL;
NW_ = NULL;
dimen_ = NULL;

}
