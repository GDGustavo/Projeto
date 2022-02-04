#include "iterd.h"

static int **NS_;
static int **NE_;
static int **NN_;
static int **NW_;

//Creating a new function to find the gender of  states {(q, ds, p)} in the primitive Bases;
// gen_ = 0; gen. South;
// gen_ = 1; gen. East;
// gen_ = 2; gen. North;
// gen_ = 3; gen. West;
// gen_ =-1; error;

int genre(int q, int ds, int p){

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

void iter0_l(double W, int **dimen_){

int nq = 2*(0+2)+1; // charge: {-2,-1, 0, 1, 2} => {0,1,2,3,4}
int ns = (0+2)+1;  // double spin: {0, 1 2} => {0,1,2}

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

// Completing the Bases for N= 0
// In computation basis the endress in the memory for each basis = (p-1)
//Numeric Basis            Analytical Basis
//{|[Qn][Sn]|}             {|Q 2S p|} => Q = Qn - 2 and S = Sn/2;
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
// p_max(Q,S) = NS(Q,S) + NE(Q,S) + NN(Q,S) + NW(Q,S); Total Number of primitive Basis for each (Q,S) sector;
// W = {W_1 or W_2};

// H0(-2, 0) = [0];                                                                    {|-2,0,p|}, p=1     (NElem=1)
// H0(-1, 1) = [Ed, raiz(2)V ; raiz(2)V, 2W];                                          {|-1,1,p|}, p=1,2   (NElem=4)
// H0( 0, 0) = [(2Ed+U0), 0, -2V; 0, 4W, -2V; -2V, -2V, (2W+Ed)]; 		      {| 0,0,p|}, p=1,2,3 (NEelm=9)
// H0( 0, 2) = [2W + Ed];                                                              {| 0,2,p|}, p=1     (NEelm=1)
// H0(+1, 1) = [(2W+2Ed+U), -raiz(2)V; -raiz(2)V, (4W+Ed)];                            {|+1,1,p|}, p=1,2   (NEelm=4)
// H0(+2, 0) = [4W + 2Ed + U];                                                         {|+2,0,p|}, p=1     (NEelm=1)

// Hamiltonian in the Linear Lower-Triangular form (Number of elements = 1 + 2 + ... +p_max = (1+p_max)p_max/2 )

// H0(-2, 0) = [0];                                                                    {|-2,0,p|}, p=1     (NElem=1)
// H0(-1, 1) = [Ed, raiz(2)V, 2W];                                                     {|-1,1,p|}, p=1,2   (NElem=3)
// H0( 0, 0) = [(2Ed+U), 0, 4W, -2V, -2V, (2W+Ed)];                                    {| 0,0,p|}, p=1,2,3 (NElem=6)
// H0( 0, 2) = [2W + Ed];                                                              {| 0,2,p|}, p=1     (NElem=1)
// H0(+1, 1) = [(2W+2Ed+U), -raiz(2)V, (4W+Ed)];                                       {|+1,1,p|}, p=1,2   (NElem=3)
// H0(+2, 0) = [4W + 2Ed + U];                                                         {|+2,0,p|}, p=1     (NElem=1)

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

double lamb_ = lamb();									// Lambda
double D_0   = (1-pow(lamb_, (float) -1))*(pow(lamb_,(float) 1/2))/log(lamb_);			// D_0
double E_f   = (double) E_uv();								// Fundamental Energy 

eigen_start(0); 						

// Finding the Bases that diagonalize the Hamiltonian and the eigenvalues for each Base.
// Solving the Hamiltonian for N= 0 using the functions givens;
std::cout << "Diagonalization process for  'Left Side' N = 0 is starting now...." << std::endl<< std::endl;
for (int q=0; q < nq; q++) {
	for (int ds=0; ds < ns; ds++) {
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		dimen_[q][ds] = dim;							// Saving the Dimension for each (Q,dS)
		if(dim > 0) {
			eigen_erg_alloc_memory(q,ds,(long) dim);
			eigen_vect_alloc_memory(q,ds,(long) dim*dim);
			double *eigen_values = new double[dim];	            		// Matrix to calculate E-Values
			double **eigen_vectors = new double*[dim];	  			// Matrix to calculate E-Vectors
			for (int k=0; k<dim; k++) {
				eigen_vectors[k] = new double[dim];
			}
			int Nel_max = (1+dim)*dim/2;                                		// Maximum number of elements in LTM
			double *Hamiltonian = new double[Nel_max];             		// Change the size
			for (int k=0; k<Nel_max; k++) {					// Saving the values to diagolize
				Hamiltonian[k]= H0_[q][ds][k]/D_0;				// Att, scaled hamiltonian.
			}
			delete[] H0_[q][ds];
			
			int ret = givens(dim, dim , Hamiltonian , eigen_values, eigen_vectors, 1); 	// Solving the H
			
			delete[] Hamiltonian;
			Hamiltonian = NULL;

			for (long k=0; k<dim; k++) {
				eigen_erg_write(q,ds,k,eigen_values[k]);	        			// Saving the eigenvalues
				if (E_f > eigen_values[k]){
					E_f = eigen_values[k];
				}
			}

			delete[] eigen_values;
			eigen_values = NULL; 
			
			for (int i=0; i<dim; i++) {
				for (int j=0; j<dim; j++) {
					long k =i*dim +j;
					eigen_vect_write(q,ds,k,eigen_vectors[i][j]);		// Saving the eigenvectors
				}
				delete[] eigen_vectors[i];	
			}
			delete[] eigen_vectors;
			eigen_vectors = NULL;
		}// end if dim>0
	} //end for ds
	delete[] H0_[q]; 
}//end for q
delete[] H0_;
H0_ = NULL;

// Printing the eigen energies and eigen vectors
std::cout << std::endl << "Fundamental Energy (Non Escaled) for N = 0:" << '\t' << D_0*E_f << std::endl<< std::endl;
for (int q=0; q < nq; q++) {
	for (int ds=0; ds < ns; ds++) {
		int dim = dimen_[q][ds];
		if(dim > 0) {
			std::cout << "["<< (q - 2) << ";" << ds;
			std::cout << "] Sector"<<'\t'<<"dim ="<< dim << std::endl;
			std::cout << "Eigen values (Non Scaled): " << std::endl;

			for (long k=0; k<dim; k++) {
				double Energy = eigen_erg_read(q,ds,k) - E_f;
				eigen_erg_write(q,ds,k,Energy);	        	
				std::cout << D_0*Energy << ";";          
			}

			std::cout << std::endl << "Eigen vectors matrix:" << std::endl;
			for (int i=0; i<dim; i++) {
				std::cout << '\t';
				for (int j=0; j<dim; j++) {
					long k =i*dim +j;
					double vector = eigen_vect_read(q,ds,k);
					std::cout << vector << ";" <<'\t';
				}
				std::cout << std::endl;
			}
			std::cout << std::endl<< std::endl<< std::endl;
		}// end if dim>0
	} //end for ds
}//end for q



mel_start(0);									// Starting the Matrix |Q' S' r'|f_N^+|Q S r|

// Matrix problem the Vector V_d is save in the LINE from matrix U: [V_d(i) = sum_j u(i,j)V_0(j);] ok. 
// Non null elements off Matrix |Q' dS' p'|(f0^+)|Q dS p| for N=0
//             Non null Elements in Analitical Basis
//        |-1  1 2 E|f0+|-2 0 1 S| = 1             direction ne
//        | 0  2 1 E|f0+|-1 1 1 S| = 1             direction ne
//        | 0  0 3 W|f0+|-1 1 1 S| = 1             direction nw
//        | 0  0 2 N|f0+|-1 1 2 E| = -sqrt(2)      direction nw
//        |+1  1 1 E|f0+| 0 0 1 S| = 1             direction ne
//        |+1  1 2 N|f0+| 0 0 3 W| = sqrt(1/2)     direction ne
//        |+1  1 2 N|f0+| 0 2 1 E| = -sqrt(3/2)    direction nw
//        |+2  0 1 N|f0+|+1 1 1 E| = -sqrt(2)      direction nw
// I don't need to save this in a matrix.
// I just need to create a function to simulate this matrix, since we already know how to get in this. OK.

// Saving the matrix elements 
// Mel(Q,S,Nk) = Mel(Q,S,i,j), where NR = N1*N2 and Nk in {0,1,2,.....,N1*N2-1}
// {0,1,2,....,N1*N2-1} = {[0,1, ... N2-1], [N2,N2+1,....,2*N2-1], .....} 
// If i have a position Nk beteween 0 and N2*N1-1 then,
// line i = integer_less(Nk/N2); 0<i<N1-1; 
// Column j = Nk - N2*integer_less(Nk/N2);   0<j<N2-1;  

//|Q+1 dS+1 r2|(f0^+)|Q dS r1| = sum_{p2,p1} U*(Q+1;dS+1)[r2,p2]U(Q;S)[r1,p1]*|Q+1 dS+1 p2|(f0^+)|Q dS p1|
for(int q=0;q<(nq-1);q++){
	for(int ds=0;ds<(ns-1);ds++){
		int dim = NS_[q][ds] + NE_[q][ds] + NN_[q][ds] + NW_[q][ds];
		int dim2 = NS_[q+1][ds+1] + NE_[q+1][ds+1] + NN_[q+1][ds+1] + NW_[q+1][ds+1]; 
		if ((dim>0)&&(dim2>0)){
			mel_ne_alloc_memory(q,ds,(long) dim*dim2);						// 
			for(long k=0; k < (long) dim*dim2; k++){
				double sum = 0;
				double term = 0;
				int r2 = k/dim;      // line
				int r1 = k - r2*dim; // Collum
					for(int p1=0; p1<dim; p1++){
						for(int p2=0; p2<dim2; p2++){
							term = 0; 
							int g1 = genre(q,ds,p1+1);    	 	// What is the genre p1? 
							int g2 = genre(q+1,ds+1,p2+1);		// What is the genre p2?
							if ((g1==0)&&(g2==1)){term = 1;}
							if ((g1==3)&&(g2==2)){term = sqrt(ds+1)/sqrt(ds+2);}
							  double aux2 = eigen_vect_read(q+1,ds+1,(long) r2*dim2 +p2);
							  double aux1 = eigen_vect_read(q,ds, (long) r1*dim +p1);
							sum = sum + term*aux2*aux1;
						}
					}
				mel_ne_write(q,ds,k,sum);
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
			mel_nw_alloc_memory(q,ds,(long) dim*dim2);						// 
			for(long k=0;k< (long) dim*dim2;k++){
				double sum = 0; 
				double term = 0;
				int r2 = k/dim;      // line
				int r1 = k - r2*dim; // Collum
					for(int p1=0; p1<dim; p1++){
						for(int p2=0; p2<dim2; p2++){
							term = 0; 
							int g1 = genre(q,ds,p1+1);     		// What is the genre p1? 
							int g2 = genre(q+1,ds-1,p2+1);		// What is the genre p2?
							if ((g1==0)&&(g2==3)){ term = 1;}
							if ((g1==1)&&(g2==2)){term = -sqrt(ds+1)/sqrt(ds);}
							  double aux2 = eigen_vect_read(q+1,ds-1,(long) r2*dim2 +p2);
							  double aux1 = eigen_vect_read(q,ds,(long) r1*dim +p1);
							sum = sum + term*aux2*aux1;	
						}
					}
				mel_nw_write(q,ds,k,sum);
			}
		}
	}
}

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

// End Function ITER0(); OK. Working very well!

/*
double t0 = 0.619806;
double t1 = 0.271711;

Eigen::MatrixXd H;
H = Eigen::MatrixXd::Ones(4,4);
for(int i=0; i< 4; i++){
	for (int j = 0; j<4; j++){ 
			H(i,j) = 0;
	}
}

H(0,0) = E_d();
H(0,1) = sqrt(2)*V_0();

H(1,2) = t0;
H(1,1) = 2*W_1();
H(1,0) = sqrt(2)*V_0();

H(2,1) = t0;
H(2,3) = t1; 

H(3,2) = t1;

std::cout << H << std::endl;
Eigen::EigenSolver<Eigen::MatrixXd> s(H);
std::cout << s.eigenvalues() << std::endl;
std::cout << s.eigenvectors() << std::endl;
//*/


}
