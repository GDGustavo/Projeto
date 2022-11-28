#include <cstddef>

// Start both Matrix.  
void mel_start(int N);

// MEL_NE Functions.  
void mel_ne_alloc_memory(int q, int ds, long nk);
void mel_ne_write(int q, int ds, long k, double value);
double mel_ne_read(int q, int ds, long k);
void mel_ne_delete(int N, int **dimen);

// MEL_NW Functions.
void mel_nw_alloc_memory(int q, int ds, long k);
void mel_nw_write(int q, int ds, long k, double value);
double mel_nw_read(int q, int ds, long k);
void mel_nw_delete(int N, int **dimen);

// EIGEN Functions. 
void eigen_start(int N);
void eigen_delete(int N, int **dimen);
void eigen_vect_delete(int N, int **dimen);
void eigen_erg_alloc_memory(int q, int ds, long nk);
void eigen_vect_alloc_memory(int q, int ds, long nk);
void eigen_erg_write(int q, int ds, long k, double value);
double eigen_erg_read(int q, int ds, long k);
void eigen_vect_write(int q, int ds, long k, double value);
double eigen_vect_read(int q, int ds, long nk);
