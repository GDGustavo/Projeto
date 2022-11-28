#include <cstddef>

// Start both Matrix.  
void mel2_start(int N);

// MEL_NE Functions.  
void mel2_ne_alloc_memory(int q, int ds, long nk);
void mel2_ne_write(int q, int ds, long k, double value);
double mel2_ne_read(int q, int ds, long k);
void mel2_ne_delete(int N, int **dimen);

// MEL_NW Functions.
void mel2_nw_alloc_memory(int q, int ds, long nk);
void mel2_nw_write(int q, int ds, long k, double value);
double mel2_nw_read(int q, int ds, long k);
void mel2_nw_delete(int N, int **dimen);

// EIGEN Functions. 
void eigen2_start(int N);
void eigen2_delete(int N, int **dimen);
void eigen2_vect_delete(int N, int **dimen);
void eigen2_erg_alloc_memory(int q, int ds, long nk);
void eigen2_vect_alloc_memory(int q, int ds, long nk);
void eigen2_erg_write(int q, int ds, long k, double value);
double eigen2_erg_read(int q, int ds, long k);
void eigen2_vect_write(int q, int ds, long k, double value);
double eigen2_vect_read(int q, int ds, long k);
