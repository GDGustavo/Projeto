#ifndef _GIV_
#define _GIV_
/*
upon calling:
dim = matrix dimension
for nrootx > 0, compute nroot=nrootx eigenvectors associated with
lowest nrootx eigenvalues; return nroot
for nrootx < 0, compute no eigenvectors; return 0
for nrootx == 0, compute eigenvalues associated with eigenvectors below the value stored in root[0]

a[] = matrix to diagonalize, stored in lower-triangular form:
a[0]
a[1] a[2]
a[3] a[4] a[5] ... ;
the diagonalization destroys a;
root[] allocated space for dim eigenvalues;
if nrootx < 0, root[0] is cut-off eigenvalue
vect[][] allocated space for dim eigenvectors of dimension dim
check == 1 to check diagonalization

on return:
root[] stores eigenvalues in ascending order:
root[0], root[1], ..., root[dim]
vect[][] stores eigenvectors:
vect[0][0] vect[0][1] ... vect[0][dim]
vect[1][1] vect[1][2] ... vect[1][dim]
...
vect[nroot][0] vect[nroot][1] ...vect[nroot][dim]

returns nroot
*/
extern int givens(int n, int nrootx, double *a, 
		  double *root, double *vect[], int check =  0); 
#endif
