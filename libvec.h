#ifdef _INCLUDE_HPUX_SOURCE
extern "C" void vec_$dcopy(double *from, double *to, int *size);
extern "C" double vec_$ddot(double *, double *, int *);
extern "C" void vec_$dadd_vector(double *, double *, int *, double *);
extern "C" void vec_$dmult_add(double *, double *, int *size,
			       double *factor, double *result);
#endif
#ifdef _LINUX
extern "C" void dcopy_(int *size, double *from, int *incf, double *to, int *inct);
#endif /* _LINUX */
