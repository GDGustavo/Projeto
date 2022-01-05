#include <iostream>
#include "mel.h"
#include "mel_2.h"


void projection_start(int N);
void projection_alloc_memory(int q, int ds, long nk);
void projection_write(int q, int ds, long k, double value);
double projection_read(int q, int ds, long k);
void projection_delete(int N, int **dimen);
void save_projection(int N, int **dimen_);
double save_projection_read(int q, int ds, long k);
void save_projection_delete(int N, int **dimen_);
