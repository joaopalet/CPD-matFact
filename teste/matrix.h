#ifndef MATRIX_LIB
#define MATRIX_LIB

// ---------------------
// MATRICES OPERATIONS
// ---------------------

void multiply_non_zeros(double *L, double *RT, double *B, int *A, int nNonZero, int nFeatures, int nItems);
void copy_matrix(double *m1, double *m2, int n, int m);
void multiply_matrix(double *X, double *Y, double *Z, int n, int m, int p);
double *transpose_matrix(double *M, int n, int m);

// -------------------
// MATRICES CREATION
// -------------------

int *create_compact_matrix(int n);
double *create_matrix_double(int r, int c);

void free_matrix_int(int *matrix);
void free_matrix_double(double *matrix);

// -------------------
// MATRICES PRINTING
// -------------------

void print_matrix_int(int *M, int r, int c);
void print_matrix_double(double *M, int r, int c);

#endif