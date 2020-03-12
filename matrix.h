#ifndef MATRIX_LIB
#define MATRIX_LIB

// ---------------------
// MATRICES OPERATIONS
// ---------------------

void multiply_matrix(double **X, double **Y, double **Z, int n, int m, int p);
double **transpose_matrix(double **M, int n, int m);

// -------------------
// MATRICES CREATION
// -------------------

int **create_matrix_int(int r, int c);
double **create_matrix_double(int r, int c);

void free_matrix_int(int **matrix, int rows);
void free_matrix_double(double **matrix, int rows);

// -------------------
// MATRICES PRINTING
// -------------------

void print_matrix_int(int **M, int r, int c);
void print_matrix_double(double **M, int r, int c);

#endif