#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

#define POS(i,j,columns) ((i)*(columns)+(j))

// ---------------------
// MATRICES OPERATIONS
// ---------------------

void multiply_non_zeros(double *L, double *RT, double *B, int *A, int nElements, int nFeatures, int nItems, int block_low) {
    for (int n = 0; n < nElements; n++) {
        int i = A[POS(n,0,3)];
        int j = A[POS(n,1,3)];
        int i2 = i - block_low;

        double sum = 0;
        for (int  k= 0; k < nFeatures; k++) {
            sum += L[POS(i2,k,nFeatures)] * RT[POS(j,k,nFeatures)];
        }
        B[POS(i2,j,nItems)] = sum;
    }
}

void multiply_matrix(double *X, double *Y, double *Z, int n, int m, int p) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0;
            for (int k = 0; k < p; k++)
                sum += X[POS(i,k,p)] * Y[POS(j,k,p)];
            Z[POS(i,j,m)] = sum;
        }
    }
}

double *transpose_matrix(double *M, int n, int m) {
    double *MT = create_matrix_double(m, n);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
    
            MT[POS(j,i,n)] = M[POS(i,j,m)];
        }

    return MT;
}

// -------------------
// MATRICES CREATION
// -------------------

int *create_compact_matrix(int n) {
    int *M = (int *)calloc(3 * n, sizeof(int));
    return M;
}

double *create_matrix_double(int r, int c) {
    double *M = (double *)calloc(c * r, sizeof(double)); 
    return M;
}

void free_matrix_int(int *matrix) {
    free(matrix);
}

void free_matrix_double(double *matrix) {
    free(matrix);
}

// -------------------
// MATRICES PRINTING
// -------------------

void print_matrix_int(int *M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%d ", M[POS(i,j,c)]);
        printf("\n");
    }
    printf("\n");
}

void print_matrix_double(double *M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%f ", M[POS(i,j,c)]);
        printf("\n");
    }
    printf("\n");
}