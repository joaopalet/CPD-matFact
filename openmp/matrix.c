#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

// ---------------------
// MATRICES OPERATIONS
// ---------------------

void multiply_non_zeros(double **L, double **RT, double **B, int **A, int nNonZero, int nFeatures) {
    #pragma omp for
    for (int n = 0; n < nNonZero; n++) {
        int i = A[n][0];
        int j = A[n][1];

        double sum = 0;
        for (int  k= 0; k < nFeatures; k++) {
            sum += L[i][k] * RT[j][k];
        }
        B[i][j] = sum;
    }
}

void multiply_matrix(double **X, double **Y, double **Z, int n, int m, int p) {
    #pragma omp for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0;
            for (int k = 0; k < p; k++)
                sum += X[i][k] * Y[j][k];
            Z[i][j] = sum;
        }
    }
}

double **transpose_matrix(double **M, int n, int m) {
    double **MT = create_matrix_double(m, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            MT[j][i] = M[i][j];

    return MT;
}

// -------------------
// MATRICES CREATION
// -------------------

int **create_compact_matrix(int n) {
    int **M = (int **)malloc(n * sizeof(int *)); 

    for (int i= 0; i < n; i++) 
        M[i] = (int *) calloc(3, sizeof(int));
    return M;
}

double **create_matrix_double(int r, int c) {
    double **M = (double **)malloc(r * sizeof(double *)); 

    for (int i= 0; i < r; i++) 
        M[i] = (double *) calloc(c, sizeof(double));

    return M;
}

void free_matrix_int(int **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void free_matrix_double(double **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// -------------------
// MATRICES PRINTING
// -------------------

void print_matrix_int(int **M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%d ", M[i][j]);
        printf("\n");
    }
    printf("\n");
}

void print_matrix_double(double **M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%f ", M[i][j]);
        printf("\n");
    }
    printf("\n");
}