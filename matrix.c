#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

// ---------------------
// MATRICES OPERATIONS
// ---------------------

void copy_matrix(double **m1, double **m2, int n, int m) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            m2[i][j] = m1[i][j];
}

void multiply_matrix(double **X, double **Y, double **Z, int n, int m, int p) {
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

int **create_matrix_int(int r, int c) {
    int **M = (int **)malloc(r * sizeof(int *)); 

    for (int i= 0; i < r; i++) 
        M[i] = (int *) calloc(c, sizeof(int));
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