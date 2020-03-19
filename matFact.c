// Serial Implementation - Group 5
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

// Macros
#define RAND01 ((double) random() / (double) RAND_MAX)

// Function declarations
void random_fill_LR();
void read_input(char **argv);
void create_matrix_structures();
void free_matrix_structures();
void update();
void loop();

// Global Variables
int iterations, nFeatures, nUsers, nItems, nNonZero;
int **A;
double **L, **R, **RT, **B, **Lnew, **RTnew, **mAux;
double alpha;



// Main
int main(int argc, char **argv) {
    

    if(argc != 2) {
        printf("Invalid number of arguments provided: matFact [fileInput]");
        return -1;
    }

    // Read all the input and create the necessary structures
    read_input(argv);

    // Fill the matrixes randomly
    random_fill_LR();

    // Get better performance because of cache 
    // by working in the same chunks of memory -> cache hits 
    RT = transpose_matrix(R, nFeatures, nItems);

    multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);

    print_matrix_int(A, nNonZero, 3);
    print_matrix_double(L, nUsers, nFeatures);
    print_matrix_double(R, nFeatures, nItems);
    print_matrix_double(RT, nItems, nFeatures);
    print_matrix_double(B, nUsers, nItems);

    loop();
    multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);

    print_matrix_double(L, nUsers, nFeatures);
    print_matrix_double(RT, nItems, nFeatures);
    print_matrix_double(B, nUsers, nItems);

    free_matrix_structures();
    
    return 0;
}


void read_input(char **argv) {
    FILE *file_pointer;
    file_pointer = fopen(argv[1], "r");
    fscanf(file_pointer, "%d", &iterations);
    fscanf(file_pointer, "%lf", &alpha);
    fscanf(file_pointer, "%d", &nFeatures);
    fscanf(file_pointer, "%d", &nUsers);
    fscanf(file_pointer, "%d", &nItems);
    fscanf(file_pointer, "%d", &nNonZero);
    
    create_matrix_structures();

    for (int i = 0; i < nNonZero; i++) {
        int n, m;
        double v;

        fscanf(file_pointer, "%d", &n);
        fscanf(file_pointer, "%d", &m);
        fscanf(file_pointer, "%lf", &v);

        A[i][0] = n;
        A[i][1] = m;
        A[i][2] = v;
    }
    
    fclose(file_pointer);
}

void create_matrix_structures() {
    A = create_compact_matrix(nNonZero);
    B = create_matrix_double(nUsers, nItems);
    L = create_matrix_double(nUsers, nFeatures);
    R = create_matrix_double(nFeatures, nItems);
    Lnew = create_matrix_double(nUsers, nFeatures);
    RTnew = create_matrix_double(nItems, nFeatures);
}

void free_matrix_structures() {
    free_matrix_int(A, nNonZero);
    free_matrix_double(B, nUsers);
    free_matrix_double(L, nUsers);
    free_matrix_double(R, nFeatures);
    free_matrix_double(Lnew, nUsers);
    free_matrix_double(RTnew, nFeatures);
    free_matrix_double(RT, nItems);
}

void random_fill_LR() {
    srandom(0);

    for (int i = 0; i < nUsers; i++)
        for (int j = 0; j < nFeatures; j++)
            L[i][j] = RAND01 / (double) nFeatures;

    for (int i = 0; i < nFeatures; i++)
        for (int j = 0; j < nItems; j++)
            R[i][j] = RAND01 / (double) nFeatures;
}

void update(){
    int i, j, n, k;

    copy_matrix(L, Lnew, nUsers, nFeatures);
    copy_matrix(RT, RTnew, nItems, nFeatures);

    for (n = 0; n < nNonZero; n++) {
        i = A[n][0];
        j = A[n][1];

        for (k = 0; k < nFeatures; k++) {
            Lnew[i][k] -= alpha * ( 2 * ( A[n][2] - B[i][j] ) * ( -RT[j][k] ) );
            RTnew[j][k] -= alpha * ( 2 * ( A[n][2] - B[i][j] ) * ( -L[i][k] ) );
        }
    }

    mAux = Lnew;    Lnew = L;       L = mAux;
    mAux = RTnew;   RTnew = RT;     RT = mAux;
}

void loop() {
    for (int i = 0; i < iterations; i++) {
        update();
        //multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);
        multiply_non_zeros(L, RT, B, A, nNonZero, nFeatures);
    }
}