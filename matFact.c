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

// Global Variables
int iterations, nFeatures, nUsers, nItems, nNonZero;
int **A;
int **nNonZeroPositions;
double **L, **R, **RT, **B, **Lnew, **RTnew, **mAux;
double alpha;


void update(){
    int i, j, n, k;

    copy_matrix(L, Lnew, nUsers, nFeatures);
    copy_matrix(RT, RTnew, nItems, nFeatures);

    for (n = 0; n < nNonZero; n++) {
        i = nNonZeroPositions[n][0];
        j = nNonZeroPositions[n][1];

        for (k = 0; k < nFeatures; k++) {
            double lSum = 0, rSum = 0;

            for (int jSum = 0; jSum < nItems; jSum++)
                lSum += 2 * ( A[i][j] - B[i][j] ) * (-RT[j][k] );
                
            Lnew[i][k] = L[i][k] - (alpha * lSum);

            for (int iSum = 0; iSum < nUsers; iSum++)
               rSum += 2 * ( A[i][j] - B[i][j] ) * (-L[i][k] );

            RTnew[j][k] = RT[j][k] - (alpha * rSum);
        }

        mAux = Lnew;
        Lnew = L;
        L = mAux;

        mAux = RTnew;
        RTnew = RT;
        RT = mAux;
    }
}

void loop() {
    for (int i = 0; i < iterations; i++) {
        multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);
        update();
    }

    multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);
}

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

    print_matrix_int(A, nUsers, nItems);
    print_matrix_double(L, nUsers, nFeatures);
    print_matrix_double(R, nFeatures, nItems);
    print_matrix_double(RT, nItems, nFeatures);
    print_matrix_double(B, nUsers, nItems);
    print_matrix_int(nNonZeroPositions, nNonZero, 2);

    loop();

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

        A[n][m] = v;
        nNonZeroPositions[i][0] = n;
        nNonZeroPositions[i][1] = m;
    }
    
    fclose(file_pointer);
}

void create_matrix_structures() {
    A = create_matrix_int(nUsers, nItems);
    B = create_matrix_double(nUsers, nItems);
    L = create_matrix_double(nUsers, nFeatures);
    R = create_matrix_double(nFeatures, nItems);
    Lnew = create_matrix_double(nUsers, nFeatures);
    RTnew = create_matrix_double(nItems, nFeatures);
    nNonZeroPositions = create_matrix_int(nNonZero, 2);
}

void free_matrix_structures() {
    free_matrix_int(A, nUsers);
    free_matrix_double(B, nUsers);
    free_matrix_double(L, nUsers);
    free_matrix_double(R, nFeatures);
    free_matrix_double(Lnew, nUsers);
    free_matrix_double(RTnew, nFeatures);
    free_matrix_int(nNonZeroPositions, nNonZero);
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