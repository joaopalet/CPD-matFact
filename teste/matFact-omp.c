// Serial Implementation - Group 5

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

// Macros
#define RAND01 ((double) random() / (double) RAND_MAX)
#define max(x,y) ( (x) > (y) ? (x) : (y) )

#define POS(i,j,columns) ((i)*(columns)+(j))

// Function declarations
void random_fill_LR();
void read_input(char **argv);
void create_matrix_structures();
void free_matrix_structures();
void update();
void loop();
void print_recomendations();

// Global Variables
int iterations, nFeatures, nUsers, nItems, nNonZero, numThreads, max;
int *A;
double *L, *R, *RT, *B;
double **privLsum, **privRTsum;
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

    loop();

    print_recomendations();

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

    max = max(nUsers, nItems);
    
    create_matrix_structures();

    for (int i = 0; i < nNonZero; i++) {
        int n, m;
        double v;

        fscanf(file_pointer, "%d", &n);
        fscanf(file_pointer, "%d", &m);
        fscanf(file_pointer, "%lf", &v);

        A[POS(i,0,3)] = n;
        A[POS(i,1,3)] = m;
        A[POS(i,2,3)] = v;
    }
    
    fclose(file_pointer);
}

void create_matrix_structures() {
    A = create_compact_matrix(nNonZero);
    B = create_matrix_double(nUsers, nItems);
    L = create_matrix_double(nUsers, nFeatures);
    R = create_matrix_double(nFeatures, nItems);
}

void free_matrix_structures() {
    free_matrix_int(A);
    free_matrix_double(B);
    free_matrix_double(L);
    free_matrix_double(R);
    free_matrix_double(RT);
}

void random_fill_LR() {
    srandom(0);

    for (int i = 0; i < nUsers; i++)
        for (int j = 0; j < nFeatures; j++)
            L[POS(i,j,nFeatures)] = RAND01 / (double) nFeatures;

    for (int i = 0; i < nFeatures; i++)
        for (int j = 0; j < nItems; j++)
            R[POS(i,j,nItems)] = RAND01 / (double) nFeatures;
}

void update() {
    int i, j, n, k;

    int myId = omp_get_thread_num();

    #pragma omp for private(i, j, k)
    for (n = 0; n < nNonZero; n++) {
        i = A[POS(n,0,3)];
        j = A[POS(n,1,3)];

        for (k = 0; k < nFeatures; k++) {
            privLsum[myId][POS(i,k,nFeatures)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i,j,nItems)] ) * ( -RT[POS(j,k,nFeatures)] ) );
            privRTsum[myId][POS(j,k,nFeatures)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i,j,nItems)] ) * ( -L[POS(i,k,nFeatures)] ) );
        }
    }

    for (int nt = 0; nt < numThreads; nt++) {
        #pragma omp for
        for (int i = 0; i < max; i++)
            for (int j = 0; j < nFeatures; j++) {
                if (i < nUsers) {
                    L[POS(i,j,nFeatures)] -= privLsum[nt][POS(i,j,nFeatures)];
                    privLsum[nt][POS(i,j,nFeatures)] = 0;
                }
                if (i < nItems) {
                    RT[POS(i,j,nFeatures)] -= privRTsum[nt][POS(i,j,nFeatures)];
                    privRTsum[nt][POS(i,j,nFeatures)] = 0;
                }
        }
    }
}

void loop() {
    #pragma omp parallel if (nFeatures > 50 || nNonZero > 50)
    {
        numThreads = omp_get_num_threads();
        #pragma omp single
        {
            privLsum  = (double **) malloc(numThreads * sizeof(double**));
            privRTsum = (double **) malloc(numThreads * sizeof(double**));
        }

        #pragma omp for
        for (int nt = 0; nt < numThreads; nt++) {
            privLsum[nt]  = create_matrix_double(nUsers, nFeatures);
            privRTsum[nt] = create_matrix_double(nItems, nFeatures);
        }

        for (int i = 0; i < iterations; i++) {
            multiply_non_zeros(L, RT, B, A, nNonZero, nFeatures, nItems);
            update();
        }
        multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);

        #pragma omp for
        for (int nt = 0; nt < numThreads; nt++) {
            free_matrix_double(privLsum[nt]);
            free_matrix_double(privRTsum[nt]);
        }
    }
}

void print_recomendations() {
    int index = 0;

    for (int user = 0; user < nUsers; user++) {
        double max = -1;
        int recomendation;

        for (int item = 0; item < nItems; item++) {
            if (index < nNonZero && A[POS(index,0,3)] == user && A[POS(index,1,3)] == item) {
                index++;
                continue;
            }

            if (B[POS(user,item,nItems)] > max) {
                max = B[POS(user,item,nItems)];
                recomendation = item;
            }
        }
        printf("%d\n", recomendation);
    }
}
