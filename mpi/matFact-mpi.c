// Serial Implementation - Group 5

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matrix.h"

// Macros
#define RAND01 ((double) random() / (double) RAND_MAX)
#define max(x,y) ( (x) > (y) ? (x) : (y) )

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) \
 (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) \
(((p)*((index)+1)-1)/(n))

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
int **A;
double **L, **R, **RT, **B, **Lsum, **RTsum;
double alpha;

// Main
int main(int argc, char **argv) {

    if(argc != 2) {
        printf("Invalid number of arguments provided: matFact [fileInput]");
        return -1;
    }

    MPI_Init(&argc, &argv);

    int id, nproc;

	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if(!id) {
        // Read all the input and create the necessary structures
        read_input(argv);

        // Fill the matrixes randomly
        random_fill_LR();

        // Get better performance because of cache 
        // by working in the same chunks of memory -> cache hits 
        RT = transpose_matrix(R, nFeatures, nItems);
    }

    // TODO: distribuir matrizes pelos processos (começar por broadcast, todos têm as matrizes inteiras)

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
    Lsum = create_matrix_double(nUsers, nFeatures);
    RTsum = create_matrix_double(nItems, nFeatures);
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
            L[i][j] = RAND01 / (double) nFeatures;

    for (int i = 0; i < nFeatures; i++)
        for (int j = 0; j < nItems; j++)
            R[i][j] = RAND01 / (double) nFeatures;
}

void update() {
    int i, j, n, k;

    // TODO: processo 0 faz broadcast de RT

    for (n = 0; n < nNonZero; n++) {
        i = A[n][0];
        j = A[n][1];

        for (k = 0; k < nFeatures; k++) {
            Lsum[i][k] += alpha * ( 2 * ( A[n][2] - B[i][j] ) * ( -RT[j][k] ) );
            RTsum[j][k] += alpha * ( 2 * ( A[n][2] - B[i][j] ) * ( -L[i][k] ) );
        }
    }

    // TODO: fazer um reduce do Lsum e RTsum para o processo 0

    for (int i = 0; i < max; i++)
        for (int j = 0; j < nFeatures; j++) {
            if (i < nUsers) {
                L[i][j] -= Lsum[i][j];
                Lsum[i][j] = 0;
            }
            if (i < nItems) {
                RT[i][j] -= RTsum[i][j];
                RTsum[i][j] = 0;
            }
        }
}

void loop() {

    for (int i = 0; i < iterations; i++) {
        multiply_non_zeros(L, RT, B, A, nNonZero, nFeatures);

        update();
    }

    // TODO: voltar a agrupar tudo no processo 0

    multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);

}

void print_recomendations() {
    int index = 0;

    for (int user = 0; user < nUsers; user++) {
        double max = -1;
        int recomendation;

        for (int item = 0; item < nItems; item++) {
            if (index < nNonZero && A[index][0] == user && A[index][1] == item) {
                index++;
                continue;
            }

            if (B[user][item] > max) {
                max = B[user][item];
                recomendation = item;
            }
        }
        printf("%d\n", recomendation);
    }
}
