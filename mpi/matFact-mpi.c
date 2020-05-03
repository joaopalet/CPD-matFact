// Serial Implementation - Group 5

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matrix.h"

// Macros
#define RAND01 ((double) random() / (double) RAND_MAX)
#define max(x,y) ( (x) > (y) ? (x) : (y) )
#define POS(i,j,columns) ((i)*(columns)+(j))

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) \
 (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) \
(((p)*((index)+1)-1)/(n))

// Function declarations
void random_fill_LR();
void read_input(FILE *file_pointer);
void create_matrix_structures();
void free_matrix_structures();
void update();
void loop();
void print_recomendations();

// Global Variables
int iterations, nFeatures, nUsers, nItems, nNonZero, numThreads, max, block_size, id, nproc, nElements = 0;
int *A;
double *L, *R, *RT, *B, *Lsum, *RTsum, *RTsumcopy;
double alpha;

// Main
int main(int argc, char **argv) {

    if(argc != 2) {
        printf("Invalid number of arguments provided: matFact [fileInput]");
        return -1;
    }

    MPI_Init(&argc, &argv);
    MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    FILE *file_pointer;
    file_pointer = fopen(argv[1], "r");
    fscanf(file_pointer, "%d", &iterations);
    fscanf(file_pointer, "%lf", &alpha);
    fscanf(file_pointer, "%d", &nFeatures);
    fscanf(file_pointer, "%d", &nUsers);
    fscanf(file_pointer, "%d", &nItems);
    fscanf(file_pointer, "%d", &nNonZero);
    if (id) fclose(file_pointer);

    block_size = BLOCK_SIZE(id, nproc, nUsers);

    create_matrix_structures();

    if(!id) {
        // Read all the input and create the necessary structures
        read_input(file_pointer);

        // Fill the matrixes randomly
        random_fill_LR();

        // Get better performance because of cache 
        // by working in the same chunks of memory -> cache hits 
        RT = transpose_matrix(R, nFeatures, nItems);

        // Scatter L
        for (int i = 1; i < nproc; i++)
            MPI_Send(&L[POS(BLOCK_LOW(i, nproc, nUsers), 0, nFeatures)], 
                     BLOCK_SIZE(i, nproc, nUsers) * nFeatures, 
                     MPI_DOUBLE, i, i, MPI_COMM_WORLD);
    
    } else {
        MPI_Recv(&nElements, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &status);
        printf("I am process %d. Received message 1.\n", id);
        MPI_Recv(A, nElements * 3, MPI_INT, 0, id, MPI_COMM_WORLD, &status);
        printf("I am process %d. Received message 2.\n", id);

        // Receive L
        MPI_Recv(L, BLOCK_SIZE(id, nproc, nUsers) * nFeatures, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &status);
    }

    // Broadcast RT
    MPI_Bcast(RT, nItems * nFeatures, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("I am processor %d. This is my A:\n", id);
    print_matrix_int(A, nElements, 3);
    
    printf("I am processor %d. This is my L:\n", id);
    print_matrix_double(L, BLOCK_SIZE(id, nproc, nUsers), nFeatures);

    printf("I am processor %d. This is my R:\n", id);
    print_matrix_double(RT, nItems, nFeatures);

    printf("\n\n\n");

    MPI_Barrier(MPI_COMM_WORLD);

    loop();

    // print_recomendations();

    free_matrix_structures();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}



void read_input(FILE *file_pointer) {

    int *buffer = create_compact_matrix(block_size);

    int proc_number = 0;
    int elem_number = 0;
    int row_number = 0;
    int prev_n = 0;

    for (int i = 0; i < nNonZero; i++) {
        int n, m;
        double v;

        int high = BLOCK_HIGH(proc_number, nproc, nUsers);

        fscanf(file_pointer, "%d", &n);
        fscanf(file_pointer, "%d", &m);
        fscanf(file_pointer, "%lf", &v);
        
        if (n <= high) {
            if (proc_number == 0) {
                A[POS(i,0,3)] = n;
                A[POS(i,1,3)] = m;
                A[POS(i,2,3)] = v;

                nElements++;
            } else {
                buffer[POS(elem_number,0,3)] = n;
                buffer[POS(elem_number,1,3)] = m;
                buffer[POS(elem_number,2,3)] = v;

                elem_number++;
            }

            if (i == nNonZero - 1)
            {
                printf("Sending non zeros to process %d (message 1)\n", proc_number);
                MPI_Send(&elem_number, 1, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
                printf("Sending non zeros to process %d (message 2)\n", proc_number);
                MPI_Send(buffer, elem_number * 3, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
            }
        } else {
            if (proc_number) {
                printf("Sending non zeros to process %d (message 1)\n", proc_number);
                MPI_Send(&elem_number, 1, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
                printf("Sending non zeros to process %d (message 2)\n", proc_number);
                MPI_Send(buffer, elem_number * 3, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
            
                elem_number = 0;
            }

            buffer[POS(elem_number,0,3)] = n;
            buffer[POS(elem_number,1,3)] = m;
            buffer[POS(elem_number,2,3)] = v;
            elem_number++;

            proc_number++;
        }
    }
    
    fclose(file_pointer);
}

void create_matrix_structures() {
    A = create_compact_matrix(block_size);

    if(id) {
        L = create_matrix_double(block_size, nFeatures);
        RT = create_matrix_double(nItems, nFeatures);
        B = create_matrix_double(block_size, nItems); 
    } else {
        B = create_matrix_double(nUsers, nItems); 
        L = create_matrix_double(nUsers, nFeatures);
        R = create_matrix_double(nFeatures, nItems);
        RTsumcopy = create_matrix_double(nItems, nFeatures);
    }

    Lsum = create_matrix_double(block_size, nFeatures);
    RTsum = create_matrix_double(nItems, nFeatures);
}

void free_matrix_structures() {
    free_matrix_int(A);
    free_matrix_double(B);
    free_matrix_double(L);
    free_matrix_double(Lsum);
    free_matrix_double(RT);
    free_matrix_double(RTsum);
    if (!id)
    {
        free_matrix_double(R);
        free_matrix_double(RTsumcopy);
    }
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

    for (n = 0; n < nElements; n++) {
        i = A[POS(n,0,3)];
        j = A[POS(n,1,3)];
        

        for (k = 0; k < nFeatures; k++) {
            Lsum[POS(i,k,nFeatures)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i,j,nItems)] ) * ( -RT[POS(j,k,nFeatures)] ) );
            RTsum[POS(j,k,nFeatures)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i,j,nItems)] ) * ( -L[POS(i,k,nFeatures)] ) );
        }
    }

    // if (!id)
    // {
    //     print_matrix_double(RTsum, nItems, nFeatures);
    // }

    // Fazer um reduce do Lsum e RTsum para o processo 0
    MPI_Reduce(RTsum, RTsumcopy, nFeatures * nItems, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!id)
    {
        double * aux = RTsumcopy;
        RTsumcopy = RTsum;
        RTsum = aux; 
    }

    MPI_Bcast(RTsum, nItems * nFeatures, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < max(block_size, nItems); i++)
        for (int j = 0; j < nFeatures; j++) {
            if (i < block_size) {
                L[POS(i,j,nFeatures)] -= Lsum[POS(i,j,nFeatures)];
                Lsum[POS(i,j,nFeatures)] = 0;
            }
            if (i < nItems) {
                RT[POS(i,j,nFeatures)] -= RTsum[POS(i,j,nFeatures)];
                RTsum[POS(i,j,nFeatures)] = 0;
                if (!id)
                    RTsumcopy[POS(i,j,nFeatures)] = 0;
            }
        }
}

void loop() {
    
    for (int i = 0; i < iterations; i++) {
        multiply_non_zeros(L, RT, B, A, nElements, nFeatures, nItems);

        update();
    }
    

    if (id)
    {   
        MPI_Send(L, BLOCK_SIZE(id, nproc, nUsers) * nFeatures, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
    } else {
        // print_matrix_double(L, nUsers, nFeatures); // ANTES

        MPI_Status status;
        for (int i = 1; i < nproc; i++) {
            // Receive L
            // MPI_recv(buf, 32, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&L[POS(BLOCK_LOW(i, nproc, nUsers), 0, nFeatures)], BLOCK_SIZE(i, nproc, nUsers) * nFeatures, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        }
        
        // print_matrix_double(L, nUsers, nFeatures); // DEPOIS
        // print_matrix_double(RT, nItems, nFeatures);

        multiply_matrix(L, RT, B, nUsers, nItems, nFeatures);
        
        print_matrix_double(B, nUsers, nItems);
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
