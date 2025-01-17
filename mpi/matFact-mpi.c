// MPI Implementation - Group 5


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
int *A, *recomendations;
double *L, *R, *RT, *B, *Lsum, *RTsum;
double alpha;


// Main
int main(int argc, char **argv) {
    
    double start, end;

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

    } else {
        MPI_Recv(&nElements, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &status);

        A = create_compact_matrix(nElements);

        MPI_Recv(A, nElements * 3, MPI_INT, 0, id, MPI_COMM_WORLD, &status);

        // Receive L
        MPI_Recv(L, BLOCK_SIZE(id, nproc, nUsers) * nFeatures, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &status);
    }

    // Broadcast RT
    MPI_Bcast(RT, nItems * nFeatures, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    loop();

    print_recomendations();

    free_matrix_structures();

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    return 0;
}


void read_input(FILE *file_pointer) {

    int *buffer = create_compact_matrix(((block_size + 2) * nItems) * 3);

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

            if (i == nNonZero - 1) {
                MPI_Send(&elem_number, 1, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
                MPI_Send(buffer, elem_number * 3, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
            }
        } else {
            if (proc_number) {
                MPI_Send(&elem_number, 1, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
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
    
    free_matrix_int(buffer);
    fclose(file_pointer);
}


void create_matrix_structures() {
    if(id) {
        RT = create_matrix_double(nItems, nFeatures);
    } else {
        A = create_compact_matrix(nNonZero);
        R = create_matrix_double(nFeatures, nItems);
    }

    L = create_matrix_double(block_size, nFeatures);
    Lsum = create_matrix_double(block_size, nFeatures);
    RTsum = create_matrix_double(nItems, nFeatures);
    B = create_matrix_double(block_size, nItems); 

}


void free_matrix_structures() {
    free(A);
    free(B);
    free(L);
    free(Lsum);
    free(RT);
    free(RTsum);
    if (!id) {
        free(R);
    }
    free(recomendations);
}


void random_fill_LR() {
    srandom(0);

    double *buffer = create_matrix_double(block_size + 1, nFeatures);

    for (int i = 0; i < block_size; i++)
        for (int j = 0; j < nFeatures; j++)
            L[POS(i,j,nFeatures)] = RAND01 / (double) nFeatures;

    for (int i = block_size, p = 1; i < nUsers; i++) {
        for (int j = 0; j < nFeatures; j++)
            buffer[POS(i - BLOCK_LOW(p, nproc, nUsers),j,nFeatures)] = RAND01 / (double) nFeatures; 

        if (BLOCK_OWNER(i + 1,nproc,nUsers) > p) {
                MPI_Send(buffer, BLOCK_SIZE(p, nproc, nUsers) * nFeatures, MPI_DOUBLE, p, p, MPI_COMM_WORLD);
                p++;
        }
    }

    for (int i = 0; i < nFeatures; i++)
        for (int j = 0; j < nItems; j++)
            R[POS(i,j,nItems)] = RAND01 / (double) nFeatures;
}


void update() {
    int i, j, n, k, i2;

    for (n = 0; n < nElements; n++) {
        i = A[POS(n,0,3)];
        j = A[POS(n,1,3)];
        i2 = i - BLOCK_LOW(id, nproc, nUsers);

        for (k = 0; k < nFeatures; k++) {
            Lsum[POS(i2,k,nFeatures)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i2,j,nItems)] ) * ( -RT[POS(j,k,nFeatures)] ) );
            RTsum[POS(j,k,nFeatures)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i2,j,nItems)] ) * ( -L[POS(i2,k,nFeatures)] ) );
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, RTsum, nFeatures * nItems, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < max(block_size, nItems); i++)
        for (int j = 0; j < nFeatures; j++) {
            if (i < block_size) {
                L[POS(i,j,nFeatures)] -= Lsum[POS(i,j,nFeatures)];
                Lsum[POS(i,j,nFeatures)] = 0;
            }
            if (i < nItems) {
                RT[POS(i,j,nFeatures)] -= RTsum[POS(i,j,nFeatures)];
                RTsum[POS(i,j,nFeatures)] = 0;
            }
        }
}


void loop() {
    
    for (int i = 0; i < iterations; i++) {
        multiply_non_zeros(L, RT, B, A, nElements, nFeatures, nItems, BLOCK_LOW(id, nproc, nUsers));
        update();
    }

    multiply_matrix(L, RT, B, block_size, nItems, nFeatures);
}


void print_recomendations() {
    recomendations = (int *) malloc((block_size + 1) * sizeof(int));

    for (int user = 0, i = 0, index = 0; user < block_size; user++) {
        double max = -1;
        int recomendation;

        for (int item = 0; item < nItems; item++) {
            if (index < nElements && ( A[POS(index,0,3)] - BLOCK_LOW(id, nproc, nUsers) ) == user && A[POS(index,1,3)] == item) {
                index++;
                continue;
            }

            if (B[POS(user,item,nItems)] > max) {
                max = B[POS(user,item,nItems)];
                recomendation = item;
            }
        }
        recomendations[i++] = recomendation;
    }

    if (id) {   
        MPI_Send(recomendations, block_size, MPI_INT, 0, id, MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < block_size; i++)
                printf("%d\n", recomendations[i]);

        MPI_Status status;
        for (int i = 1; i < nproc; i++) {
            int proc_block_size = BLOCK_SIZE(i, nproc, nUsers);
            // Process 0 receives recomendations
            MPI_Recv(recomendations, BLOCK_SIZE(i, nproc, nUsers), MPI_INT, i, i, MPI_COMM_WORLD, &status);
            for (int i = 0; i < proc_block_size; i++)
                printf("%d\n", recomendations[i]);
        }
    }
}
