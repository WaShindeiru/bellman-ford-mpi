#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>

#include "mpi.h"

#define INF INT_MAX / 2

int N; //number of vertices
int *mat; // the adjacency matrix

void abort_with_error_message(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    abort();
}

int convert_dimension_2D_1D(int x, int y, int n) {
    return x * n + y;
}

int print_result(bool has_negative_cycle, int *dist) {
    FILE *outputf = fopen("output.txt", "w");
    if (!outputf) {
        perror("Failed to open output file");
        return -1;
    }

    if (!has_negative_cycle) {
        for (int i = 0; i < N; i++) {
            // if (dist[i] > INF)
            //     dist[i] = INF;
            fprintf(outputf, "%d\n", dist[i]);
        }
        fflush(outputf);
    } else {
        fprintf(outputf, "FOUND NEGATIVE CYCLE!\n");
    }
    
    fclose(outputf);
    return 0;
}

int read_file(const char* filename) {
    FILE* inputf = fopen(filename, "r");
    if (inputf == NULL) {
        abort_with_error_message("ERROR OCCURRED WHILE READING INPUT FILE");
    }
    
    if (fscanf(inputf, "%d", &N) != 1) {
        abort_with_error_message("ERROR READING MATRIX SIZE FROM FILE");
    }
    
    assert(N < (1024 * 1024 * 20));
    
    mat = (int *)malloc(N * N * sizeof(int));
    if (mat == NULL) {
        abort_with_error_message("MEMORY ALLOCATION FAILED");
    }

    int o;
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (fscanf(inputf, "%d", &mat[convert_dimension_2D_1D(i, j, N)]) != 1) {
                abort_with_error_message("ERROR READING MATRIX ELEMENT FROM FILE");
            }
        }
    }

    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         if (fscanf(inputf, "%d", &o) != 1) {
    //             abort_with_error_message("ERROR READING MATRIX ELEMENT FROM FILE");
    //         } else {
    //             fprintf(stdout, "%d\n", o);
    //         }
    //     }
    // }
    
    fclose(inputf);
    return 0;
}

// int read_file(string filename) {
//     std::ifstream inputf(filename, std::ifstream::in);
//     if (!inputf.good()) {
//         abort_with_error_message("ERROR OCCURRED WHILE READING INPUT FILE");
//     }
//     inputf >> N;
//     //input matrix should be smaller than 20MB * 20MB (400MB, we don't have too much memory for multi-processors)
//     assert(N < (1024 * 1024 * 20));
//     mat = (int *) malloc(N * N * sizeof(int));
//     for (int i = 0; i < N; i++)
//         for (int j = 0; j < N; j++) {
//             inputf >> mat[convert_dimension_2D_1D(i, j, N)];
//         }
//     return 0;
// }

void bellman_ford(int my_rank, int p, MPI_Comm comm, int n, int *mat, int *dist, bool *has_negative_cycle, int source) {
    int loc_n; // need a local copy for N
    int loc_start, loc_end;
    int *loc_mat; //local matrix
    int *loc_dist; //local distance

    //step 1: broadcast N
    if (my_rank == 0) {
        loc_n = n;
    }
    MPI_Bcast(&loc_n, 1, MPI_INT, 0, comm);

    //step 2: find local task range
    int ave = loc_n / p;
    loc_start = ave * my_rank;
    loc_end = ave * (my_rank + 1);
    if (my_rank == p - 1) {
        loc_end = loc_n;
    }

    //step 3: allocate local memory
    loc_mat = (int *) malloc(loc_n * loc_n * sizeof(int));
    loc_dist = (int *) malloc(loc_n * sizeof(int));

    //step 4: broadcast matrix mat
    if (my_rank == 0)
        memcpy(loc_mat, mat, sizeof(int) * loc_n * loc_n);
    MPI_Bcast(loc_mat, loc_n * loc_n, MPI_INT, 0, comm);

    //step 5: bellman-ford algorithm
    for (int i = 0; i < loc_n; i++) {
        loc_dist[i] = INF;
    }
    loc_dist[source] = 0;
    MPI_Barrier(comm);

    if (my_rank == 0) {
        printf("[before] loc dist: ");
        for (int i = 0; i < loc_n; i++) {
            printf("%d ", loc_dist[i]);
        }
        printf("\n");
    }

    bool loc_has_change;
    int loc_iter_num = 0;
    for (int iter = 0; iter < loc_n - 1; iter++) {
        loc_has_change = false;
        loc_iter_num++;
        for (int u = loc_start; u < loc_end; u++) {
            for (int v = 0; v < loc_n; v++) {
                int weight = loc_mat[convert_dimension_2D_1D(u, v, loc_n)];
                if (weight < INF) {
                    if (loc_dist[u] + weight < loc_dist[v]) {
                        loc_dist[v] = loc_dist[u] + weight;
                        loc_has_change = true;
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &loc_has_change, 1, MPI_CXX_BOOL, MPI_LOR, comm);
        if (!loc_has_change)
            break;
        MPI_Allreduce(MPI_IN_PLACE, loc_dist, loc_n, MPI_INT, MPI_MIN, comm);
        
        if (my_rank == 0) {
            printf("[%d] loc dist: ", iter);
            for (int i = 0; i < loc_n; i++) {
                printf("%d ", loc_dist[i]);
            }
            printf("\n");
        }
    }

    //do one more step
    if (loc_iter_num == loc_n - 1) {
        loc_has_change = false;
        for (int u = loc_start; u < loc_end; u++) {
            for (int v = 0; v < loc_n; v++) {
                int weight = loc_mat[convert_dimension_2D_1D(u, v, loc_n)];
                if (weight < INF) {
                    if (loc_dist[u] + weight < loc_dist[v]) {
                        loc_dist[v] = loc_dist[u] + weight;
                        loc_has_change = true;
                        break;
                    }
                }
            }
        }
        MPI_Allreduce(&loc_has_change, has_negative_cycle, 1, MPI_CXX_BOOL, MPI_LOR, comm);
    }

    //step 6: retrieve results back
    if(my_rank == 0)
        memcpy(dist, loc_dist, loc_n * sizeof(int));

    //step 7: remember to free memory
    free(loc_mat);
    free(loc_dist);

}

int main(int argc, char **argv) {

    int test = INT_MAX;

    if (argc <= 2) {
        abort_with_error_message("INPUT FILE WAS NOT FOUND!");
    }
    char* filename = argv[1];
    int source = atoi(argv[2]);

    int *dist;
    bool has_negative_cycle = false;

    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int p;//number of processors
    int my_rank;//my global rank
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    if (my_rank == 0) {
        printf("Int max: %d \n", test);
    }

    //only rank 0 process do the I/O
    if (my_rank == 0) {
        assert(read_file(filename) == 0);
        dist = (int *) malloc(sizeof(int) * N);
    }

    //time counter
    double t1, t2;
    MPI_Barrier(comm);
    t1 = MPI_Wtime();

    //bellman-ford algorithm
    bellman_ford(my_rank, p, comm, N, mat, dist, &has_negative_cycle, source);
    MPI_Barrier(comm);

    //end timer
    t2 = MPI_Wtime();

    if (my_rank == 0) {
        printf("t1: %d \n", t1);
        printf("t2: %d \n", t2);
    }

    if (my_rank == 0) {
        fprintf(stdout, "Time (s): %d\n", t2 - t1);
        print_result(has_negative_cycle, dist);
        free(dist);
        free(mat);
    }

    MPI_Finalize();
    return 0;
}