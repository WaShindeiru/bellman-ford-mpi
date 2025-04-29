#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>

#include "mpi.h"

#define INF INT_MAX / 2

typedef struct AdjListNode AdjListNode;

struct AdjListNode {
    int dest;
    int length;
    AdjListNode* next;
};

typedef struct {
    AdjListNode *head;
    int size;
} AdjList;

typedef struct {
    int V;
    AdjList* array;
} Graph;

void freeGraph(Graph* graph) {
    if (!graph) return;
    
    if (graph->array) {
        for (int v = 0; v < graph->V; v++) {
            AdjListNode* curr = graph->array[v].head;
            while (curr) {
                AdjListNode* next = curr->next;
                free(curr);
                curr = next;
            }
        }
        free(graph->array);
    }
    free(graph);
}

AdjListNode* newAdjListNode(int dest, int length) {
    AdjListNode* newNode = (AdjListNode*) malloc(sizeof(AdjListNode));
    newNode->dest = dest;
    newNode->length = length;
    newNode->next = NULL;

    return newNode;
}

Graph* createGraph(int V) {
    Graph* graph = (Graph*) malloc(sizeof(Graph));
    graph->V = V;
    graph->array = (AdjList*) malloc(V * sizeof(AdjList));

    int i;
    for (i = 0; i < V; ++i) {
        graph->array[i].head = NULL;
        graph->array[i].size = 0;
    }

    return graph;
}

void addEdge(Graph* graph, int src, int dest, int length) {
    AdjListNode* newNode = newAdjListNode(dest, length);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
    graph->array[src].size++;
}

Graph* copyGraph(Graph* graph) {
    int V = graph->V;

    Graph* result = createGraph(V);
    for (int v = 0; v < V; v++) {
        int size = graph->array[v].size;
        AdjListNode* pCrawl = graph->array[v].head;

        while (pCrawl) {
            addEdge(result, v, pCrawl->dest, pCrawl->length);
            pCrawl = pCrawl->next;
        }
    }

    return result;
}

void printGraph(Graph* graph) {
    int v;
    for (v = 0; v < graph->V; v++) {
        int size = graph->array[v].size;
        AdjListNode* pCrawl = graph->array[v].head;
        printf("\n Adjacency list of vertex %d, size: %d \n head ", v, size);
        while (pCrawl) {
            printf("-> dest: %d, length: %d \n", pCrawl->dest, pCrawl->length);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}

int N; //number of vertices
Graph* graph; // adjacency list

void abort_with_error_message(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    abort();
}

int print_result(bool has_negative_cycle, int *dist, int* path, int dest) {
    FILE *outputf = fopen("output_newest.txt", "w");
    if (!outputf) {
        perror("Failed to open output file");
        return -1;
    }

    // int* temp = (int *) malloc(N * sizeof(int));

    // int i;
    // int index;
    // for (i=N-1, index=dest; index != -1; index=path[index], i--) {
    //     temp[i] = index;
    // }

    // fprintf(outputf, "path: ");
    // for (i=0; i<N; i++) {
    //     fprintf(outputf, "%d, ", temp[i]);
    // }
    // fprintf(outputf, "\n");

    fprintf(outputf, "Node/Parent:");
    int i;
    for (i=0; i<N; i++) {
        fprintf(outputf, " (%d, %d)", i, path[i]);
    }
    fprintf(outputf, "\n");

    if (!has_negative_cycle) {
        for (int i = 0; i < N; i++) {
            if (dist[i] > INF)
                dist[i] = INF;
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
    char temp[100];
    FILE* inputf = fopen(filename, "r");
    if (inputf == NULL) {
        abort_with_error_message("ERROR OCCURRED WHILE READING INPUT FILE");
    }
    
    if (fscanf(inputf, "%d", &N) != 1) {
        abort_with_error_message("ERROR READING MATRIX SIZE FROM FILE");
    }
    
    assert(N < (1024 * 1024 * 20));
    
    graph = createGraph(N);
    if (graph == NULL) {
        abort_with_error_message("MEMORY ALLOCATION FAILED");
    }

    int o;

    int dest;
    int length;
    
    for (int i = 0; i < N; i++) {
        if (fscanf(inputf, "%d", &o) != 1) {
            abort_with_error_message("ERROR READING MATRIX SIZE FROM FILE");
        }

        for (int j = 0; j < o; j++) {
            if (fscanf(inputf, "%d,%d", &dest, &length) != 2) {
                sprintf(temp, "ERROR READING MATRIX ELEMENT FROM FILE, i: %d, j: %d", i, j);
                abort_with_error_message(temp);
            }

            addEdge(graph, i, dest, length);
        }
    }
    
    fclose(inputf);
    return 0;
}

typedef struct {
    int dist;
    int rank;
} helper;

void bellman_ford(int my_rank, int p, MPI_Comm comm, int n, Graph* graph, int *dist, int* path, bool *has_negative_cycle, int source) {
    int loc_n; // need a local copy for N
    int loc_start, loc_end;
    Graph* loc_graph;  //local graph
    int* loc_dist; //local distance
    int* loc_path;
    helper* dist_helper;
    int* all_path_values;

    //step 1: broadcast N
    if (my_rank == 0) {
        loc_n = n;
    }
    MPI_Bcast(&loc_n, 1, MPI_INT, 0, comm);

    //step 3: allocate local memory
    loc_dist = (int *) malloc(loc_n * sizeof(int));
    loc_path = (int *) malloc(loc_n * sizeof(int));
    dist_helper = (helper *) malloc(loc_n * sizeof(helper));
    all_path_values = malloc(p * loc_n * sizeof(int));
    
    if (my_rank == 0) {
        loc_graph = graph;
    } else {
        loc_graph = createGraph(loc_n);
    }

    int local_size;
    int dest;
    int length;
    AdjListNode* pCrawl;

    //broadcast adjacency list
    for (int i = 0; i < loc_n; i++) {
        if (my_rank == 0) {
            local_size = graph->array[i].size;
        }
        MPI_Bcast(&local_size, 1, MPI_INT, 0, comm);

        if (my_rank == 0) {
            pCrawl = graph->array[i].head;
        }

        for (int j = 0; j < local_size; j++) {
            if (my_rank == 0) {
                dest = pCrawl->dest;
                length = pCrawl->length;
                pCrawl = pCrawl->next;
            }

            MPI_Bcast(&dest, 1, MPI_INT, 0, comm);
            MPI_Bcast(&length, 1, MPI_INT, 0, comm);

            if (my_rank != 0) {
                addEdge(loc_graph, i, dest, length);
            }
        }
    }

    if (my_rank == 0) {
        loc_graph = copyGraph(graph);
    }

    //step 2: find local task range
    int ave = loc_n / p;
    loc_start = ave * my_rank;
    loc_end = ave * (my_rank + 1);
    if (my_rank == p - 1) {
        loc_end = loc_n;
    }

    //step 5: bellman-ford algorithm
    for (int i = 0; i < loc_n; i++) {
        loc_dist[i] = INF;
    }

    for (int i = 0; i < loc_n; i++) {
        loc_path[i] = -1;
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
        printf("Enter iteration: %d\n", iter);
        loc_has_change = false;
        loc_iter_num++;
        for (int u = loc_start; u < loc_end; u++) {

            local_size = loc_graph->array[u].size;
            pCrawl = loc_graph->array[u].head;

            for (int v = 0; v < local_size; v++) {
                dest = pCrawl->dest;
                length = pCrawl->length;
                pCrawl = pCrawl->next;

                if (length < INF) {
                    if (loc_dist[u] + length < loc_dist[dest]) {
                        loc_dist[dest] = loc_dist[u] + length;
                        loc_path[dest] = u;
                        loc_has_change = true;
                    }
                }
            }

            // for (int v = 0; v < loc_n; v++) {
            //     int weight = loc_mat[convert_dimension_2D_1D(u, v, loc_n)];
            //     if (weight < INF) {
            //         if (loc_dist[u] + weight < loc_dist[v]) {
            //             loc_dist[v] = loc_dist[u] + weight;
            //             loc_has_change = true;
            //         }
            //     }
            // }
        }

        MPI_Allreduce(MPI_IN_PLACE, &loc_has_change, 1, MPI_CXX_BOOL, MPI_LOR, comm);
        if (!loc_has_change) {
            printf("no change, iteration: %d\n", iter);
            break;
        }

        for (int i=0; i<loc_n; i++) {
            dist_helper[i].dist = loc_dist[i];
            dist_helper[i].rank = my_rank;
        }
        
        MPI_Allreduce(MPI_IN_PLACE, dist_helper, loc_n, MPI_2INT, MPI_MINLOC, comm);

        for (int i=0; i<loc_n; i++) {
            loc_dist[i] = dist_helper[i].dist;
        }

        MPI_Allgather(loc_path, loc_n, MPI_INT, all_path_values, loc_n, MPI_INT, comm);

        for (int i=0; i<loc_n; i++) {
            loc_path[i] = all_path_values[dist_helper[i].rank * loc_n + i];
        }
        
        // if (my_rank == 0) {
        //     printf("[%d] loc dist: ", iter);
        //     for (int i = 0; i < loc_n; i++) {
        //         printf("%d ", loc_dist[i]);
        //     }
        //     printf("\n");
        // }
        printf("Exit iteration: %d\n", iter);
    }

    //todo: refactor
    //do one more step
    // no
    // if (loc_iter_num == loc_n - 1) {
    //     loc_has_change = false;
    //     for (int u = loc_start; u < loc_end; u++) {
    //         for (int v = 0; v < loc_n; v++) {
    //             int weight = loc_mat[convert_dimension_2D_1D(u, v, loc_n)];
    //             if (weight < INF) {
    //                 if (loc_dist[u] + weight < loc_dist[v]) {
    //                     loc_dist[v] = loc_dist[u] + weight;
    //                     loc_has_change = true;
    //                     break;
    //                 }
    //             }
    //         }
    //     }
    //     MPI_Allreduce(&loc_has_change, has_negative_cycle, 1, MPI_CXX_BOOL, MPI_LOR, comm);
    // }

    //step 6: retrieve results back
    if(my_rank == 0) {
        memcpy(dist, loc_dist, loc_n * sizeof(int));
        memcpy(path, loc_path, loc_n * sizeof(int));
    }

    //step 7: remember to free memory
    // todo: free memory
    // free(loc_mat);
    
    free(loc_dist);
    freeGraph(loc_graph);
    free(dist_helper);
    free(all_path_values);
    loc_dist = NULL;
}

int main(int argc, char **argv) {
    if (argc <= 2) {
        abort_with_error_message("INPUT FILE WAS NOT FOUND!");
    }

    char* filename = argv[1];
    int source = atoi(argv[2]);
    int dest = atoi(argv[3]);

    int *dist;
    int* path;
    bool has_negative_cycle = false;

    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int p;//number of processors
    int my_rank;//my global rank
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    //only rank 0 process do the I/O
    if (my_rank == 0) {
        assert(read_file(filename) == 0);
        dist = (int *) malloc(sizeof(int) * N);
        path = (int* ) malloc(sizeof(int) * N);
    }

    double t1, t2;
    MPI_Barrier(comm);
    t1 = MPI_Wtime();

    //bellman-ford algorithm
    bellman_ford(my_rank, p, comm, N, graph, dist, path, &has_negative_cycle, source);
    MPI_Barrier(comm);

    //end timer
    t2 = MPI_Wtime();

    if (my_rank == 0) {
        printf("t1: %d \n", t1);
        printf("t2: %d \n", t2);
    }

    if (my_rank == 0) {
        fprintf(stdout, "Time (s): %d\n", t2 - t1);
        print_result(has_negative_cycle, dist, path, source);
        free(dist);
        freeGraph(graph);
        free(path);
    }

    if (my_rank != 0) {
        fprintf(stdout, "debug\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}