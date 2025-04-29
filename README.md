# bellman-ford-mpi

Parallel Bellman-Ford

```bash
mpicc -o belman-ford-better belman-ford-better.c
mpiexec -f nodes -n 4 ./belman-ford-better input_simple.txt 0 4
```
