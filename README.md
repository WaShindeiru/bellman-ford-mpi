# bellman-ford-mpi

Parallel Bellman-Ford

```bash
export MPICH_TARGET_DIR=/opt/nfs/mpich-4.2.0
cd build
cmake ..
make
cd ..
mpiexec -f nodes -n 80 ./build/bellman_ford ./input/input_300.txt 3 241
```

<!-- ```bash
mpicc -o belman-ford-better belman-ford-better.c
mpiexec -f nodes -n 4 ./belman-ford-better input_simple.txt 0 4
``` -->
