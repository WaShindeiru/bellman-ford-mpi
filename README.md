# bellman-ford-mpi

Parallel Bellman-Ford

### How to run:

#### Make:

```bash
source /opt/nfs/config/source_mpich420.sh
source /opt/nfs/config/source_cuda121.sh
make build
make run MPI_FLAGS="-f nodes -n 80" RUN_FLAGS="./input/input_300.txt 3 241"
make clean
```

#### CMake:

```bash
source /opt/nfs/config/source_mpich420.sh
source /opt/nfs/config/source_cuda121.sh
export MPICH_TARGET_DIR=/opt/nfs/mpich-4.2.0
mkdir build
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
