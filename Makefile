CC = mpicc
CFLAGS = -std=c17 -Iinclude -I$(MPI_HOME)/include
LDFLAGS = -L$(MPI_HOME)/lib -lmpi

SRC = src/bellman-ford.c src/graph.c
OBJ = $(SRC:.c=.o)
TARGET = bellman_ford

MPI_HOME ?= /opt/nfs/mpich-4.2.0

MPI_FLAGS ?= -f nodes -n 80
RUN_FLAGS ?= ./input/input_300.txt 3 241

.PHONY: all build run clean

all: build

build: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

run: build
	mpiexec $(MPI_FLAGS) ./$(TARGET) $(RUN_FLAGS) ./output.txt

clean:
	rm -f $(OBJ) $(TARGET) ./output.txt

help:
	@echo "Usage:"
	@echo "  make build    - Compile the project"
	@echo "  make run      - Build and execute with mpiexec"
	@echo "  make clean    - Remove generated files"
	@echo "  make all      - Default target (same as 'build')"
	@echo ""
	@echo "Current configuration:"
	@echo "  MPI_HOME = $(MPI_HOME)"
	@echo "  Run flags = $(RUN_FLAGS)"
	@echo "  MPI args = -n $(MPI_FLAGS)"