# Compilers
CC = gcc
MPICC = mpicc

# Compiler flags
CFLAGS = -Wall -Wextra -fopenmp -I./include
MPI_FLAGS = -Wall -Wextra -I./include


# Source and object files
SRC_DIR = src
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:.c=.o)

# Executables
EXEC_SERIAL = stencil_serial
EXEC_PARALLEL = stencil_parallel

# Default target
all: $(EXEC_SERIAL) $(EXEC_PARALLEL)

# Serial version (compiled with gcc)
$(EXEC_SERIAL): $(SRC_DIR)/stencil_template_serial.o
	$(CC) $(CFLAGS) -o $@ $^

# Parallel version (compiled with mpicc)
$(EXEC_PARALLEL): $(SRC_DIR)/stencil_template_parallel.o
	$(MPICC) $(MPI_FLAGS) -o $@ $^ -fopenmp

# Rule for serial object file
$(SRC_DIR)/stencil_template_serial.o: $(SRC_DIR)/stencil_template_serial.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule for parallel object file (MPI)
$(SRC_DIR)/stencil_template_parallel.o: $(SRC_DIR)/stencil_template_parallel.c
	$(MPICC) $(MPI_FLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJS) $(EXEC_SERIAL) $(EXEC_PARALLEL)

.PHONY: all clean
