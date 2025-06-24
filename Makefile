# Compilers
CC = gcc
MPICC = mpicc

# Compiler flags
CFLAGS = -g -Wall -Wextra -fopenmp -I./include
MPI_FLAGS = -g -Wall -Wextra -fopenmp -I./include
NO_OMP_FLAGS = -Wall -Wextra -I./include

# CFLAGS = -fopenmp -I./include
# MPI_FLAGS = -fopenmp -I./include
# NO_OMP_FLAGS = -I./include
#
# Source directory
SRC_DIR = src

# Source files
SRC_SERIAL       = $(SRC_DIR)/stencil_template_serial.c
SRC_SERIAL_FR    = $(SRC_DIR)/stencil_template_serial_fr.c
SRC_PARALLEL     = $(SRC_DIR)/stencil_template_parallel.c

# Object files
OBJ_SERIAL       = $(SRC_SERIAL:.c=.o)
OBJ_SERIAL_FR    = $(SRC_SERIAL_FR:.c=.o)
OBJ_PARALLEL     = $(SRC_PARALLEL:.c=.o)

# Executables
EXEC_SERIAL      = stencil_serial
EXEC_SERIAL_FR   = stencil_serial_fr
EXEC_PARALLEL    = stencil_parallel

# Default target: build everything
all: $(EXEC_SERIAL) $(EXEC_SERIAL_FR) $(EXEC_PARALLEL)

# Build only serial versions (OpenMP + fallback)
openmp: $(EXEC_SERIAL) $(EXEC_SERIAL_FR)

# Serial with OpenMP
$(EXEC_SERIAL): $(OBJ_SERIAL)
	$(CC) $(CFLAGS) -o $@ $^

# Fallback serial without OpenMP
$(EXEC_SERIAL_FR): $(OBJ_SERIAL_FR)
	$(CC) $(NO_OMP_FLAGS) -o $@ $^

# Parallel (with MPI + OpenMP)
$(EXEC_PARALLEL):
	$(MPICC) $(MPI_FLAGS) -o $@ $(SRC_PARALLEL)

# Individual object compilation
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(SRC_DIR)/*.o $(EXEC_SERIAL) $(EXEC_SERIAL_FR) $(EXEC_PARALLEL)

.PHONY: all clean openmp

# # Compilers
# CC = gcc
# MPICC = mpicc

# # Compiler flags
# CFLAGS = -Wall -Wextra -fopenmp -I./include
# MPI_FLAGS = -Wall -Wextra -I./include


# # Source and object files
# SRC_DIR = src
# SRCS = $(wildcard $(SRC_DIR)/*.c)
# OBJS = $(SRCS:.c=.o)

# # Executables
# EXEC_SERIAL = stencil_serial
# EXEC_PARALLEL = stencil_parallel

# # Default target
# all: $(EXEC_SERIAL) $(EXEC_PARALLEL)

# # Serial version (compiled with gcc)
# $(EXEC_SERIAL): $(SRC_DIR)/stencil_template_serial.o
# 	$(CC) $(CFLAGS) -o $@ $^

# # Parallel version (compiled with mpicc)
# $(EXEC_PARALLEL): $(SRC_DIR)/stencil_template_parallel.o
# 	$(MPICC) $(MPI_FLAGS) -o $@ $^ -fopenmp

# # Rule for serial object file
# $(SRC_DIR)/stencil_template_serial.o: $(SRC_DIR)/stencil_template_serial.c
# 	$(CC) $(CFLAGS) -c $< -o $@

# # Rule for parallel object file (MPI)
# $(SRC_DIR)/stencil_template_parallel.o: $(SRC_DIR)/stencil_template_parallel.c
# 	$(MPICC) $(MPI_FLAGS) -c $< -o $@

# # Clean target
# clean:
# 	rm -f $(OBJS) $(EXEC_SERIAL) $(EXEC_PARALLEL)

# .PHONY: all clean
