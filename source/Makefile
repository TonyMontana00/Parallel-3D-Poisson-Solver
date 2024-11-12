include Makefile.in

.PHONY: default clean

default: cgp3d.x cgp3d_parallel.x

# Compiler and flags
CC = gcc
CFLAGS = -O3 -DCLOCK=CLOCK_MONOTONIC
LIBS = -lm

# Object files
OBJS = cgp3d.o pcg.o params.o timing.o
PAR_OBJS = cgp3d_par.o pcg_par.o params_par.o timing_par.o

# Sequential version (no parallelization)
cgp3d.x: $(OBJS)
	$(CC) -no-pie -o $@ $(OBJS) $(LIBS)

# Parallel version with -DPARALLEL
cgp3d_parallel.x: $(PAR_OBJS)
	$(CC) -fopenmp -no-pie -o $@ $(PAR_OBJS) $(LIBS)

# Compilation rule for sequential object files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Compilation rule for parallel object files
%_par.o: %.c
	$(CC) $(CFLAGS) -fopenmp -DPARALLEL $(INCLUDE) -c $< -o $@

# Clean up the binaries and object files
clean:
	rm -f *.o *~ cgp3d.x cgp3d_parallel.x
