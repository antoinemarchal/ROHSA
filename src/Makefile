#fortran compiler
FC=gfortran
# FC=mpif90

# GFLAGS=-Wall -Wextra -pedantic -std=f2008 -g -O2 -fno-second-underscore -Wno-compare-reals -fcheck=all -fimplicit-none -fbacktrace -gdwarf-2 -fall-intrinsics -Wno-unused-function -fcheck=bounds -I/usr/local/include
GFLAGS= -O -Wall -fbounds-check -g -Wno-uninitialized -I/usr/local/include

CFLAGS=-openmp -llapack -lcfitsio -lm -g -O2

OUT=ROHSA

all: blas.o lbfgsb.o linpack.o timer.o mod_constants.o mod_read.o mod_array.o mod_optimize.o mod_functions.o mod_start.o  main.o
	$(FC) $^ $(CFLAGS) -o $(OUT)

main.o: main.f90
	$(FC) $(GFLAGS) -c $^

# L-BFGS-B package
blas.o: blas.f
	$(FC) $(GFLAGS) -c $^
lbfgsb.o: lbfgsb.f
	$(FC) $(GFLAGS) -c $^
linpack.o: linpack.f
	$(FC) $(GFLAGS) -c $^
timer.o: timer.f
	$(FC) $(GFLAGS) -c $^

# My files
mod_constants.o: mod_constants.f90
	$(FC) $(GFLAGS) -c $^

mod_read.o: mod_read.f90
	$(FC) $(GFLAGS) -c $^

mod_array.o: mod_array.f90
	$(FC) $(GFLAGS) -c $^

mod_optimize.o: mod_optimize.f90
	$(FC) $(GFLAGS) -c $^

mod_functions.o: mod_functions.f90
	$(FC) $(GFLAGS) -c $^

mod_start.o: mod_start.f90
	$(FC) $(GFLAGS) -c $^

clean:
	\rm -rf *.o $(OUT) *.mod *~

run: all
	./$(OUT)

rerun: clean all run
