# Compiler and compilation flags
cc=gcc-11
cxx=g++-11 -fopenmp
cflags=-Wall -ansi -pedantic -march=native -O3

# MPI compiler
mpicxx=mpicxx -Wno-long-long

# Flags for linking to PNG library
png_iflags=-DHAS_PNG -I/opt/local/include
png_lflags=-L/opt/local/lib -lpng

# Flags for FFTW library
fftw_iflags=-I/opt/local/include
fftw_lflags=-L/opt/local/lib -lfftw3

# Flags for Eigen library
eigen_iflags=-isystem/opt/local/include/eigen3

# LAPACK flags for linear algebra
lp_lflags=-framework Accelerate
