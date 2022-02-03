CXX=g++-11
CC=gcc-11
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

IDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

OBJS1 = p_grid.o swirl.o particle.o common.o swirl_io.o wall.o swirl_param.o dat_store.o AYaux.o AYlinalg.o AYmat.o AYvec.o AYrng.o AYtens.o AYdecomps.o AYdata.o RYdat2AYdat.o

%.o:%.c
	$(CC) $(CFLAGS) $(IDIR) -c $<

%.o:%.cc
	$(CXX) $(CFLAGS) $(IDIR) -c $<

all: swirl_binary_test

swirl_binary_test: swirl_binary_test.cc $(OBJS1)
	$(CXX) $(CFLAGS) $(IDIR) -o $@ $^ $(LINK)

clean:
	rm -f swirl_binary_test
	rm -f *.o
clean_dat:
	rm -f *.aysml
	rm -f *.rysml
	rm -f *.aydat
clean_datdir:
	rm -f ./dat_dir/*.aysml
	rm -f ./dat_dir/*.rysml
	rm -f ./dat_dir/*.aydat

clean_all: clean clean_dat
