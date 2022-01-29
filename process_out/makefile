CXX=g++-11
CC=gcc-11
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

CDEPS = nrutil.h auxiliary_functions.h knuth_lcg.h

CPPDEPS = AYlinalg.hh

OBJS1=main1.o nrutil.o auxiliary_functions.o knuth_lcg.o AYlinalg.o AYmat.o AYvec.o AYrng.o AYtens.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: process_circ6

process_circ6: $(OBJS1)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

clean:
	rm -f process_circ6
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_datdir:
	rm -f ./dat_dir/*.dat
	rm -f ./dat_dir/*.aysml
	rm -f ./dat_dir/*.aydat

clean_all: clean clean_dat
