.PHONY: clean cygnus ofp K
DEBUG=-DNDEBUG
#########################################################
CC=mpicc
CFLAGS=-std=gnu99 -O3 -Wall -mcmodel=medium -Wno-unknown-pragmas -march=native $(DEBUG)
OPENMP=-fopenmp
LFLAGS=-lm
#########################################################
ifeq ($(ENV), cygnus)
CC=mpiicc
CFLAGS=-O2 -std=gnu99 $(DEBUG)
OPENMP=-qopenmp
LFLAGS=
else ifeq ($(ENV), ofp)
CC=mpiicc
CFLAGS=-axMIC-AVX512 -O3 -std=gnu99 $(DEBUG)
OPENMP=-qopenmp
LFLAGS=
else ifeq ($(ENV), K)
CC=mpifccpx
# With "-Kfast", a value of low_ASPL may be changed. 
CFLAGS=-D_KCOMPUTER -O3 -Kdalign,lib,mfunc,ns,omitfp,prefetch_conditional,rdconv -x- -noansi -Xg $(DEBUG)
OPENMP=-Kopenmp
LFLAGS=
endif
#########################################################
SRCS := main.c sa.c timer.c evaluation.c common.c
OBJS := $(SRCS:%.c=%.o)

gg: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.c common.h
	$(CC) $(CFLAGS) -c $<

common_omp.o: common.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

sa_omp.o: sa.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

evaluation_omp.o: evaluation.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

main_omp.o: main.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

timer_omp.o: timer.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

gg_openmp: common_omp.o sa_omp.o evaluation_omp.o main_omp.o timer_omp.o
	$(CC) $(CFLAGS) $(OPENMP) -o $@ $^ $(LFLAGS)

openmp: gg_openmp
all: gg openmp
#########################################################
clean:
	rm -rf *.o *~

clean_all: clean
	rm -f gg gg_openmp
########################################################
test: gg
	rm -f result.txt
	for s in $$(seq 1 10); do\
	  ./gg -f data/n512d4.random.edges  -s $$s -n 1000 ;\
          ./gg -f data/n256d4.random.edges  -s $$s -n 1000 -g 2;\
          ./gg -f data/n128d4.random.edges  -s $$s -n 1000 -g 4;\
          ./gg -f data/n64d4.random.edges   -s $$s -n 1000 -g 8;\
          ./gg -f data/n32d4.random.edges   -s $$s -n 1000 -g 16;\
          ./gg -f data/n16d4.random.edges   -s $$s -n 1000 -g 32;\
          ./gg -f data/n8d4.random.edges    -s $$s -n 1000 -g 64;\
          ./gg -f data/n114d29.random.edges -s $$s -n 1000 -v 16 -g 15;\
	done >> result.txt
	diff result.txt data/result.org | grep -v sec | sort | uniq
