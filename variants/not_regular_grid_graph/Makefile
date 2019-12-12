.PHONY: clean cygnus ofp
DEBUG=-DNDEBUG
#DEBUG=-D_DEBUG_MSG
#########################################################
CC=mpicc
CFLAGS=-std=gnu99 -O3 -Wall -Wno-unknown-pragmas -march=native $(DEBUG)
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
all: gg gg_openmp
#########################################################
clean:
	rm -rf *.o *~

clean_all: clean
	rm -rf gg gg_openmp

test: gg
	rm -f result.txt
	for s in $$(seq 1 10); do\
	  ./gg -f data/w5h5d4.random.edges   -s $$s -n 1000 -r 2 -W 2;\
          ./gg -f data/w5h10d4.random.edges  -s $$s -n 1000 -r 2 -W 2;\
          ./gg -f data/w10h10d4.random.edges -s $$s -n 1000 -r 2 -W 2;\
          ./gg -f data/w10h10d6.random.edges -s $$s -n 1000 -r 3 -W 2;\
          ./gg -f data/w10h10d8.random.edges -s $$s -n 1000 -r 4 -W 2;\
          ./gg -f data/w20h20d4.random.edges -s $$s -n 1000 -r 2 -W 2;\
	  ./gg -f data/w5h5d4.random.edges   -s $$s -n 1000 -r 2 -W 2 -g 2;\
	  ./gg -f data/w5h5d4.random.edges   -s $$s -n 1000 -r 2 -W 2 -g 4;\
	done >> result.txt
	diff result.txt data/result.txt | grep -v sec | grep -v 2019 | grep -Ev "[0-9]c" | sort | uniq
	rm -f result.txt
