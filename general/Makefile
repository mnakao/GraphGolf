.PHONY: clean cygnus ofp
DEBUG=-DNDEBUG
#########################################################
CC=mpicc
CFLAGS=-std=gnu99 -O3 -Wall -mcmodel=medium -Wno-unknown-pragmas -march=native $(DEBUG)
OPENMP=-fopenmp
LFLAGS=-lm
NVCC=nvcc
NVCC_FLAG=-D__C2CUDA__
NVCC_LDFLAG=-lcuda
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
#################################################
gg_cuda: main_cuda.o common.o sa.o timer.o matrix_cuda.o
	$(NVCC) -O3 $^ -o $@ $(NVCC_LDFLAG) -ccbin $(CC)

main_cuda.o: main.c common.h
	$(CC) $(CFLAGS) $< -c -o $@ $(NVCC_FLAG)

matrix_cuda.o: matrix.cu common.h
	$(NVCC) -O3 $< -c -o $@ $(NVCC_FLAG) -ccbin $(CC)
#########################################################
cuda: gg_cuda
openmp: gg_openmp
all: gg openmp
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
          ./gg -f data/n114d29.random.edges -s $$s -n 100 -v 16 -g 15;\
	done >> result.txt
	diff result.txt data/result.txt | grep -v sec | grep -v 2019 | grep -Ev "[0-9]c" | sort | uniq
	rm -f result.txt
