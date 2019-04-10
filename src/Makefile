.PHONY: clean coma ofp K
DEBUG=-DNDEBUG
#########################################################
CC=mpicc
CFLAGS=-std=gnu99 -O3 -Wall -mcmodel=medium $(DEBUG)
OPENMP=-fopenmp
LFLAGS=-lm
#########################################################
ifeq ($(ENV), coma)
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
CFLAGS=-D_ENV_K -O3 -Kopenmp,dalign,lib,mfunc,ns,omitfp,prefetch_conditional,rdconv -x- -noansi -Xg $(DEBUG)
LFLAGS=
endif
#########################################################
SRCS := main.c sa.c timer.c bfs.c common.c
OBJS := $(SRCS:%.c=%.o)

gg: $(OBJS)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $^

%.o: %.c common.h
	$(CC) $(CFLAGS) -c $<

common_omp.o: common.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

sa_omp.o: sa.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

bfs_omp.o: bfs.c common.h
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $<

gg_thread: common_omp.o sa_omp.o bfs_omp.o main.o timer.o
	$(CC) $(CFLAGS) $(OPENMP) $(LFLAGS) -o $@ $^

all: gg gg_thread
#########################################################
# make K DATA="..\/data\/n72d4.random.edges" N=10000 G=2 K=1 T=100 C=1
K:
	make clean; make ENV=K
	@for s in `seq 1 $(K)`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/K/submit.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" | sed "s/@T@/$(T)/" | sed "s/@C@/$(C)/" | sed "s/@S@/$$s/" > submit.sh; \
	ksub --rsc-list rscgrp=micro ../scripts/K/submit.sh; \
	done

# make fruits DATA="..\/data\/n16d5.random.edges" N=1000000 G=2 K=1
fruits:
	cp ../scripts/ana-5-0.5.sh ana.sh
	make clean; make ENV=fruits
	@for t in 5 4.770477 4.551491 4.342557 4.143214 3.953022 3.77156 3.598428 3.433244 3.275643 \
	3.125276 2.981812 2.844933 2.714338 2.589737 2.470857 2.357433 2.249216 2.145967 2.047458 \
	1.953470 1.863797 1.778240 1.696611 1.618729 1.544422 1.473526 1.405884 1.341348 1.279774 \
	1.221027 1.164976 1.111498 1.060475 1.011795 0.965349 0.921035 0.878755 0.838416 0.799929 \
	0.763209 0.728174 0.694748 0.662856 0.632428 0.603396 0.575698 0.549271 0.524057 0.5 SA1 SA2; do \
	for s in `seq 1 $(K)`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/fruits/submit.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" | sed "s/@T@/$$t/" | sed "s/@S@/$$s/" > submit.sh; \
	sbatch submit.sh; \
	done; \
	done

# make fruits DATA="..\/data\/n16d5.random.edges" N=1000000 G=2 W=332.9740153 c=0.667616
fruits-cooling:
	cp ../scripts/fruits/ana_cooling.sh ana.sh
	make clean; make ENV=fruits
	for C in 1 10 100 1000 10000; do \
	for s in `seq 1 10`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/fruits/submit_cooling.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" > /tmp/a; \
	sed "s/@W@/$$W/" /tmp/a | sed "s/@c@/$$c/" | sed "s/@C@/$$C/" | sed "s/@S@/$$s/" > submit.sh; \
	sbatch submit.sh; rm /tmp/a; \
	done; \
	done

#make coma DATA="..\/data\/n72d4.random.edges" N=1000000 G=1
#make coma DATA="..\/data\/n36d4.random.edges" N=1000000 G=2
#make coma DATA="..\/data\/n24d4.random.edges" N=1000000 G=3
#make coma DATA="..\/data\/n18d4.random.edges" N=1000000 G=4
#make coma DATA="..\/data\/n12d4.random.edges" N=1000000 G=6
#make coma DATA="..\/data\/n9d4.random.edges" N=1000000 G=8
#make coma DATA="..\/data\/n8d4.random.edges" N=1000000 G=9
#make coma DATA="..\/data\/n6d4.random.edges" N=1000000 G=12
coma: 
	cp ../scripts/coma/ana-5-0.5.sh ana.sh
	make clean; make ENV=coma
	@for i in `seq 0 9`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/coma/submit-$$i.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" > submit-$$i.sh; \
	sbatch submit-$$i.sh; \
	done

#make ofp DATA="..\/data\/n6d4.random.edges" N=1000000 G=12
ofp:
	cp ../scripts/ana.sh .
	make clean; make ENV=ofp
	@for i in `seq 0 16`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/ofp/submit-$$i.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" > submit-$$i.sh; \
	pjsub submit-$$i.sh; \
	done

clean:
	rm -rf *.o *~

