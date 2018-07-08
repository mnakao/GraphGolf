.PHONY: clean coma
CC=mpicc
#DEBUG=-DNDEBUG

CFLAGS=-std=gnu99 -O3 -Wall -mcmodel=medium -fopenmp $(DEBUG)
LFLAGS=-lm
ifeq ($(ENV), coma)
CC=mpiicc
CFLAGS=-O3 -qopenmp -std=gnu99 $(DEBUG)
LFLAGS=
else ifeq ($(ENV), ofp)
CC=mpiicc
CFLAGS=-axMIC-AVX512 -O3 -std=gnu99 -qopenmp
LFLAGS=
endif
SRCS := main.c sa.c timer.c top_down_bfs.c common.c
OBJS := $(SRCS:%.c=%.o)

all: a.out

a.out: $(OBJS)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $^

%.o: %.c common.h
	$(CC) $(CFLAGS) -c $<

run: a.out
	nohup bash run.sh &

# make fruits DATA="..\/data\/n16d5.random.edges" N=1000000 G=2
fruits:
	cp ../scripts/ana-5-0.5.sh ana.sh
	make clean; make ENV=fruits
	@for t in 5 4.770477 4.551491 4.342557 4.143214 3.953022 3.77156 3.598428 3.433244 3.275643 \
	3.125276 2.981812 2.844933 2.714338 2.589737 2.470857 2.357433 2.249216 2.145967 2.047458 \
	1.953470 1.863797 1.778240 1.696611 1.618729 1.544422 1.473526 1.405884 1.341348 1.279774 \
	1.221027 1.164976 1.111498 1.060475 1.011795 0.965349 0.921035 0.878755 0.838416 0.799929 \
	0.763209 0.728174 0.694748 0.662856 0.632428 0.603396 0.575698 0.549271 0.524057 0.5 SA1 SA2; do \
	for s in `seq 1 10`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/fruits/submit.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" | sed "s/@T@/$$t/" | sed "s/@S@/$$s/" > submit.sh; \
	sbatch submit.sh; \
	done; \
	done

# make coma DATA="..\/data\/n16d5.random.edges" N=1000000 G=2
coma: 
	cp ../scripts/ana-5-0.5.sh ana.sh
	make clean; make ENV=coma
	@for i in `seq 0 9`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/coma/submit-$$i.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" > submit-$$i.sh; \
	sbatch submit-$$i.sh; \
	done

ofp:
	cp ../scripts/ana.sh .
	make clean; make ENV=ofp; \
	@for i in `seq 0 16`; do \
	sed "s/@DATA@/$(DATA)/" ../scripts/ofp/submit-$$i.sh | sed "s/@N@/$(N)/" | sed "s/@G@/$(G)/" > submit-$$i.sh; \
	pjsub submit-$$i.sh; \
	done

check: a.out
	OMP_CANCELLATION=true ./a.out -f data/n32d5k4l228.20170305-fimfgn.edges
	OMP_CANCELLATION=true ./a.out -f data/n256d18k3l218.20170305-fimf4z.edges
	OMP_CANCELLATION=true ./a.out -f data/n576d30k3l214.20170305-fimfev.edges
	OMP_CANCELLATION=true ./a.out -f data/n1344d30k3l248.20170305-fimfao.edges

clean:
	rm -rf *.o *~

