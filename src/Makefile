.PHONY: clean coma
CC=mpicc
ifeq ($(CC), mpicc)
CFLAGS=-std=gnu99 -O3 -Wall -mcmodel=medium -fopenmp
LFLAGS=-lm
else ifeq ($(CC), mpiicc)
CFLAGS=-O3 -qopenmp -std=gnu99
LFLAGS=
endif
SRCS := main.c sa.c timer.c top_down_bfs.c common.c
OBJS := $(SRCS:%.c=%.o)

all: GraphGolf

GraphGolf: $(OBJS)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $^

%.o: %.c common.h
	$(CC) $(CFLAGS) -c $<

run: GraphGolf
	nohup bash run.sh &

coma: 
	@for i in `seq 0 9`; do \
	sbatch ../scripts/coma/submit-$$i.sh; \
	done

check: GraphGolf
	OMP_CANCELLATION=true ./GraphGolf -f data/n32d5k4l228.20170305-fimfgn.edges
	OMP_CANCELLATION=true ./GraphGolf -f data/n256d18k3l218.20170305-fimf4z.edges
	OMP_CANCELLATION=true ./GraphGolf -f data/n576d30k3l214.20170305-fimfev.edges
	OMP_CANCELLATION=true ./GraphGolf -f data/n1344d30k3l248.20170305-fimfao.edges

clean:
	rm -rf *.o *~

