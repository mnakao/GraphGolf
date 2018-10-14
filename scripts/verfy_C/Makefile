CC=mpicc
CFLAGS=-std=gnu99 -O2 -Wall -mcmodel=medium -fopenmp $(DEBUG)
SRCS := main.c
OBJS := $(SRCS:%.c=%.o)

all: a.out

a.out: $(OBJS)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $^

%.o: %.c common.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o *~
