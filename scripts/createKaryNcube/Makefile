CC=gcc
all: a.out
a.out: main.o
	$(CC) main.o -o a.out
main.o: main.c
	$(CC) main.c -c
clean:
	rm -f a.out *.o
