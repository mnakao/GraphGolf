CC=gcc
all: mesh torus
mesh: mesh.o
	$(CC) mesh.o -o mesh
mesh.o: mesh.c
	$(CC) mesh.c -c
torus: torus.o
	$(CC) torus.o -o torus
torus.o: torus.c
	$(CC) torus.c -c
clean:
	rm -f torus mesh *.o
