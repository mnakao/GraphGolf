# Dependent software
* C compiler (gcc, intel, and so on)
* MPI library (mpich, openmpi, and so on)

# Options
* -f : Input graph (When using a simple graph, add -N option)
* -o : Output graph
* -s : Random seed (>= 0, integer, default : 0)
* -n : Number of evaluations (> 0, integer, default : 10,000)
* -w : Maximum temperature (> 0, real, default : 100.0)
* -c : Minimum temperature (> 0, real, default : 0.217147, which is the value that accepts the smallest alteration in SA with a probability of 0.01%)
* -g : Number of groups (> 0, integer, default : 1)
* -v : Number of vertices added to the input graph (> 0, integer, default : 0)
       With this option, the degree of the graph increases by one
* -h : Output help
* -B : Use breadth-first search for APSP calculation
* -D : Estimate optimal maximum and minimum temperatures
* -H : Hill climbing method is used instead of SA
* -L : Use an adjacency matrix with saving memory for APSP calculation
* -M : Assume that the input graph is symmetric
* -N : Support not-simple graph
* -O : Omits graph verification that occurs before and after the execution of SA

# Examples of execution
    $ mpirun -np 1 ./gg -f ./data/n16d4.random.edges

    $ mpirun -np 1 ./gg -f ./data/n200d29.random.edges -g 15 -v 19 -o n3019d30.edges

* While the number of vertices in the input graph is 200, the degree in the graph is 29
* Since the number of groups is 15, and 19 vertices are added, the number of vertices in the output graph is 3019
* Also, since the degree is increased by 1, the degree of the output graph is 30

To deal with the graph output above as an input graph, add the "-M" option.

    $ mpirun -np 1 ./gg -f n3019d30.edges -g 15 -v 19 -M -o n3019d30-2.edges
