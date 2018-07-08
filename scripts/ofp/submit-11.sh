#!/bin/bash

#------ pjsub option --------#
#PJM -L rscgrp=regular-flat
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=64
#PJM -L elapse=10:00:00
#PJM -g xg18i003
#PJM -j
#------- Program execution -------#

threads=$(nproc --all)
data=@DATA@
N=@N@
G=@G@

for t in 0.928318 0.886124 0.845849 0.807403 0.770706 0.735676
do
for i in $(seq 1 10)
do
    ./a.out -f $data -g $G -n $N -s $i -w $t -c $t > log.$t.$i.txt &

    while [ $(jobs|wc -l) -ge $threads ]; do
	sleep 1
    done
done
done

wait

echo $SECONDS

