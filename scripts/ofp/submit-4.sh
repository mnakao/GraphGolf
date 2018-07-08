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

for t in 6.549098 6.251432 5.967294 5.696072 5.437176 5.190048
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

