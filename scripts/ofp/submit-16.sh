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

for t in 0.229951 0.219500 0.209523 0.200000 A Y
do
for i in $(seq 1 10)
do
    if [ $t = "A" ]; then
      ./GraphGolf -f $data -g $G -n $N -s $i -w 100 -c 0.5 > log.a.$i.txt &
    elif [ $t = "Y" ]; then
      ./GraphGolf -f $data -g $G -n $N -s $i -y > log.y.$i.txt &
    else
      ./GraphGolf -f $data -g $G -n $N -s $i -w $t -c $t > log.$t.$i.txt &
    fi

    while [ $(jobs|wc -l) -ge $threads ]; do
	sleep 1
    done
done
done

wait

echo $SECONDS

