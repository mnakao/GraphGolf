#!/bin/bash
#SBATCH -J job-name
#SBATCH -N 1                     # Num of nodes         (<--kaeru)
#SBATCH -n 1                     # Num of MPI processes (<--kaeru)
#SBATCH --ntasks-per-node=1      # Num of MPI processes per node
#SBATCH --ntasks-per-socket=1    # Num of MPI processes per socket
#SBATCH --cpus-per-task=1        # Num of threads per MPI process
#SBATCH -t 0:10:00
module purge
module load intel intelmpi mkl
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $SLURM_SUBMIT_DIR

threads=$(nproc --all)
data="../data/n16d5.random.edges"
N=1000000
G=2

for t in 4.954153 4.728979 4.514039 4.308869 4.113025 3.926081 3.747635 3.577299 3.414705 3.259502 
do
for i in $(seq 1 10)
do
    ./GraphGolf -f $data -g $G -n $N -s $i -w $t -c $t > log.$t.$i.txt &

    while [ $(jobs|wc -l) -ge $threads ]; do
	sleep 1
    done
done
done

wait
echo $SECONDS
