#!/bin/bash
#SBATCH -N 1                     # Num of nodes         (<--kaeru)
#SBATCH -n 1                     # Num of MPI processes (<--kaeru)
#SBATCH --ntasks-per-node=1      # Num of MPI processes per node
#SBATCH --ntasks-per-socket=1    # Num of MPI processes per socket
#SBATCH --cpus-per-task=1        # Num of threads per MPI process
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $SLURM_SUBMIT_DIR

data=@DATA@
N=@N@
G=@G@
T=@T@
K=@K@

for i in $(seq 1 $K)
do
./GraphGolf -f $data -g $G -n $N -w $T -c $T -s $i > log.$T.$i.txt
done

echo $SECONDS
