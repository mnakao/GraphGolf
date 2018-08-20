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
W=@W@
c=@c@
S=@S@
C=@C@

./a.out -f $data -g $G -n $N -w $W -c $c -C $C -s $S > log.C${C}.$S.txt

echo $SECONDS
