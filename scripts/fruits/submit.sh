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
S=@S@

if [ ${T} = "SA1" ]; then
  ./a.out -f $data -g $G -n $N -w 100 -c 0.2 -s $S > log.SA1.$S.txt
elif [ ${T} = "SA2" ]; then
  ./a.out -f $data -g $G -n $N -w 100 -c 0.2 -s $S -C 10000 > log.SA2.$S.txt
else
  ./a.out -f $data -g $G -n $N -w $T -c $T -s $S > log.$T.$S.txt
fi

echo $SECONDS
