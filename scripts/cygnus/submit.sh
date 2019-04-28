#!/bin/bash
#--------------------------------
#PBS -q fpga
#PBS -A XMPF
#PBS -l elapstim_req=12:00:00
#PBS -v NQSV_MPI_VER=2.3.1/intel
#PBS -b 8
#--------------------------------
cd $PBS_O_WORKDIR
module purge
module load intel/19.0.3 mvapich/2.3.1/intel

data=@DATA@
N=@N@
G=@G@
T=@T@
s=@S@

mpiexec ${NQSII_MPIOPTS} -np 192 ./gg -f $data -n $N -g $G -w $T -s $s -o $G.$s.edges > log.$G.$s.txt



