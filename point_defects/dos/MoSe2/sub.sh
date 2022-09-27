#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J BKB
# Queue (Partition):
#SBATCH --partition=medium
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --mem-per-cpu=4GB
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=12:00:00

# Run the program:

module load cmake/3.22     
module load  mkl/2021.3     
module load  intel/21.3.0  
module load  impi/2021.3


export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compiler/2021.3.0/linux/compiler/lib/intel64_lin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64


AIMSPATH=/u/alaa/band/FHIaims/build
EXE=aims.x
export aimsbin=${AIMSPATH}/${EXE}
srun $aimsbin > out_$date



