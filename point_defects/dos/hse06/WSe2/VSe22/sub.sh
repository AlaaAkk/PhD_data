#!/bin/bash -l
#Standard output and error:
#SBATCH -o ./out.%j
#SBATCH -e ./err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J Alaa
#
#SBATCH --nodes=4            # Request 1 (or more) node(s)
#SBATCH --constraint="gpu"    #    providing GPUs.
#SBATCH --ntasks-per-node=72  # Launch 72 tasks per node
#SBATCH --gres=gpu:a100:4     # Request all 4 GPUs of each node
#SBATCH --nvmps               # Launch NVIDIA MPS to enable concurrent access to the GPUs from multiple processes efficiently
#
#SBATCH --mail-type=all
#SBATCH --mail-user=akkoush@fhi-berlin.mpg.de
#SBATCH --time=24:00:00
export OMP_NUM_THREADS=1
module purge

module load cmake/3.22     
module load mkl/2021.4     
module load intel/21.4.0   
module load impi/2021.4




export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compiler/2021.4.0/linux/compiler/lib/intel64_lin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries/linux/lib/intel64:$MKL_HOME/lib/intel64/:$HOME/.local/lib
#FHI
AIMSPATH=/u/alaa/FHIaims/build
EXE=aims.x


export aimsbin=${AIMSPATH}/${EXE} 
srun $aimsbin > aims.out

