#! /bin/bash -l

##SBATCH -o ./gpu_test.out.%j
##SBATCH -e ./gpu_test.err.%j
#SBATCH -J MoS2_cluster
#SBATCH --partition=p.ada
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --gpus-per-node=a100:4
#SBATCH --nvmps
#SBATCH --time=2:00:00
export OMP_NUM_THREADS=1


# CHECK THAT THESE MATCH THE MODULE VERSIONS YOU COMPILED WITH
module load intel/21.3.0 
module load mkl/2021.3 
module load impi/2021.3 
module load cuda/11.6
module load anaconda/3/2021.11 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compiler/2021.3.0/linux/compiler/lib/intel64_lin/:$MKL_HOME/lib/intel64/:$HOME/.local/lib:$CUDA_HOME/lib64

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

AIMS=/u/alaa/local/build/aims.x 

srun $AIMS > aims.out 

