#!/bin/bash
#SBATCH --job-name=L11-25a
#SBATCH --partition=im2080-gpu,im1080-gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=error.%j.out
#SBATCH --export=ALL
#SBATCH -t 48:00:00
#SBATCH --oversubscribe

module load anaconda3
module load cuda/11.6.0

export OPENMM_PLUGIN_DIR='/home/haz519/.conda/envs/atm/lib/plugins'
export PYTHONPATH='/home/haz519/.conda/envs/atm/lib/python3.10/site-packages'
export LD_LIBRARY_PATH='/home/haz519/.conda/envs/atm/lib:$LD_LIBRARY_PATH'
export OPENMM_CUDA_COMPILER='/share/Apps/compilers/opt/spack/linux-centos8-x86_64/gcc-8.3.1/cuda/11.6.0-kmhdalj/bin/nvcc'
cd $SLURM_SUBMIT_DIR

CURJOB=TOSED
#/home/haz519/.conda/envs/atm/bin/python3.10 3pwh-L11-25a_mintherm.py
#/home/haz519/.conda/envs/atm/bin/python3.10 3pwh-L11-25a_mdlambda.py
echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
/home/haz519/.conda/envs/atm/bin/python3.10 /home/haz519/workdir/ATM/AToM-OpenMM-master-memb/rbfe_explicit.py 3pwh-L11-25a_asyncre.cntl
NEXTJOB=$((CURJOB + 1))
sed "0,/TOSED/{s/TOSED/${NEXTJOB}/}" slurm.sh > job${NEXTJOB}.sh
sbatch job${NEXTJOB}.sh
