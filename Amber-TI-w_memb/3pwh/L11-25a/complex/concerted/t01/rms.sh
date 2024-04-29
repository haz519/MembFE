#!/bin/csh
##SBATCH --reservation=woi216_166
#SBATCH --partition=im2080,im1080
#SBATCH --qos=nogpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 48:00:00


cpptraj -i rms_lig.inp
