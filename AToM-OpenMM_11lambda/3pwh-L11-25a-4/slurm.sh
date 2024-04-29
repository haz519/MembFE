#!/bin/bash

/home/haz519/.conda/envs/atm/bin/python3.10 3pwh-L11-25a_mintherm.py
/home/haz519/.conda/envs/atm/bin/python3.10 3pwh-L11-25a_mdlambda.py
echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
/home/haz519/.conda/envs/atm/bin/python3.10 /home/haz519/workdir/ATM/AToM-OpenMM-master/rbfe_explicit.py 3pwh-L11-25a_asyncre.cntl

