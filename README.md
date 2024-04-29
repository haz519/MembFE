Amber-TI contains files for relative binding free energy simulations:

1. According to the difference between ligands (LA and LB), change the scmask1 and scmask2 (in the .tmpl files) to be consistent with the system.
2. Run 1_setup.sh to generate 12 lambda windows.
3. Submit the L11-25a.sh script to perform AMBER-TI.
4. To extend the simulations, change the jidx to a different number in the submit.sh file. The number of MD steps (nstlim in ti.in.tmpl file) might also need to be adjusted to reach the desired simulation length.
5. Run L11-25a.sh for the TI simulation extension.
6. Run 3_analysis.sh to compute the relative binding free energy and plot convergence

AToM-OpenMM contains files for relative binding free energy simulations:
1. Adjust the values of LAMBDA, LAMBDA1, and LAMBDA2 based on the desired number of lambdas for the calculations.
2. Run slurm.sh to initiate minimization, thermalization, and relaxation process.
3. Then the AToM simulations will be performed.
4. After reaching the desired simulation length, execute ./analyze.sh to compute the relative binding free energy and plot convergence
