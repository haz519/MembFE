from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from atmmetaforce import *

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = "3zpq-L16-L19"

displ = [ -20.0, -18.0, -40.0 ]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [ 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885 ]
lig2_atoms = [ 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 4899, 4900, 4901, 4902, 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916 ]
refatoms_lig1 = [ 29, 17, 15 ]
refatoms_lig2 = [ 27, 17, 15 ]
rcpt_cm_atoms = [ 1304, 1328, 1372, 1384, 1430, 2362, 2644, 2680, 2742, 2763, 2802, 2813, 2862, 3964, 3984, 4358, 4422 ]
restrained_atoms = [ 4, 11, 21, 36, 55, 74, 85, 102, 119, 143, 158, 168, 175, 192, 203, 222, 241, 258, 268, 287, 303, 319, 338, 357, 376, 392, 402, 409, 423, 439, 458, 474, 493, 503, 513, 532, 539, 550, 564, 581, 605, 624, 641, 655, 674, 688, 702, 721, 741, 760, 774, 785, 804, 814, 825, 835, 847, 866, 882, 898, 905, 924, 943, 959, 983, 989, 1009, 1016, 1026, 1040, 1059, 1075, 1091, 1115, 1122, 1136, 1160, 1179, 1203, 1210, 1221, 1241, 1260, 1270, 1285, 1304, 1328, 1342, 1353, 1372, 1384, 1400, 1419, 1430, 1446, 1460, 1470, 1481, 1500, 1515, 1529, 1548, 1559, 1575, 1594, 1604, 1623, 1635, 1659, 1680, 1699, 1709, 1728, 1742, 1761, 1767, 1787, 1811, 1832, 1849, 1860, 1879, 1896, 1910, 1934, 1944, 1968, 1978, 2000, 2016, 2035, 2054, 2065, 2079, 2095, 2119, 2129, 2148, 2159, 2169, 2188, 2204, 2215, 2235, 2262, 2268, 2287, 2304, 2321, 2338, 2362, 2386, 2410, 2422, 2437, 2457, 2463, 2480, 2490, 2509, 2531, 2541, 2562, 2579, 2599, 2605, 2612, 2622, 2632, 2644, 2664, 2680, 2694, 2708, 2732, 2742, 2763, 2773, 2792, 2802, 2813, 2824, 2843, 2862, 2873, 2893, 2914, 2941, 2947, 2966, 2985, 3004, 3021, 3040, 3060, 3076, 3086, 3105, 3129, 3145, 3166, 3190, 3205, 3215, 3237, 3252, 3269, 3288, 3312, 3334, 3353, 3365, 3389, 3399, 3410, 3432, 3456, 3478, 3492, 3503, 3527, 3543, 3560, 3579, 3596, 3620, 3635, 3652, 3674, 3684, 3703, 3725, 3739, 3758, 3765, 3784, 3803, 3820, 3827, 3843, 3863, 3877, 3896, 3907, 3931, 3958, 3964, 3984, 4004, 4023, 4039, 4053, 4072, 4088, 4102, 4118, 4138, 4152, 4176, 4188, 4207, 4231, 4237, 4249, 4273, 4292, 4312, 4328, 4338, 4358, 4372, 4396, 4415, 4422, 4443, 4453, 4467, 4478, 4488, 4505, 4527, 4533, 4552, 4571, 4592, 4603, 4627, 4646, 4652, 4664, 4684, 4708, 4730, 4740, 4760, 4782, 4806, 4825, 4844 ]

#load system
prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds)

#load the ATM Meta Force facility. Among other things the initializer
#sorts Forces into groups 
atm_utils = ATMMetaForceUtils(system)

#Vsite restraints
lig1_cm_atoms = [ 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885 ]
lig2_cm_atoms = [ 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 4899, 4900, 4901, 4902, 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916 ]
kf = 25.0 * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
r0 = 5 * angstrom #radius of Vsite sphere
atm_utils.addRestraintForce(lig_cm_particles = lig1_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig1_restr_offset)

atm_utils.addRestraintForce(lig_cm_particles = lig2_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig2_restr_offset)

#alignment restraint
lig1_ref_atoms  = [ refatoms_lig1[i]+lig1_atoms[0] for i in range(3)]
lig2_ref_atoms  = [ refatoms_lig2[i]+lig2_atoms[0] for i in range(3)]
atm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                            ligb_ref_particles = lig2_ref_atoms,
                            kfdispl =   2.5 * kilocalorie_per_mole/angstrom**2,
                            ktheta =   25.0 * kilocalorie_per_mole,
                            kpsi =     25.0 * kilocalorie_per_mole,
                            offset = lig2_restr_offset)

#restrain selected atoms
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 1.5 * angstrom
atm_utils.addPosRestraints(restrained_atoms, inpcrd.positions, fc, tol)

#define the thermodynamic/alchemical state
#the system is prepared at lambda = 0 and driven to lambda=1/2
temperature = 300.0 * kelvin
lmbd = 0.0
lambda1 = lmbd
lambda2 = lmbd
alpha = 0.0 / kilocalorie_per_mole
u0 = 0.0 * kilocalorie_per_mole
w0coeff = 0.0 * kilocalorie_per_mole
umsc =  1000.0 * kilocalorie_per_mole
ubcore = 500.0 * kilocalorie_per_mole
acore = 0.062500
direction = 1

#create ATM Force (direction is 1 by default)
atmforcegroup = 2
nonbonded_force_group = 1
atm_utils.setNonbondedForceGroup(nonbonded_force_group)
atmvariableforcegroups = [nonbonded_force_group]
atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore, direction, atmvariableforcegroups )
#adds all atoms to the force with zero displacement
for at in prmtop.topology.atoms():
    atmforce.addParticle(at.index, 0., 0., 0.)
#the ligand atoms get displaced, ligand 1 from binding site to the solvent bulk
#and ligand 2 from the bulk solvent to the binding site
for i in lig1_atoms:
    atmforce.setParticleParameters(i, i, displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)
for i in lig2_atoms:
    atmforce.setParticleParameters(i, i, -displ[0] * angstrom, -displ[1] * angstrom, -displ[2] * angstrom)
atmforce.setForceGroup(atmforcegroup)
system.addForce(atmforce)

#add barostat
#barostat = MonteCarloBarostat(1*bar, temperature)
barostat = MonteCarloMembraneBarostat(1*bar,0*bar*nanometer,300*kelvin,MonteCarloMembraneBarostat.XYIsotropic,MonteCarloMembraneBarostat.ZFree,999999999)
#barostat.setFrequency(999999999)#disabled
system.addForce(barostat)

temperature = 300 * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
integrator = MTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, [(0,1), (atmforcegroup,1)])
integrator.setConstraintTolerance(0.00001)

#platform_name = 'OpenCL'
platform_name = 'CUDA'
platform = Platform.getPlatformByName(platform_name)
properties = {}
properties["Precision"] = "mixed"

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

state = simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup})
pote = state.getPotentialEnergy()

print( "LoadState ...")
simulation.loadState(jobname + '_equil.xml')

#override ATM parameters
simulation.context.setParameter(atmforce.Lambda1(), lambda1)
simulation.context.setParameter(atmforce.Lambda2(), lambda2)
simulation.context.setParameter(atmforce.Alpha(), alpha *kilojoules_per_mole)
simulation.context.setParameter(atmforce.U0(), u0 /kilojoules_per_mole)
simulation.context.setParameter(atmforce.W0(), w0coeff /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Umax(), umsc /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Ubcore(), ubcore /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Acore(), acore)
simulation.context.setParameter(atmforce.Direction(), direction)

state = simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup})
print("Potential Energy =", state.getPotentialEnergy())

print("Annealing to lambda = 1/2 ...")

stepId = 1000
totalSteps = 250000
loopStep = int(totalSteps/stepId)
simulation.reporters.append(StateDataReporter(stdout, stepId, step=True, potentialEnergy = True, temperature=True))
simulation.reporters.append(DCDReporter(jobname + ".dcd", stepId))

binding_file = jobname + '_mdlambda.out'
f = open(binding_file, 'w')

deltalambda = (0.5 - 0.0)/float(loopStep)

for i in range(loopStep):
    simulation.step(stepId)
    state = simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup})
    pot_energy = (state.getPotentialEnergy()).value_in_unit(kilocalorie_per_mole)
    pert_energy = (atmforce.getPerturbationEnergy(simulation.context)).value_in_unit(kilocalorie_per_mole)
    l1 = simulation.context.getParameter(atmforce.Lambda1())
    l2 = simulation.context.getParameter(atmforce.Lambda2())
    a = simulation.context.getParameter(atmforce.Alpha()) / kilojoules_per_mole
    umid = simulation.context.getParameter(atmforce.U0()) * kilojoules_per_mole
    w0 = simulation.context.getParameter(atmforce.W0()) * kilojoules_per_mole
    print("%f %f %f %f %f %f %f %f %f" % (temperature/kelvin,lmbd, l1, l2, a*kilocalorie_per_mole, umid/kilocalorie_per_mole, w0/kilocalorie_per_mole, pot_energy, pert_energy), file=f )
    f.flush()
    lmbd += deltalambda
    lambda1 += deltalambda
    lambda2 += deltalambda
    simulation.context.setParameter(atmforce.Lambda1(), lambda1)
    simulation.context.setParameter(atmforce.Lambda2(), lambda2)

print( "SaveState ...")
simulation.saveState(jobname + "_0.xml")

#save a pdb file for visualization
positions = simulation.context.getState(getPositions=True).getPositions()
with open(jobname + '_0.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")

