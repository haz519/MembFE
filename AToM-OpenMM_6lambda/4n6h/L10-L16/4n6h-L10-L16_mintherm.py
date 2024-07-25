from __future__ import print_function

import openmm as mm
from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from atmmetaforce import *

#the multiple-time step integrator does not have a setTemperature() method
def setTemperature(self, temperature):
    self.setGlobalVariableByName('kT', MOLAR_GAS_CONSTANT_R*temperature)
MTSLangevinIntegrator.setTemperature = setTemperature

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = "4n6h-L10-L16"

displ = [ -20.0, -16.0, -50.0 ]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [ 4842, 4843, 4844, 4845, 4846, 4847, 4848, 4849, 4850, 4851, 4852, 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 4899, 4900, 4901, 4902 ]
lig2_atoms = [ 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4937, 4938, 4939, 4940, 4941, 4942, 4943, 4944, 4945, 4946, 4947, 4948, 4949, 4950, 4951, 4952, 4953, 4954, 4955, 4956, 4957, 4958, 4959, 4960, 4961, 4962, 4963, 4964, 4965, 4966, 4967, 4968, 4969 ]
refatoms_lig1 = [ 8, 25, 47 ]
refatoms_lig2 = [ 8, 25, 47 ]
rcpt_cm_atoms = [ 1026, 1064, 1086, 1397, 1409, 2764, 3904, 4235, 4282 ]
restrained_atoms = [ 4, 11, 30, 36, 43, 53, 77, 88, 98, 109, 120, 139, 149, 168, 178, 197, 207, 226, 240, 250, 269, 290, 301, 311, 327, 338, 348, 364, 371, 390, 409, 416, 430, 446, 465, 481, 498, 518, 525, 544, 560, 584, 605, 619, 641, 658, 680, 694, 704, 718, 732, 751, 772, 791, 811, 825, 844, 854, 873, 883, 895, 905, 924, 934, 948, 959, 973, 1000, 1006, 1026, 1043, 1054, 1064, 1086, 1107, 1126, 1143, 1158, 1172, 1204, 1210, 1230, 1237, 1252, 1271, 1290, 1300, 1322, 1332, 1348, 1367, 1378, 1397, 1409, 1430, 1451, 1465, 1482, 1502, 1516, 1527, 1546, 1566, 1580, 1599, 1613, 1630, 1647, 1658, 1674, 1686, 1710, 1731, 1750, 1760, 1776, 1787, 1812, 1818, 1834, 1856, 1866, 1885, 1897, 1917, 1941, 1963, 1969, 1979, 2001, 2011, 2033, 2052, 2071, 2085, 2104, 2115, 2134, 2158, 2174, 2193, 2203, 2214, 2221, 2237, 2244, 2268, 2274, 2293, 2310, 2326, 2343, 2353, 2369, 2383, 2415, 2421, 2445, 2457, 2464, 2474, 2490, 2506, 2516, 2533, 2552, 2569, 2597, 2603, 2622, 2628, 2639, 2663, 2684, 2708, 2720, 2734, 2750, 2764, 2786, 2805, 2816, 2832, 2852, 2871, 2891, 2901, 2921, 2937, 2961, 2967, 2986, 3005, 3024, 3043, 3057, 3073, 3084, 3105, 3112, 3131, 3148, 3167, 3186, 3210, 3229, 3253, 3264, 3280, 3304, 3323, 3342, 3353, 3360, 3371, 3393, 3408, 3430, 3442, 3466, 3477, 3496, 3520, 3544, 3563, 3577, 3601, 3618, 3634, 3653, 3669, 3685, 3701, 3708, 3718, 3738, 3754, 3770, 3781, 3805, 3823, 3829, 3848, 3865, 3884, 3904, 3920, 3939, 3955, 3979, 3993, 4012, 4028, 4040, 4059, 4071, 4095, 4119, 4139, 4145, 4164, 4180, 4196, 4206, 4216, 4235, 4252, 4271, 4282, 4301, 4311, 4330, 4337, 4358, 4368, 4382, 4393, 4404, 4423, 4445, 4451, 4467, 4486, 4507, 4517, 4537, 4556, 4568, 4583, 4597, 4617, 4639, 4663, 4674, 4694, 4718, 4735, 4754, 4765, 4789, 4819, 4825, 4836  ]


#load system
prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds)

#load the ATM Meta Force facility. Among other things the initializer
#sorts Forces into groups 
atm_utils = ATMMetaForceUtils(system)

#Vsite restraints
lig1_cm_atoms = [ 4842, 4843, 4844, 4845, 4846, 4847, 4848, 4849, 4850, 4851, 4852, 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 4899, 4900, 4901, 4902 ]
lig2_cm_atoms = [ 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4937, 4938, 4939, 4940, 4941, 4942, 4943, 4944, 4945, 4946, 4947, 4948, 4949, 4950, 4951, 4952, 4953, 4954, 4955, 4956, 4957, 4958, 4959, 4960, 4961, 4962, 4963, 4964, 4965, 4966, 4967, 4968, 4969 ]
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
                            kfdispl = 2.5 * kilocalorie_per_mole/angstrom**2,
                            ktheta =  25.0 * kilocalorie_per_mole,
                            kpsi =  25.0 * kilocalorie_per_mole,
                            offset = lig2_restr_offset)

#restrain selected atoms
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 1.5 * angstrom
atm_utils.addPosRestraints(restrained_atoms, inpcrd.positions, fc, tol)

#Set up Langevin integrator
initial_temperature = 50 * kelvin
final_temperature = 300 * kelvin
temperature = initial_temperature
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
#barostat = MonteCarloBarostat(1*bar, final_temperature)
barostat = MonteCarloMembraneBarostat(1*bar, 0*bar*nanometer,300*kelvin,MonteCarloMembraneBarostat.XYIsotropic,MonteCarloMembraneBarostat.ZFree)
saved_barostat_frequency = barostat.getFrequency()
barostat.setFrequency(999999999)#disabled
system.addForce(barostat)

nonbonded_force_group = 1
atm_utils.setNonbondedForceGroup(nonbonded_force_group)

integrator = MTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, [(0,1), (1,1)])
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

print("Potential energy before minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("Energy minimizing the system ...")
simulation.minimizeEnergy()

print("Potential energy after minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

#saves checkpoint
simulation.saveState(jobname + '_min.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_min.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

print("Thermalization ...")

totalSteps = 150000
steps_per_cycle = 5000
number_of_cycles = int(totalSteps/steps_per_cycle)
delta_temperature = (final_temperature - initial_temperature)/number_of_cycles
simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, volume=True))

#MD with temperature ramp
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)
    #prepare system for new temperature
    temperature = temperature + delta_temperature
    integrator.setTemperature(temperature)

#saves checkpoint
simulation.saveState(jobname + '_therm.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_therm.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)
    
print("NPT equilibration ...")
barostat.setFrequency(saved_barostat_frequency)#enabled

#MD at constant pressure
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)

#saves checkpoint
simulation.saveState(jobname + '_npt.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_npt.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

print("NVT equilibration ...")
barostat.setFrequency(999999999)#disabled

#MD at constant volume
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)

#saves checkpoint
simulation.saveState(jobname + '_equil.xml')
#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_equil.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)
    
end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
