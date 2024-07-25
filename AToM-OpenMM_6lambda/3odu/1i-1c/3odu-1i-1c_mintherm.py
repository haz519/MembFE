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

jobname = "3odu-1i-1c"

displ = [ -20.0, -16.0, 100.0 ]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [ 7564, 7565, 7566, 7567, 7568, 7569, 7570, 7571, 7572, 7573, 7574, 7575, 7576, 7577, 7578, 7579, 7580, 7581, 7582, 7583, 7584, 7585, 7586, 7587, 7588, 7589, 7590, 7591, 7592, 7593, 7594, 7595, 7596, 7597, 7598, 7599, 7600, 7601, 7602, 7603, 7604, 7605, 7606, 7607, 7608, 7609, 7610, 7611, 7612, 7613, 7614, 7615, 7616, 7617 ]
lig2_atoms = [ 7618, 7619, 7620, 7621, 7622, 7623, 7624, 7625, 7626, 7627, 7628, 7629, 7630, 7631, 7632, 7633, 7634, 7635, 7636, 7637, 7638, 7639, 7640, 7641, 7642, 7643, 7644, 7645, 7646, 7647, 7648, 7649, 7650, 7651, 7652, 7653, 7654, 7655, 7656, 7657, 7658, 7659, 7660, 7661, 7662, 7663, 7664, 7665, 7666, 7667, 7668 ]
refatoms_lig1 = [ 5, 9, 17 ]
refatoms_lig2 = [ 5, 9, 17 ]
rcpt_cm_atoms = [ 87, 1117, 1167, 1179, 1229, 1396, 1412, 1464, 2555, 2600, 2619, 2629, 2641, 6376, 6922 ]
restrained_atoms = [ 12, 18, 28, 48, 72, 87, 102, 116, 126, 140, 160, 174, 196, 215, 235, 262, 268, 282, 301, 322, 333, 352, 371, 391, 410, 424, 431, 450, 466, 473, 487, 494, 513, 529, 548, 567, 583, 600, 607, 628, 645, 667, 689, 708, 732, 743, 760, 774, 786, 808, 829, 853, 872, 889, 908, 919, 935, 945, 957, 976, 995, 1015, 1031, 1050, 1064, 1091, 1097, 1117, 1141, 1151, 1167, 1179, 1189, 1205, 1215, 1229, 1253, 1274, 1294, 1301, 1315, 1335, 1354, 1364, 1386, 1396, 1412, 1429, 1445, 1464, 1485, 1499, 1515, 1529, 1548, 1569, 1580, 1591, 1607, 1631, 1650, 1669, 1679, 1699, 1718, 1729, 1748, 1760, 1784, 1805, 1824, 1834, 1853, 1869, 1886, 1896, 1910, 1924, 1935, 1952, 1984, 1990, 2014, 2036, 2055, 2074, 2084, 2099, 2121, 2137, 2153, 2174, 2190, 2197, 2213, 2237, 2264, 2270, 2280, 2299, 2318, 2337, 2351, 2378, 2384, 2396, 2416, 2435, 2455, 2465, 2479, 2495, 2506, 2521, 2531, 2543, 2555, 2579, 2600, 2619, 2629, 2641, 2665, 2685, 2714, 2720, 2734, 2746, 2765, 2789, 2805, 2821, 2837, 2857, 2874, 2894, 2911, 2928, 2947, 2964, 2980, 2987, 3006, 3025, 3052, 3058, 3065, 3084, 3100, 3119, 3138, 3149, 3160, 3181, 3192, 3211, 3230, 3249, 3260, 3282, 3301, 3312, 3329, 3340, 3347, 3358, 3372, 3391, 3411, 3426, 3443, 3462, 3486, 3505, 3517, 3532, 3539, 3558, 3582, 3601, 3623, 3642, 3663, 3685, 3697, 3711, 3726, 3733, 3754, 3775, 3789, 3808, 3815, 3834, 3841, 3858, 3877, 3896, 3910, 3932, 3951, 3957, 3968, 3987, 4001, 4011, 4021, 4043, 4054, 4069, 4088, 4100, 4122, 4132, 4151, 4158, 4182, 4196, 4210, 4224, 4231, 4247, 4266, 4280, 4302, 4314, 4329, 4339, 4354, 4376, 4395, 4415, 4429, 4446, 4458, 4474, 4486, 4496, 4506, 4522, 4546, 4553, 4572, 4591, 4615, 4629, 4639, 4661, 4680, 4710, 4716, 4732, 4753, 4765, 4776, 4795, 4807, 4817, 4833, 4857, 4881, 4891, 4901, 4920, 4939, 4953, 4970, 4986, 5006, 5023, 5040, 5047, 5062, 5076, 5083, 5099, 5109, 5116, 5136, 5150, 5164, 5175, 5194, 5218, 5235, 5254, 5271, 5288, 5310, 5334, 5358, 5370, 5385, 5395, 5405, 5421, 5435, 5454, 5464, 5486, 5497, 5521, 5545, 5566, 5580, 5597, 5619, 5625, 5639, 5663, 5673, 5695, 5719, 5735, 5754, 5768, 5782, 5802, 5826, 5840, 5847, 5861, 5885, 5897, 5907, 5928, 5935, 5946, 5968, 5975, 5992, 6009, 6031, 6055, 6077, 6087, 6106, 6128, 6142, 6156, 6172, 6191, 6210, 6229, 6248, 6258, 6278, 6298, 6308, 6319, 6343, 6370, 6376, 6397, 6418, 6437, 6444, 6463, 6474, 6493, 6505, 6516, 6536, 6555, 6574, 6593, 6608, 6627, 6646, 6668, 6685, 6692, 6702, 6717, 6737, 6752, 6766, 6780, 6796, 6813, 6835, 6859, 6878, 6889, 6908, 6922, 6937, 6947, 6966, 6976, 6996, 7016, 7033, 7044, 7055, 7074, 7096, 7102, 7121, 7140, 7161, 7171, 7191, 7210, 7217, 7227, 7249, 7269, 7291, 7305, 7316, 7326, 7343, 7360, 7370, 7389, 7403, 7414, 7421, 7453, 7459, 7478, 7493, 7509, 7528, 7548  ]


#load system
prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds)

#load the ATM Meta Force facility. Among other things the initializer
#sorts Forces into groups 
atm_utils = ATMMetaForceUtils(system)

#Vsite restraints
lig1_cm_atoms = [ 7564, 7565, 7566, 7567, 7568, 7569, 7570, 7571, 7572, 7573, 7574, 7575, 7576, 7577, 7578, 7579, 7580, 7581, 7582, 7583, 7584, 7585, 7586, 7587, 7588, 7589, 7590, 7591, 7592, 7593, 7594, 7595, 7596, 7597, 7598, 7599, 7600, 7601, 7602, 7603, 7604, 7605, 7606, 7607, 7608, 7609, 7610, 7611, 7612, 7613, 7614, 7615, 7616, 7617 ]
lig2_cm_atoms = [ 7618, 7619, 7620, 7621, 7622, 7623, 7624, 7625, 7626, 7627, 7628, 7629, 7630, 7631, 7632, 7633, 7634, 7635, 7636, 7637, 7638, 7639, 7640, 7641, 7642, 7643, 7644, 7645, 7646, 7647, 7648, 7649, 7650, 7651, 7652, 7653, 7654, 7655, 7656, 7657, 7658, 7659, 7660, 7661, 7662, 7663, 7664, 7665, 7666, 7667, 7668 ]
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
