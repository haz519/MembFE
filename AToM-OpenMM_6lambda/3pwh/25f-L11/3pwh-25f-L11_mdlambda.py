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

jobname = "3pwh-25f-L11"

displ = [ -20.0, -16.0, -40.0 ]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [ 4692, 4693, 4694, 4695, 4696, 4697, 4698, 4699, 4700, 4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716, 4717, 4718, 4719, 4720, 4721, 4722, 4723, 4724 ]
lig2_atoms = [ 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 4734, 4735, 4736, 4737, 4738, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4747, 4748 ]
refatoms_lig1 = [ 32, 8, 5 ]
refatoms_lig2 = [ 23, 8, 5 ]
rcpt_cm_atoms = [ 1116, 1157, 1173, 2394, 2413, 2433, 3669, 3726, 3745, 3800, 3997, 4059, 4126 ]
restrained_atoms = [ 4, 15, 31, 52, 71, 85, 101, 116, 135, 145, 164, 174, 190, 209, 219, 238, 257, 264, 278, 294, 313, 329, 340, 364, 374, 390, 414, 433, 447, 458, 472, 491, 508, 522, 538, 552, 566, 587, 607, 623, 639, 650, 669, 679, 689, 699, 711, 730, 749, 765, 772, 788, 807, 817, 844, 850, 870, 880, 899, 913, 932, 943, 957, 964, 984, 994, 1004, 1014, 1024, 1041, 1048, 1058, 1077, 1097, 1116, 1126, 1137, 1157, 1173, 1192, 1208, 1227, 1237, 1254, 1265, 1276, 1295, 1315, 1326, 1345, 1364, 1374, 1393, 1403, 1422, 1434, 1458, 1479, 1498, 1508, 1527, 1537, 1564, 1570, 1589, 1613, 1634, 1648, 1655, 1674, 1690, 1704, 1711, 1725, 1749, 1759, 1769, 1776, 1795, 1814, 1824, 1843, 1854, 1878, 1894, 1913, 1924, 1944, 1954, 1973, 1980, 1999, 2021, 2027, 2044, 2063, 2070, 2094, 2108, 2122, 2132, 2139, 2164, 2170, 2192, 2207, 2214, 2236, 2250, 2267, 2278, 2295, 2302, 2312, 2319, 2334, 2341, 2358, 2374, 2384, 2394, 2413, 2433, 2448, 2460, 2476, 2500, 2506, 2523, 2537, 2558, 2575, 2591, 2612, 2632, 2646, 2666, 2686, 2696, 2707, 2723, 2742, 2766, 2772, 2791, 2810, 2829, 2846, 2865, 2872, 2888, 2909, 2928, 2952, 2971, 2991, 3001, 3011, 3021, 3045, 3069, 3086, 3105, 3127, 3144, 3161, 3176, 3187, 3212, 3218, 3245, 3251, 3258, 3273, 3297, 3307, 3331, 3342, 3356, 3375, 3392, 3414, 3429, 3445, 3462, 3472, 3482, 3504, 3515, 3525, 3535, 3554, 3573, 3583, 3590, 3609, 3629, 3639, 3658, 3669, 3693, 3720, 3726, 3745, 3762, 3781, 3800, 3814, 3825, 3845, 3859, 3879, 3899, 3917, 3923, 3935, 3945, 3956, 3973, 3991, 3997, 4016, 4040, 4059, 4076, 4097, 4116, 4126, 4145, 4161, 4180, 4190, 4207, 4221, 4235, 4246, 4262, 4278, 4300, 4306, 4326, 4345, 4366, 4376, 4397, 4421, 4440, 4464, 4479, 4499, 4523, 4540, 4554, 4574, 4598, 4620, 4639, 4658, 4682 ]

#load system
prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds)

#load the ATM Meta Force facility. Among other things the initializer
#sorts Forces into groups 
atm_utils = ATMMetaForceUtils(system)

#Vsite restraints
lig1_cm_atoms = [ 4692, 4693, 4694, 4695, 4696, 4697, 4698, 4699, 4700, 4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716, 4717, 4718, 4719, 4720, 4721, 4722, 4723, 4724 ]
lig2_cm_atoms = [ 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 4734, 4735, 4736, 4737, 4738, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4747, 4748 ]
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

