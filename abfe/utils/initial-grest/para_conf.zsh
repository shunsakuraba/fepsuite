# GROMACS parallelization settings. 
# N* parameters are the number of replicas in Hamiltonian replica exchange.
# *_PARA parameters are the number of MPI processes for each replica.
# Thus, number of total MPI processes is e.g. NRESTR * COMPLEX_PARA for the restraint replica run.


# GREST processes
NGREST=8

# Number of MPI processes used when calculating complexes (COMPLEX_PARA).
# Recommended number of _PARA is a divisor of the number of cores / CPU.
# If you use GPU, we recommend 1 for both variables.
COMPLEX_PARA=1

# first RUN_PROD ps will be ignored as equilibration
RUN_PROD=2000

#*************************
# System setups
#*************************
# Ligand selection string, uses python "mdtraj" library's selection syntax
# Note "residue" or "resSeq" is PDB-based residue ID, while "resid" is 0-origin residue number from the beginning.
LIG_GMX="MOL"

# Receptor selection string for mdtraj (used to determine restraint)
RECEPTOR_MDTRAJ="protein"

# gREST temperature for softened atoms
TEMP_HIGH=600.0

#*************************
# Thresholding
#*************************
# Residues within this distance from ligand will be "softened" by gREST
GREST_DISTANCE=0.4

# Conformations over this threshold will be considered as a "result"
CLUSTER_THRESHOLD=0.3

#**************************
# default is exit 100 on error; this is necessary for Grid Engine family
ERRORCODE=100
