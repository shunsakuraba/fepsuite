#*************************
# GROMACS parallelization settings 
#*************************
# N* parameters are the number of replicas in Hamiltonian replica exchange.
# *_PARA parameters are the number of MPI processes for each replica.
# Thus, number of total MPI processes is e.g. NRESTR * COMPLEX_PARA for the restraint replica run.


# Prot-lig restraint relaxation. Recommended: 4
NRESTR=4

# Charge-discharge. Recommended: 8-12. For large molecules & strongly charged molecules we recommend greater values.
NCHARGE=12

# Annihilation. Recommended:12. For large molecules we recommend greater values, e.g. for ligand with 60 heavy atoms NANNIH=20 will be a good starting point.
NANNIH=12

# Number of MPI processes used when calculating ligands (LIG_PARA) and complexes (COMPLEX_PARA).
# Recommended number of _PARA is a divisor of the number of cores / CPU.
# LIG_PARA needs to be a small number (due to its small system size).
# If you use GPU, we recommend 1 for both variables.
LIG_PARA=2
COMPLEX_PARA=12

# Threads per process. Recommended values are: 1 (CPU) or (GPU).
TPP=1

# first RUN_PROD ps will be ignored as equilibration
RUN_PROD=2000

#*************************
# System setups
#*************************

# Gromacs [ moleculetype ] name of the ligand.
LIG_GMX="MOL"

# Receptor selection string (used to determine restraints). Python "mdtraj" library's selection syntax.
# Note in the mdtraj library's syntax, "residue" or "resSeq" is PDB-based residue ID (probably 1-origin), while "resid" is 0-origin residue number from the beginning.
RECEPTOR_MDTRAJ="protein"

#*************************
# Thresholding
#*************************
# If ligand RMSD during equilibration exceeds this value, calculation stops (unit: nm)
# RMSD is measured by best-fitting the receptor, but only ligand is used in RMSD calculation
# Reference structure used is the final snapshot of the equilibration
EQ_RMSD_CUTOFF=0.4

# Set this value to max (ligand diameter) if calculation fails due to the bond length exceeding the domain size
# (unit: nm)
MAX_BONDED_INTERACTION_DIST=0

# change this value to max possible ligand diamter (nm). This is usually unnecessary as the proper values are set automatically.
LIGAND_DIAMETER=0

#*************************
# Parameter optimization
#*************************
# Annihilation lambda parameters are optimized during pre-run phase, ANNIH_LAMBDA_OPT_LENGTH ps and ANNIH_LAMBDA_OPT times
ANNIH_LAMBDA_OPT=5
ANNIH_LAMBDA_OPT_LENGTH=50

#*************************
# Solvation and ionization
#*************************
# water structure file, tip3p & spc & spc-e -> spc216, tip4p -> tip4p
WATER_STRUCTURE=spc216

# water thickness used in ligand system
WATER_THICKNESS=1.0

# ionic strength in M
IONIC_STRENGTH=0.150

# positive ion name (CHARMM uses SOD)
ION_POSITIVE=NA
# negative ion name (CHARMM USES CLA)
ION_NEGATIVE=CL
