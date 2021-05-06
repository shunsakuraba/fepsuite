# GROMACS parallelization settings. 
# For all cases, *_PARA is the number of MPI procs to be run for each replica (or the entire run for pre-run).
# Thus, number of total MPI processes is e.g. NCPLX_RESTR * NCPLX_RESTR_PARA for the restraint replica run.
# If you are using GPU, it is recommended to set all *_PARA = 1.
# Recommended number of _PARA is a divisor of the number of cores / CPU.
# Number of replicas (values without _PARA) should be kept to recommended value (except when ligand is large, where you may need 20 replicas)

# Pre-run nproc
NPRERUN_PARA=36

# Prot-lig restraint. Recommended: 4
NCPLX_RESTR=4
NCPLX_RESTR_PARA=6

# Charge-discharge. Recommended: 12
NCHARGE=12
NCHARGE_PARA=6

# ligand-protein annihilation. Recommended: 12
NANN_PRO=12
NANN_PRO_PARA=6

# ligand annihilation. Recommended:12, PARA is recommended to be <=2 (due to small size)
NANN=12
NANN_PARA=2

# first RUN_PROD ps will be ignored as equilibration
RUN_PROD=2000

#*************************
# System setups
#*************************
# Ligand selection string, uses python "mdtraj" library's selection syntax
# Note "residue" or "resSeq" is PDB-based residue ID, while "resid" is 0-origin residue number from the beginning.
LIG_MDTRAJ="resname MOL"

# Receptor selection string for mdtraj (used to determine restraint)
RECEPTOR_MDTRAJ="protein"

#*************************
# Thresholding
#*************************
# If ligand RMSD during equilibration exceeds this value, calculation stops (unit: nm)
# RMSD is measured by best-fitting the receptor, but only ligand is used in RMSD calculation
# Reference structure used is the final snapshot of the equilibration
EQ_RMSD_CUTOFF=0.4

# Set this value to max (ligand diameter) if calculation fails due to bond length exceeding domain size
# (unit: nm)
MAX_BONDED_INTERACTION_DIST=0

#*************************
# Solvation and ionization
#*************************
# water structure file, tip3p & spc & spc-e -> spc216, tip4p -> tip4p
WATER_STRUCTURE=spc216

# water thickness used in ligand system
WATER_THICKNESS=1.0
WATER_THICKNESS_CHARGING=$WATER_THICKNESS

# ionic strength in M
IONIC_STRENGTH=0.150

# positive ion
ION_POSITIVE=NA
# negative ion
ION_NEGATIVE=CL

#**************************
# default is exit 100 on error; this is necessary for Grid Engine family
ERRORCODE=100
