# GROMACS parallelization settings. 
# For all cases, *_PARA is the number of MPI procs to be run for each replica (or the entire run for pre-run).
# Thus, number of total MPI processes is e.g. NCPLX_RESTR * NCPLX_RESTR_PARA for the restraint replica run.
# If you are using GPU, it is recommended to set all *_PARA = 1.
# Recommended number of _PARA is a divisor of the number of cores / CPU.
# Number of replicas should be kept to recommended value except special cases that tweaking number of replicas be beneficial.

# Pre-run nproc
NPRERUN_PARA=36

# Prot-lig restraint. Recommended: 8
NCPLX_RESTR=8
NCPLX_RESTR_PARA=6

# Charge-discharge. Recommended: 20
NCHARGE=20
NCHARGE_PARA=6

# ligand-protein annihilation. Recommended: 20
NANN_PRO=20
NANN_PRO_PARA=6

# ligand annihilation. Recommended:20, PARA is recommended to be <=8 (due to small size)
NANN=20
NANN_PARA=6

# first RUN_EQ ps will be ignored as equilibration
RUN_PROD=2000

#*************************
# System setups
#*************************
# Ligand name in GROMACS topology
LIG_TOP="MOL"

# Ligand selection string, uses python "mdtraj" library's selection syntax
LIG_MDTRAJ="resname $LIG_TOP"

# Receptor selection string in GMX ndx (only used to unpbc protein complex)
# If ndx file needs to be regenerated, implement generate_ndx function in this conf file.
RECEPTOR_NDX="Protein"

# Receptor selection string for mdtraj (used to determine restraint)
RECEPTOR_MDTRAJ="protein"

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

