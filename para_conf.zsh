


# Put this file on the parent directory.
# The file on the parent directory is sourced by the pipeline, then the file on the child directory is also sourced.

# initial structure and topology
BASECONF=conf_ionized.gro
BASETOP=topol_ionized.top

# change to yes if this mutation changes the charge
CHARGED=no

# Number of replicas to use. Recommended: 32-36 except the cases below
# 48-64 for mutations involving charge change and/or large residue such as F/Y/W
NREP=32

# Number of MPI processes per replica. For GPU, recommended = 1.
PARA=1
#PARA=32

# Number of tuning cycles for fep parameters
NTUNE=5

# Simulation time of one chunk in ps, recommended total time is 4 ns for a good trade-off. Note you can continue the run with additional stages.
SIMLENGTH=4000

# replica exchange interval, recommended: 1000
REPLICA_INTERVAL=1000

# Energy sampling interval, recommended: 100
SAMPLING_INTERVAL=100

# Increase if you are using strange topology file. If you made the reference state by nucfepgen you will need at least 1 for the lack of reference state position restraint.
BASEWARN=1

# reference coordinate upon initialization
REFINIT=$BASECONF
# Set if you want extra restraints during FEP (useful for protein complex)
REFCRD=

# Because dummy atoms in state A does not interact with those in state B, there are exclusions between these atoms,
# which results in relatively long interaction between the two. Option -rdd on mdrun is thus required with this number
DOMAIN_SHRINK=0.6

# Atom types for positive / negative monovalent ions
# AMBER ff series: Na / Cl, CHARMM: SOD/CLA, OPLS: opls_407 / opls_401
AT_POSITIVE=Na
AT_NEGATIVE=Cl

