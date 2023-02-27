#!/usr/bin/zsh

# Root directory of FEP-suite
FEPSUITE_ROOT=/path/to/fepsuite

# Root GROMACS directory. Note ABFE requires 2022.5 and later
GROMACS_DIR=/path/to/gromacs-2022.5
# uncomment these lines to manually specify the GROMACS binary
#GMX=$GROMACS_DIR/bin/gmx_mpi
#GMX_MPI=$GROMACS_DIR/bin/gmx_mpi

# Specify jobtype, currently abfe or feprest
JOBTYPE=abfe

# Specify jobsystem
JOBSYSTEM=none
# If jobsystem need additional information
#SUBSYSTEM=

# Run actual controller
source $FEPSUITE_ROOT/controller.zsh $0 $@
