#!/usr/bin/zsh

set -x

# Root directory of FEP-suite
FEPSUITE_ROOT=$HOME/work/fepsuite

# Root GROMACS directory. Note ABFE and FEP/REST currently need different versions of GROMACS
GROMACS_DIR=$HOME/opt/gromacs-2020-hrexpatch

# Specify jobtype, currently abfe or feprest
JOBTYPE=feprest

# Specify jobsystem
JOBSYSTEM=none
# If jobsystem need additional information fill variables
#SUBSYSTEM=

# Run actual controller
source $FEPSUITE_ROOT/controller.zsh $0 $@
