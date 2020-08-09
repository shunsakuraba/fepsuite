#!/bin/zsh
#QSUB -q _YOURQUEUE_
#QSUB -ug _YOURGROUP_
#QSUB -W 96:00
#QSUB -ry

export ABFE_ROOT=_ABFE_PATH_
source para_conf.zsh

if [[ -z $QSUB_PROCS ]]; then
  # submit itself
  CUR=${PWD##*/}
  ID=$1
  if [[ -n $ID ]]; then
    shift
  fi
  jobs=($@)
  if [[ -z $ID ]] || [[ -z $jobs ]]; then
    echo "Usage: $0 (runpath) (stages)"
    exit 1
  fi
  for i in $jobs; do
    PROCS=36   # Default procs
    TPP=1      # Threads / proc
    CPP=1      # Cores / proc
    CPN=36     # Cores / node
    SYSTEMCPN=36
    DEPENDS=()
    EXTRA=()
    waitcmd=()
    eval $($ABFE_ROOT/pipeline.zsh query $i)
    # edit this qstat analysis part to adapt to other job systems
    qs=$(qstat)
    deps=""
    for d in $DEPENDS; do
        if [[ -e $ID.jobid ]] && grep -q "$CUR.$d\s$ID" $ID.jobid ; then
            res=$(grep "$CUR.$d\s$ID" $ID.jobid | tail -n 1 | cut -f3)
            if grep -q $res <<< $qs; then
                deps="$deps:$res"
            fi
        fi
    done
    if [[ $deps != "" ]]; then
        waitcmd=(-W depend=afterok:"${deps#:}")
    fi
    (( NODES = (PROCS + (CPN/CPP) - 1) / (CPN/CPP) ))
    (( REQPROCS = NODES * (SYSTEMCPN / CPP) ))
    set -e 
    # edit this line to adapt to other job systems
    cmd=(qsub $waitcmd -N ${ID:gs|/|_|} -v "NAME=$CUR.$i,PROCS=$PROCS,CPN=$CPN,ID=$ID" $EXTRA -A p=${REQPROCS}:c=${CPP}:t=${TPP}:m=3413 $0)
    builtin echo $cmd
    # edit this line to adapt to other job systems
    res=$($cmd)
    echo Submitted jobid $res.
    echo "$CUR.$i\t$ID\t$res" >> $ID.jobid
    set +e
  done
  exit $?
fi


STEPNO=${NAME:e}
echo "STEP=$STEPNO"

# Site-specific settings starts from here
cd $QSUB_WORKDIR
module () {
        eval `/opt/modules/default/bin/modulecmd zsh $*`
}

mpirun_() {
    np=$1
    shift
    # edit this part based on environments
    mpiexec.hydra -np $np -perhost $CPN $@
}
# some supercomputer (e.g. Cray KNL) requires extra commands before single-process run. E.g. (aprun -n 1)
SINGLERUN=()
# ---- Site-specific settings ends here

source PATH_TO_GROMACS/bin/GMXRC
PYTHON3=python3

GMX=$(which gmx_mpi)
GMX_MPI=$(which gmx_mpi)

# actual runs

source $ABFE_ROOT/pipeline.zsh run $STEPNO

