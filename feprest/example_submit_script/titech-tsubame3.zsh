#!/bin/zsh
#$ -cwd

export FEPREST_ROOT=_FEPREST_PATH_
export PIPELINE=$FEPREST_ROOT/pipeline.zsh
export GMX_DIR=_PATH_TO_GMX_
source para_conf.zsh

if [[ -z $JOB_NAME ]]; then
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
    PROCS=1     # Default procs
    TPP=1       # Threads / proc
    CPP=1       # Cores / proc
    CPN=4       # Cores / node. Real procs are 28 but we want one GPU per run
    SYSTEMCPN=4 # Cores / node (hardware)
    T=24        # 24 hr run. Modify this part depending on your system size
    NOGPU_STAGE=
    DEPENDS=()
    EXTRA=()
    waitcmd=()
    eval $($PIPELINE query $i)
    # edit this qstat analysis part to adapt to other job systems
    qs=$(qstat)
    deps=""
    for d in $DEPENDS; do
      deps="$deps,FEP$ID.$d"
    done
    if [[ $deps != "" ]]; then
      waitcmd=(-hold_jid ${deps#,})
    fi
    (( NODES = (PROCS + (CPN/CPP) - 1) / (CPN/CPP) ))
    (( REQPROCS = NODES * (SYSTEMCPN / CPP) ))
    # Job queue request setting. Modify here if you want to change.
    case $PROCS in
        1)
            RESOURCE=q_node
            if [[ $NOGPU_STAGE = yes ]]; then
                RESOURCE=s_core
            fi
            ;;
        2)
            RESOURCE=h_node
            ;;
        *)
            RESOURCE=f_node
            ;;
    esac
    NNODE=$NODES
    set -e 
    # edit this line to adapt to other job systems
    cmd=(qsub -g _GROUP_ID_ $waitcmd $EXTRA -v "NAME=$CUR.$i,PROCS=$PROCS,CPN=$CPN,ID=$ID" -l $RESOURCE=$NNODE -l h_rt=$T:00:00 -N "FEP$ID.$i" $0) # T3 only accept names starting from alphabets
    builtin echo $cmd
    # edit this line to adapt to other job systems
    res=$($cmd | head -1 | cut -d ' ' -f 3)
    echo Submitted jobid $res.
    echo "$CUR.$i\t$ID\t$res" >> $ID.jobid
    set +e
  done
  exit $?
fi

# This hack is required because of Gridge Engine's fxxking restriction that errcode must be 100
if [[ -z $OUTSIDE_ERROR_TRAP ]]; then
    OUTSIDE_ERROR_TRAP=yes zsh $0
    if (( $? != 0 )); then
        exit 100
    fi
    exit 0
fi


STEPNO=${NAME:e}
echo "STEP=$STEPNO"

# Site-specific settings starts from here
source /etc/profile.d/modules.sh

# Modify this part based on modules you used during the compilation
module load cuda/9.2.148 intel-mpi/19.0.117

mpirun_() {
    np=$1
    shift
    # Edit this part based on environments
    # If you use openmpi you may need to specify -x LD_LIBRARY_PATH -x OMP_NUM_THREADS
    mpirun -np $np -ppn $CPN $@
}
# some supercomputer (e.g. Cray KNL) requires extra commands before single-process run. E.g. (aprun -n 1)
SINGLERUN=()
# ---- Site-specific settings ends here

source $GMX_DIR/bin/GMXRC.zsh
PYTHON3=python3

GMX=$(which gmx_mpi)
GMX_MPI=$(which gmx_mpi)

# actual runs

setopt ERR_EXIT
source $PIPELINE run $STEPNO
