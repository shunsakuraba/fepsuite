#!/bin/zsh
#$ -cwd

export ABFE_ROOT=_ABFE_PATH_
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
    PROCS=1    # Default procs
    TPP=1      # Threads / proc
    CPP=1      # Cores / proc
    CPN=4      # Cores / node. Real procs are 28 but we want one GPU per run
    SYSTEMCPN=4
    T=24       # 8 hr run. Modify this part depending on your system size
    DEPENDS=()
    EXTRA=()
    waitcmd=()
    eval $($ABFE_ROOT/pipeline.zsh query $i)
    # edit this qstat analysis part to adapt to other job systems
    qs=$(qstat)
    deps=""
    for d in $DEPENDS; do
        if [[ -e $ID.jobid ]] && grep -q "$CUR.$d\s$ID" $ID.jobid ; then
            if grep -q $res <<< $qs; then
                deps="$deps,FEP$ID.$d"
            fi
        fi
    done
    if [[ $deps != "" ]]; then
        waitcmd=(-hold_jid ${deps#,})
    fi
    (( NODES = (PROCS + (CPN/CPP) - 1) / (CPN/CPP) ))
    (( REQPROCS = NODES * (SYSTEMCPN / CPP) ))
    case $PROCS in
        1)
            RESOURCE=q_node
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

STEPNO=${NAME:e}
echo "STEP=$STEPNO"

# Site-specific settings starts from here
source /etc/profile.d/modules.sh

module load cuda/9.2.148 intel-mpi/19.0.117
module load gromacs/2019.4

mpirun_() {
    np=$1
    shift
    # edit this part based on environments
    mpirun -np $np -ppn $CPN $@
}
# some supercomputer (e.g. Cray KNL) requires extra commands before single-process run. E.g. (aprun -n 1)
SINGLERUN=()
# ---- Site-specific settings ends here

PYTHON3=python3
MOL=${ID##*-}

GMX=$(which gmx_mpi)
GMX_MPI=$(which gmx_mpi)

# actual runs

source $ABFE_ROOT/pipeline.zsh run $STEPNO

