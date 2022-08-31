#!/home/shun/lcl/bin/zsh

export ABFE_ROOT=_ABFE_
export PIPELINE=$ABFE_ROOT/pipeline.zsh
export GMX_DIR=_GMX_DIR_
source para_conf.zsh

if [[ -z $STEPNO ]]; then
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
  # source again if modified by additional conf
  if [[ -e $ID/para_conf.zsh ]]; then
    source $ID/para_conf.zsh
  fi
  if [[ all = $jobs ]]; then
    jobs=($($PIPELINE query all))
  fi
  for i in $jobs; do
    PROCS=36   # Default procs. Currently unused
    eval $($PIPELINE query $i)
    set -e 
    cmd=(NAME=$CUR STEPNO=$i PROCS=$PROCS ID=$ID $0)
    builtin echo $cmd
    echo Running job step $i.
    eval $cmd
    ERR=$?
    if (( ERR > 0 )); then
      echo "Error on step $i"
      exit $?
    fi
    echo Completed job step $i.
  done
  exit $?
fi

echo "STEP=$STEPNO"

if [[ -e $ID/para_conf.zsh ]]; then
  source $ID/para_conf.zsh
fi

# Site-specific settings starts from here

# prepare module command your used during compilation here
# module load foo bar baz...


mpirun_() {
    np=$1
    shift
    # edit this part based on environments
    mpirun -np $np $@
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

