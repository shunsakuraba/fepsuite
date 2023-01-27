#!/bin/zsh
#PBS -l walltime=9:00:00
#PBS -m n
# User customize point 1: change ^^: default queue submission settings

# User customize point 2: change paths here
export ABFE_ROOT=_PATH_TO_ABFE_
export PIPELINE=$ABFE_ROOT/pipeline.zsh
export GMX_DIR=_GMX_DIR_

# Current directory may be either $PWD or $PBS_O_WORKDIR
[[ -z $PBS_O_WORKDIR ]] && source para_conf.zsh || source $PBS_O_WORKDIR/para_conf.zsh

if [[ -z $PBS_O_WORKDIR ]]; then
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
    if [[ -e $ID/para_conf.zsh ]]; then
        source $ID/para_conf.zsh
    fi
    if [[ all = $jobs ]]; then
        jobs=($($PIPELINE query all))
    fi
    zmodload zsh/mathfunc
    qs=$(qstat)
    for i in $jobs; do
        PROCS=1   # Default MPI procs
        GPP=1     # GPUs / proc. This should be 1 or 0 (using more than 2 GPUs per proc is not beneficial in any environment)
        TPP=3     # Default Threads / proc
        CPP=1     # Default Cores / proc
        # For vnode-aware run I set this to 24. 
        SYSTEMCPN=24 # hardware limit of CPU cores / node.
        SYSTEMGPN=4  # hardware limit of GPU / node. 0 if no GPUs
        CPN=$SYSTEMCPN # software limit of CPU / node. Typically, it is not beneficial to exceed SYSTEMCPN.
        GPN=8          # software limit of GPU / node. Often, more than one processes per GPU is beneficial.
        DEPENDS=()
        EXTRA=()
        waitcmd=(-W block=true)
        eval $($PIPELINE query $i)
        # Sets dependency job id.
        # Edit this qstat analysis part to adapt to other job systems
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


        (( CPUNODES = int(ceil(PROCS * TPP / (CPN + 0.0))) ))
        if (( SYSTEMGPN > 0 )); then
            (( GPUNODES = int(ceil(PROCS * GPP / (GPN + 0.0))) ))
        else
            GPUNODES=0
        fi

        if (( CPUNODES > GPUNODES )); then
            NODES=$CPUNODES
        else
            NODES=$GPUNODES
        fi

        # User customize point 3: edit this part based on your queue setting.
        # 
        if (( PROCS == 1 )); then
            QUEUE=gpu1
            SYSTEMGPN=1
        elif (( NODES <= 16 )); then
            QUEUE=gpu16
        else
            echo "Unimplemented node number (PORCS=$PROCS)"
            exit 1
        fi

        (( REQPROCS = NODES * SYSTEMCPN ))
        (( REQGPU_PER_NODE = GPP * CPN / CPP ))
        if (( REQGPU_PER_NODE > SYSTEMGPN )); then
            # overbooking: more than one process per GPU
            # pass SYSTEMGPN to queue system, but inside the pipeline we still use requested number of GPUs
            REQGPU_PER_NODE=$SYSTEMGPN
        fi
        set -e 
        # edit this line to adapt to other job systems. -v parts are used in pipeline software. -l select part should be modified for your running environments.
        cmd=(qsub $waitcmd -q $QUEUE -N $(basename $(realpath $PWD/..)).$CUR.$i -v "NAME=$CUR.$i,ID=$ID,STEPNO=$i,PROCS=$PROCS,CPN=$CPN,PPN=$(( (PROCS+NODES-1)/NODES )),GPP=$GPP,REQGPU_PER_NODE=$REQGPU_PER_NODE" $EXTRA $PROJ -l select=$NODES:ncpus=$CPN:ngpus=$REQGPU_PER_NODE:ompthreads=$TPP:mpiprocs=$(( int(ceil(PROCS/(NODES + 0.0))) )) $0)
        builtin echo $cmd
        res=$($cmd)
        echo Submitted jobid $res.
        echo "$CUR.$i\t$ID\t$res" >> $ID.jobid
        set +e
        qs+="$res\n"
    done
    exit $?
fi

echo "STEP=$STEPNO"

# In PBS-based systems cd to workdir first
cd $PBS_O_WORKDIR
# Site-specific settings starts from here
# source again if there is an additional conf file
[[ -e $ID/para_conf.zsh ]] && source $ID/para_conf.zsh

# User customize point 3: load module and libraries if necessary


mpirun_() {
    np=$1
    shift
    # User customize point 4: edit this part based on your environments
    mpirun -perhost $PPN -np $PROCS $@
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
