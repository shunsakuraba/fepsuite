#!/bin/zsh
#PJM -L elapse=9:00:00

# Customization point 1: ^^ elapse time

# Customization point 2: Directory
export FEPREST_ROOT=_PATH_TO_FEPREST_
export PIPELINE=${FEPREST_ROOT}/pipeline.zsh
export GMX_DIR=_PATH_TO_GMX_

source para_conf.zsh

if [[ -z $PJM_O_WORKDIR ]]; then
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
    qs=$(pjstat)
    for i in $jobs; do
        PROCS=1   # Default MPI procs
        GPP=1     # GPUs / proc. This should be 1 or 0 (using more than 2 GPUs per proc is not beneficial in any environment)
        TPP=4     # Default Threads / proc
        CPP=1     # Default Cores / proc
        # For vnode-aware run I set this to 24. 
        SYSTEMCPN=20 # hardware limit of CPU cores / node.
        SYSTEMGPN=4  # hardware limit of GPU / node. 0 if no GPUs
        CPN=$SYSTEMCPN # software limit of CPU / node. Typically, it is not beneficial to exceed SYSTEMCPN.
        GPN=16         # software limit of GPU / node. Often, more than one processes per GPU is beneficial.
        DEPENDS=()
        EXTRA=()
        waitcmd=()
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
            waitcmd=(--sparam "jnam=$ID,sd=ec!=0:after,sn=$i")
        else
            waitcmd=(--sparam "jnam=$ID,sn=$i")
        fi

        if (( PROCS >= 32 )); then
            # for replica run
            CPN=$SYSTEMCPN
            TPP=1
            GPN=16
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

        # set queue. edit this block to your specific envieonments
        if (( PROCS == 1 )); then
            QUEUE=cx-share
            SYSTEMGPN=1
	    EXTRA=(-L gpu=1)
        elif (( NODES <= 8 )); then
            QUEUE=cx-small
            EXTRA=(-L node=$NODES)
        else
            echo "Unimplemented node number (PROCS=$PROCS)"
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
        # edit this line to adapt to other job systems.
        cmd=(pjsub $waitcmd -L rscgrp=$QUEUE -N $(basename $(realpath $PWD/..)).$CUR.$i -x "NAME=$CUR.$i,WID=$ID,STEPNO=$i,PROCS=$PROCS,CPN=$CPN,PPN=$(( (PROCS+NODES-1)/NODES )),GPP=$GPP,REQGPU_PER_NODE=$REQGPU_PER_NODE" --step $EXTRA --mpi proc=$(( int(ceil(PROCS/(NODES + 0.0))) )) $0)
        builtin echo $cmd
        res=$($cmd)
	case $res in
            \[INFO\]\ PJM\ *\ submitted.)
                res=${res%% submitted.}
                res=${res##*Job }
                ;;
            *)
                echo "Unexpected pjsub result: \"$res\""
                exit 1
                ;;
	esac
        echo Submitted jobid $res.
        echo "$CUR.$i\t$ID\t$res" >> $ID.jobid
        set +e
        qs+="$res\n"
    done
    exit $?
fi

echo "STEP=$STEPNO"
export ID=$WID

# Site-specific settings starts from here
cd $PJM_O_WORKDIR
# source again if there is an additional conf file
[[ -e $ID/para_conf.zsh ]] && source $ID/para_conf.zsh

source /etc/profile.d/modules.sh
# use module files according to your compile-time setting, we used below
module load cuda/11.3.1 gcc/8.4.0
module load openmpi_cuda/4.1.2

mpirun_() {
    np=$1
    shift
    # edit this part based on environments
    mpirun -machinefile $PJM_O_NODEINF -npernode $PPN -n $PROCS $@
}
# some supercomputer (e.g. Cray KNL) requires extra commands before single-process run. E.g. (aprun -n 1)
SINGLERUN=(mpirun -n 1)
# ---- Site-specific settings ends here

source $GMX_DIR/bin/GMXRC
PYTHON3=python3

GMX=$(which gmx_mpi)
GMX_MPI=$(which gmx_mpi)

env > $ID/env.txt
# actual runs

setopt ERR_EXIT
source $PIPELINE run $STEPNO
