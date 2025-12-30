# Hokkaido University GrandChariot2
# Combination of PBS Pro + OpenMPI
# How to compile:
# (1) load modules one-by-one (because gcc/11 enables ompi-cuda, load in this order)
# % module load gcc/11 cuda/12.8.0
# % module load ompi-cuda/4.1.6-12.8.0
# % module list
# Currently Loaded Modulefiles:
# 1) gcc/11   2) cuda/12.8.0   3) ompi-cuda/4.1.6-12.8.0
#
# (2) compile
# % ((( download and untargz GROMACS. For FEP/REST, also apply the patch )))
# % cd PATH_TO_GROMACS
# % mkdir build; cd build
# % cmake -DCMAKE_INSTALL_PREFIX=WHERE_TO_INSTALL_GROMACS -DGMX_GPU=CUDA -DGMX_CUDA_TARGET_COMPUTE="90" -DGMX_BUILD_OWN_FFTW -DGMX_MPI=on  ..
# % make -j8 && make install
# (3) install pip-packages
# (4) update run.zsh path
# (5) update run.zsh with
# SUBSYSTEM=GPU
# GROUP=c000000 <- your group ID
#

job_prelude() {
    cd $PBS_O_WORKDIR
    source /etc/profile.d/modules.sh
    module purge
    module load gcc/11
    module load cuda/12.8.0
    module load ompi-cuda/4.1.6-12.8.0
    PYTHON3=python3
    export OMP_NUM_THREADS=$TPP
}

job_prelude_after_gmx() {
    # do nothing
}

job_init_queue_stat() {
    # Load "live" job id list.
    # If we use "dead" ID number in -W depend=afterok:(ID), the simulation silently fails. 
    # Live ID list should be global because we want to use this throughout the process
    typeset -ag LIVE_ID_LIST
    LIVE_ID_LIST=($(qstat | cut -d' ' -f1 | tail -n +3))
}

job_set_preferred_resource() {
    # Although we can separate 1 node into 2 vnodes, it seems not worth the hassle
    case $SUBSYSTEM in
        CPU)
            HW_CPN=64 # 2 sockets, 32 cores/socket
            HW_GPN=0
            CPN=64
            GPN=0
            GPP=0
            ;;
        GPU)
            HW_CPN=64 # 2 sockets, 32 cores/socket
            HW_GPN=4
            CPN=64
            GPN=16
            GPP=1
            # ignores CPU_ONLY_STAGE. This is because in Grand Chariot 2 CPU / GPU points are separate.
            ;;
        *)
            echo "Specify SUBSYSTEM in run.zsh" 1>&2
            exit 1
            ;;
    esac
}

job_submit() {
    # Export variables
    local -a exports
    exports=(PROCS TPP PPN GPP ID STEPNO)
    local -a key_vars
    key_vars=()
    for e in $exports; do
        key_vars+="$e=${(P)e}"
    done

    # Default time
    local timelimit="8:00:00"

    # Set dependency.
    local -a deps
    deps=()
    for d in $DEPENDS; do
        local jid
        jid=$(controller_get_jobid $d)
        # -W depend=afterok:(ID) only accepts live ID in the jobsystem. Few days after the job completion the ID is invalidated.
        # Note in the pathological case, the ID may turn dead while running this script, but we can't do anything (if anyone knows the way please teach me!)
        # Expansion I: index of last occurence (0 not found)
        # Expansion e: exact match
        if (( $LIVE_ID_LIST[(Ie)$jid] )); then
            deps+=$jid
        fi
    done
    local -a dependency
    if (( ${#deps} > 0 )); then
        dependency=(-W depend=afterok:${(j|:|)deps})
    else
        dependency=()
    fi

    # set queue
    local queue
    local -a extra
    extra=(-m n) # disables mail notification
    if (( GPP > 0 )); then
        if (( JOB_NODES <= 1 )); then
            queue=sg
        elif (( JOB_NODES <= 8 )); then
            queue=lg
        else
            echo "Too large number of nodes (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
            exit 1
        fi
    else
        if (( JOB_NODES <= 32 )); then
            queue=sc
        elif (( JOB_NODES <= 256 )); then
            queue=lc
        else
            echo "Too large number of nodes (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
            exit 1
        fi
    fi

    if [[ -z $GROUP ]]; then
        echo "GROUP is not set. Run \"show_token\" first, then set GROUP in run.zsh according to the appropriate group ID." 1>&2
        exit 1
    fi

    local cpugpureq
    if [[ $SUBSYSTEM = GPU ]]; then
        cpugpureq=":ngpus=$JOB_GPN"
    else
        cpugpureq=":nsockets=2"
    fi

    # Run actual code and get jobid
    local -a cmd
    # select: num nodes
    # mpiprocs: mpi procs per node
    # ompthreads: threads per mpi proc
    cmd=(qsub -q $queue -v ${(j:,:)key_vars} -N $JOB_NAME $dependency -l select=$JOB_NODES:mpiprocs=$JOB_PPN:ompthreads=$TPP$cpugpureq -l walltime=$timelimit -W group_list=$GROUP $extra $BASEFILE)
    echo $cmd
    # PBS Pro only emits job id like "12345.sjms", that's great
    JOBID=$($cmd)
    if [[ $? != 0 ]]; then
        echo "Submission failed" 1>&2
        exit 1
    fi
    LIVE_ID_LIST+=$JOBID
}

job_mpirun() {
    local N=$1
    shift

    # MPI_PROC_PER_NODE should be correct even if N is different from total MPI process
    mpirun -n $N --map-by ppr:${MPI_PROC_PER_NODE}:node $@
}

job_singlerun() {
    $@
}

job_get_mode() {
    if [[ -n $PBS_O_WORKDIR ]]; then
        echo "run"
    else
        echo "submit"
    fi
}


