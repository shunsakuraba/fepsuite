# JAEA SGI-8600
# Combination of PBS Pro + Intel MPI (Hydra)

job_prelude() {
    cd $PBS_O_WORKDIR
    source /etc/profile.d/modules.sh
    module purge
    module load intel/2021.4.0 mpt/cur gnu/7.4.0 cuda/11.4 
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
            HW_CPN=40 # 2 sockets, 20 cores/socket
            HW_GPN=0
            CPN=40
            GPN=0
            GPP=0
            ;;
        GPU)
            HW_CPN=48 # 2 sockets, 24 cores/socket
            HW_GPN=4
            CPN=48
            GPN=16
            if (( MULTI == 1 )); then
                # try to submit into shared queue
                HW_CPN=24
                HW_GPN=1
                CPN=24
                GPN=1
            fi
            if [[ -z $CPU_ONLY_STAGE ]]; then
                GPP=1
            fi
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
        if (( PROCS == 1 )) && (( GPP == 1 )); then
            # use sg8 but shrinked size so that it will be run on shared node
            queue=sg8
            JOB_PROCS=1
            JOB_PPN=1
        elif (( JOB_NODES <= 8 )); then
            queue=sg8
        elif (( JOB_NODES <= 32 )); then
            queue=mg32
        elif (( JOB_NODES <= 72 )); then
            queue=lg72
        else
            echo "Too large number of nodes (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
            exit 1
        fi
        if [[ -n $NJOB ]]; then
            # overwrite queue name
            queue=ng72
        fi
    else
        if [[ -n $NJOB ]]; then
            queue=nc192
        elif (( JOB_NODES <= 16 )); then
            queue=sc16
        elif (( JOB_NODES <= 64 )); then
            queue=mc32
        elif (( JOB_NODES <= 192 )); then
            queue=lc192
        else
            echo "Too large number of nodes (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
            exit 1
        fi
    fi

    # check whether project name is set
    if [[ -z $PROJECT ]] && [[ -z $IGNORE_NO_PROJECT ]]; then
        echo "You typically want to set \$PROJECT so that your project points, not personal points, are used" 1>&2
        exit 1
    else
        extra+=(-P $PROJECT)
    fi

    # GPU extra code, PBS requests you to specify "per node" resources, here GPUs per node
    local gpureq
    if (( GPP > 0 )); then
        gpureq=":ngpus=$JOB_GPN"
    else
        gpureq=""
    fi

    # In JAEA qsub ncpus being HW_CPN means exclusive job
    local ncpus_pbs
    (( ncpus_pbs = JOB_PPN * TPP ))
    # JAEA considers half node (vnode) as a unit, if job is below the threshold run on shared
    if (( JOB_NODES > 1 )) || (( ncpus_pbs > HW_CPN / 2 ));  then
        # if using GPU, GPU <= 2 fits vnode
        if (( GPP == 0 )) || (( JOB_GPN > HW_GPN / 2 )) ; then
            # in most cases it ends up here, using exclusive resource
            (( ncpus_pbs = HW_CPN ))
        fi
    fi

    # Run actual code and get jobid
    local -a cmd
    # select: num nodes ("chunks")
    # ncpus: total cpu cores per 1 node
    # mpiprocs: mpi procs per node
    # ompthreads: threads per mpi proc
    cmd=(qsub -q $queue -v ${(j:,:)key_vars} -N $JOB_NAME $dependency -l select=$JOB_NODES:ncpus=$ncpus_pbs:mpiprocs=$JOB_PPN:ompthreads=$TPP$gpureq -l walltime=$timelimit $extra $BASEFILE)
    echo $cmd
    # PBS Pro only emits job id like "12345.s86pbs01", that's great
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

    # omplace use OMP_NUM_THREADS information and $PPN can be ignored
    mpirun -n $N omplace $@
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


