# Nagoya University supercomputer flow type II (cx)
# Combination of PJM (Fujitsu job manager) + OpenMPI

job_prelude() {
    cd $PJM_O_WORKDIR
    source /etc/profile.d/modules.sh
    module purge
    module load gcc/10.3.0 cuda/11.4.2
    module load openmpi_cuda/4.1.2
    PYTHON3=python3
    ID=$WID
    export OMP_NUM_THREADS=$TPP
}

job_prelude_after_gmx() {
    # do nothing
}

job_init_queue_stat() {
    # do nothing, Fujitsu pjm is a good boy and we don't need queue state information
}

job_set_preferred_resource() {
    HW_CPN=40 # 2 sockets, 20 cores/socket
    HW_GPN=4
    CPN=40
    GPN=16

    if (( MULTI == 1 )); then
        # try to submit into cx-share
        CPN=10
        HW_CPN=10
        GPN=1
        HW_GPN=1
    fi

    if (( PPM > 1 )); then
        echo "nagoya-cx.zsh is designed to run with GPU, and for GPU using more than one process per -multidir is not welcomed." 1>&2
        exit 1
    fi

    if [[ -z $CPU_ONLY_STAGE ]]; then
        GPP=1
    fi
}

job_submit() {
    # Export variables
    # rename ID to WID because it conflicts with PJ variables
    WID=$ID
    local -a exports
    exports=(PROCS TPP PPN GPP WID STEPNO)
    local -a key_vars
    key_vars=()
    for e in $exports; do
        key_vars+="$e=${(P)e}"
    done

    # Default time
    local timelimit="8:00:00"

    # Set dependency. Instead of (logically correct) direct dependency on all previous jobs, use step execution (= step-by-step execution of stages)
    local -a depparams
    depparams=()
    local jnam
    jnam=$ID
    case $jnam in
        [0-9]*)
            jnam=J$jnam # avoid starting from numbers
            ;;
    esac
    depparams+="jnam=$jnam"
    depparams+="sn=$STEPNO"
    if [[ -n $DEPENDS ]] && (( ${#DEPENDS} > 0 )); then
        depparams+="sd=ec!=0:after"
    fi
    local dependency
    dependency=(--step --sparam ${(j:,:)depparams})

    # set queue
    local extra
    if (( PROCS == 1 )) && (( GPP == 1 )); then
        queue=cx-share
        JOB_PROCS=1
        JOB_PPN=1
        extra=(-L gpu=1)
    elif (( JOB_NODES <= 8 )); then
        queue=cx-small
        extra=(-L node=$JOB_NODES)
    elif (( JOB_NODES <= 16 )); then
        queue=cx-middle
        extra=(-L node=$JOB_NODES)
    elif (( JOB_NODES <= 64 )); then
        queue=cx-large
        extra=(-L node=$JOB_NODES)
    else
        echo "Unimplemented node number (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
        exit 1
    fi

    # Run actual code and get jobid
    local -a cmd
    cmd=(pjsub -L rscgrp=$queue -x ${(j:,:)key_vars} -N $JOB_NAME $dependency --mpi proc=$JOB_PROCS -L elapse=$timelimit $extra $BASEFILE)
    echo $cmd
    local jobidinfo
    jobidinfo=$($cmd)
    if [[ $? != 0 ]]; then
        echo "Submission failed" 1>&2
        exit 1
    fi
    # Output is like: "[INFO] PJM 0000 pjsub Job 932267_1 submitted.".
    JOBID=(${(@s: :)jobidinfo})
    JOBID=${JOBID[6]} # this value is used in the controller
}

job_mpirun() {
    local N=$1
    shift

    # --mca orte_tmpdir_base is required when submitting to cx-shared (ticket 00041609)
    mpirun -machinefile $PJM_O_NODEINF --mca orte_tmpdir_base /tmp/${PJM_JOBID} -npernode $PPN -n $N $@
}

job_singlerun() {
    $@
}

job_get_mode() {
    if [[ -n $PJM_O_WORKDIR ]]; then
        echo "run"
    else
        echo "submit"
    fi
}


