# Run without batch job system

job_prelude() {
    PYTHON3=python3
}

job_prelude_after_gmx() {
    # do nothing
}

job_init_queue_stat() {
    # do nothing, without batch job system
}

job_set_preferred_resource() {
    # Set as "unlimited"
    CPN=9999
    HW_CPN=9999
    GPN=9999
    HW_GPN=9999
}

job_submit() {
    # Export variables
    OMP_NUM_THREADS=$TPP
    local exports=(PROCS OMP_NUM_THREADS PPN GPN ID STEPNO)
    for e in $exports; do
        export $e
    done

    # Run actual code
    local cmd=(RUN=yes $BASEFILE)
    echo $cmd
    local jobidinfo=$($cmd)
    if [[ $? != 0 ]]; then
        echo "Submission failed" 1>&2
        exit 1
    fi
    # jobid does not exist
    JOBID="(not_used)"
}

job_mpirun() {
    local N=$1
    shift

    # OpenMPI version
    # For other MPIs, you can skip -npernode $PPN
    mpirun -npernode $PPN -n $N $@
}

job_singlerun() {
    $@
}

job_get_mode() {
    if [[ -n $RUN ]]; then
        echo "run"
    else
        echo "submit"
    fi
}


