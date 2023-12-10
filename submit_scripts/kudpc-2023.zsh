# Kyoto University supercomputer (2023 spring-)
# Combination of slurm (cutomized) + Intel MPI

if [[ -z $GROUP ]]; then
    GROUP=$(groups | tr ' ' '\n' | grep '^gr')
fi
if [[ -z $SUBSYSTEM ]]; then
    echo "Specify SUBSYSTEM (currently allowed choices: CL/D/G)" 1>&2
    exit 1
fi

job_prelude() {
    # if module is not set, probably /etc/profile is not read (happens because zsh is non-default shell on KUDPC)
    whence module > /dev/null || source /etc/profile

    # FIXME: this should depend on $SUBSYSTEM
    module load PrgEnvNvidia/2022

    # In KUDPC, default pip3 = pip3.6, and thus default python should be 3.6
    # In the past it was 3.8 but they seem fixed it
    PYTHON3=python3
}

job_prelude_after_gmx() {
    # do nothing
}

job_init_queue_stat() {
    # do nothing, slurm is a good boy and we don't need queue state information
}

job_set_preferred_resource() {
    case $SUBSYSTEM in
        D)
            HW_CPN=40 # 2 sockets, 20 cores/socket
            HW_GPN=0
            CPN=40
            GPN=0
        ;;
        G)
            HW_CPN=64 # 2 sockets, 32 cores/socket
            HW_GPN=4
            CPN=64
            GPN=16
        ;;
        CL)
            echo "Hardware in cloud system is unknown" 1>&2
            exit 1
        ;;
        *)
            echo "Unknown \$SUBSYSTEM ($SUBSYSTEM). Must be one of D/G"
            exit 1
        ;;
    esac
}

job_submit() {
    # Export variables
    local -a exports
    exports=(PROCS PPN GPP ID STEPNO)
    local -a key_vars
    key_vars=()
    for e in $exports; do
        key_vars+="$e=${(P)e}"
    done

    # Default time
    local timelimit="8:00:00"

    # Set dependency
    local -a depstr
    depstr=()
    if (( ${#DEPENDS} > 0 )); then
        local deps=()
        for d in $DEPENDS; do
            deps+=$(controller_get_jobid $d)
        done
        depstr=(--dependency=afterok:${(j|:|)deps})
    fi

    # set queue
    local queue="${GROUP}${SUBSYSTEM:l}"

    # Run actual code and get jobid
    local mem
    (( mem = 4750 * TPP ))
    local cmd=(sbatch -p $queue --export=${(j:,:)key_vars} -J $JOB_NAME $depstr --rsc p=${JOB_PPN}:t=${TPP}:c=${TPP}:m=${mem}M -t $timelimit $BASEFILE)
    echo $cmd
    local jobidinfo=$($cmd)
    if [[ $? != 0 ]]; then
        echo "Submission failed" 1>&2
        exit 1
    fi
    # Output is like: "Submitted batch job 34987".
    local jobid=(${(@s: :)jobidinfo})
    JOBID=${jobid[4]} # output is sent back to controller
    return 0
}

job_mpirun() {
    local N=$1
    shift
    # KUDPC Gardenia does not export envs and thus LD_LIBRARY_PATH is not shared!
    srun --export ALL -n $N --ntasks-per-node $PPN $@
}

job_singlerun() {
    # Even the sequential run needs "srun"
    # same above
    srun --export ALL -n 1 $@
}

job_get_mode() {
    if [[ -n $SLURM_JOB_NAME ]]; then
        echo "run"
    else
        echo "submit"
    fi
}

