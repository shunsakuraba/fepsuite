# Kyoto University supercomputer (2023 spring-)

if [[ -z $GROUP ]]; then
    GROUP=$(groups | tr ' ' '\n' | grep '^gr')
fi
if [[ -z $SUBSYSTEM ]]; then
    echo "Specify SUBSYSTEM (currently allowed choices: CL/D/G)" 2>&1
    exit 1
fi

job_prelude() {
    module purge
    module load Sys$SUBSYSTEM/2022 PrgEnvGCC/2022
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
            echo "Hardware in cloud system is unknown" 2>&1
            exit 1
        ;;
    esac
}

job_submit() {
    # Export variables
    local exports=(PROCS PPN ID STEPNO)
    local key_vars=()
    for e in $exports; do
        key_vars+="$e=${$e}"
    done

    # Default time
    local timelimit=8:00:00

    # Set dependency
    local depstr=()
    if (( ${#DEPENDS} > 0 )); then
        local deps=()
        for d in $DEPENDS; do
            deps+=controller_get_jobid $d
        done
        depstr=(--dependency=afterok:${(j|:|)deps})
    fi

    # set queue
    # local queue=${GROUP}${(SUBSYSTEM:l}
    # During testing period
    local queue
    case $SUBSYSTEM in
        D)
            queue=s32d
        ;;
        G)
            queue=eg
        ;;
        CL)
            queue=so
        ;;
    esac

    # Run actual code and get jobid
    local mem
    (( mem = 4750 * CPP ))
    local cmd=(sbatch -p $queue --export=${(j:,:)key_vars} -J $JOB_NAME $depstr --rsc p=$JOB_RANK:t=$CPP:c=$CPP:m=${mem}M -t $timelimit $BASEFILE)
    echo $cmd
    local jobidinfo=$($cmd)
    if [[ $? != 0 ]]; then
        echo "Submission failed" 2>&1
        exit 1
    fi
    # Output is like: "sbatch: Submitted batch job 34987".
    local jobid=(${(@s: :)jobidinfo})
    echo ${jobid[5]} # output is sent back to controller
}

job_mpirun() {
    local N=$1
    shift

    srun -n $N --ntasks-per-node $PPN $@
}

job_singlerun() {
    # Even the sequential run needs "srun"
    srun -n 1 $@
}

job_get_mode() {
    if [[ -n $SBATCH_JOB_NAME ]]; then
        echo "run"
    else
        echo "submit"
    fi
}

