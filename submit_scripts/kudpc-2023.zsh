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
    # if module is not set, probably /etc/profile is not read. This is because bashrc and zshrc behave differently and /etc/profile is not read on non-interactive zsh.
    if ! whence module > /dev/null; then
        source /etc/profile

        case $SUBSYSTEM in
            G)
                module load slurm/2022 SysG/2022 PrgEnvNvidia/2022
                module list
                ;;
            *)
                echo "SUBSYSTEM $SUBSYSTEM not implemented yet"
                exit 1
                ;;
        esac
    fi

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
            # Although Gardenia has 2 sockets/node & 32 cores / socket, the resource allocation algorithm uses 1GPU & 16 core as the unit.
            # Thus we consider "1 GPU" as 1 node.
            HW_CPN=16
            HW_GPN=1
            CPN=16
            GPN=4 # oversubscribe 4 GPUs
            if [[ -z $CPU_ONLY_STAGE ]]; then
                GPP=1
            fi
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
    local -a rsc
    if [[ $SUBSYSTEM = G ]] && [[ $JOB_GPU > 1 ]]; then
        # submit based on GPU number
        # Do not use -g because of t= issue
        rsc=(--rsc p=$((JOB_GPU * HW_CPN / TPP)):t=${TPP}:c=${TPP}:m=${mem}M)
    else
        rsc=(--rsc p=${JOB_PROCS}:t=${TPP}:c=${TPP}:m=${mem}M)
    fi
    local cmd=(sbatch -p $queue --export=${(j:,:)key_vars} -J $JOB_NAME $depstr $rsc -t $timelimit $BASEFILE)
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
    local ntpn=$PPN
    # KUDPC Gardenia does not export envs and thus LD_LIBRARY_PATH is not shared!
    if [[ $SUBSYSTEM = "G" ]]; then
        # PPN is hacked to be 1/4 node slice
        (( ntpn = $PPN * 4 ))
    fi
    srun --export ALL -n $N --ntasks-per-node $ntpn $@
}

job_singlerun() {
    # Even the sequential run needs "srun"
    # same above
    local -a gpus
    if [[ $SUBSYSTEM = "G" ]]; then
        gpus=(--gpus 1)
    fi
    srun $gpus --export ALL -n 1 -N 1 $@
}

job_get_mode() {
    if [[ -n $SLURM_JOB_NAME ]]; then
        echo "run"
    else
        echo "submit"
    fi
}

