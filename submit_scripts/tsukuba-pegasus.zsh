# Tsukuba Univ. Pegasus
# Combination of NQSV (PBS-like) + Open MPI
#
# How to compile:
# 
# module purge
# module load intelpython/2022.3.1
# module load openmpi/4.1.8/gcc11.4.0-cuda11.8.0
# (for FEP/REST only, apply patch here)
# mkdir build; cd build
# cd build
# cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/gromacs -DGMX_GPU=cuda -DGMX_MPI=on -DGMX_BUILD_OWN_FFTW=on ..

job_prelude() {
    cd $PBS_O_WORKDIR
    source /etc/profile.d/modules.sh
    module purge
    module load intelpython/2022.3.1
    module load openmpi/4.1.8/gcc11.4.0-cuda11.8.0
    PYTHON3=python3
    export OMP_NUM_THREADS=$TPP
}

job_prelude_after_gmx() {
    # do nothing
}

job_init_queue_stat() {
    # Load "live" job id list.
    # We can't use "dead" ID number in "--after (ID)"
    # Live ID list should be global because we want to use this throughout the process
    typeset -ag LIVE_ID_LIST
    LIVE_ID_LIST=($(qstat | cut -d' ' -f1 | tail -n +3))
}

job_set_preferred_resource() {
    # Pegasus's H100 has 7 or 8 GPC. Proper GPN value is debatable, but clamming 2+ jobs on each GPC is at least necessary.
    # Thus I'm setting GPN to be 16...
    HW_CPN=48
    HW_GPN=1
    CPN=48
    GPN=16
    GPP=1
}

job_submit() {
    # Export variables
    local -a exports
    exports=(PROCS TPP PPN GPP ID STEPNO)
    local -a key_vars
    # NQSV_MPI_VER is mandatory. This should match with module name loaded above as "openmpi/xxxxxxx"
    key_vars=(NQSV_MPI_VER=4.1.8/gcc11.4.0-cuda11.8.0)
    # for openmpi this is required
    key_vars+=(NQSV_PBS_NODEFILE_R115=ON)
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
        # --after ID only accepts "live" ID in the jobsystem. After the job completion the ID is invalidated.
        # Note in the pathological case, the ID may turn dead while running this script, but we can't do anything (if anyone knows the way please teach me!)
        # Expansion I: index of last occurence (0 not found)
        # Expansion e: exact match
        if (( $LIVE_ID_LIST[(Ie)$jid] )); then
            deps+=$jid
        fi
    done
    local -a dependency
    if (( ${#deps} > 0 )); then
        dependency=(--after ${(j|,|)deps})
    else
        dependency=()
    fi

    # set queue
    local queue
    local -a extra
    extra=()
    #extra=(-m n) # disables mail notification
    # This is a bug in the pegasus wrapper script, but -m n requests -M
    #extra+=(-M nobody@example.com)
    # --> this mail can be disabled if not specified
    # In the current version this is not necessary
    # Assumes standard or HPCI use. You may be requested to use gpu queue.
    if (( JOB_NODES <= 31 )); then
        queue=gen_S
    elif (( JOB_NODES <= 63 )); then
        queue=gen_M
    elif (( JOB_NODES <= 150 )); then
        queue=gen_L
    else
        echo "Too large number of nodes (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
        exit 1
    fi

    # Set group name
    if [[ -z $GROUP ]]; then
        GROUP=$(id -g -n)
    fi
    extra+=(-A $GROUP)

    # use OpenMPI
    extra+=(-T openmpi)


    # Run actual code and get jobid
    local -a cmd


    # default qsub wrapper has several bugs and will submit without.
    cmd=(/system/tool/nqsv/src/qsub_wrapper2 -q $queue -v ${(j:,:)key_vars} -N $JOB_NAME $dependency -b $JOB_NODES -l elapstim_req=$timelimit $extra $BASEFILE --cancel-after)
    echo $cmd 1>&2
    RET=($($cmd))
    if [[ $? != 0 ]]; then
        echo "Submission failed" 1>&2
        exit 1
    fi
    echo $RET 1>&2
    JOBID=$(cut -f2 -d' ' <<< $RET)
    LIVE_ID_LIST+=$JOBID
}

job_mpirun() {
    local N=$1
    shift

    # --bind-to core needed because no thread pinning without
    # $NQSV_MPIOPTS_A doesn't work on zsh
    NQSV_MPIOPTS_A=(${(@s: :)NQSV_MPIOPTS})
    mpirun -x OMP_NUM_THREADS -n $N $NQSV_MPIOPTS_A -map-by ppr:$PPN:node --bind-to core $@
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


