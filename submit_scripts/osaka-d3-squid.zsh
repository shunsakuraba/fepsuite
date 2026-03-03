# Osaka Univ. D3 center SQUID
# Combination of NQSV (PBS-like) + Open MPI
#
# How to compile:
# 
# module purge
# module load BasePy/2025
# module load BaseGCC/2025 cuda/12.6
# (for FEP/REST only, apply patch here)
# mkdir build; cd build
# cd build
# cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/gromacs -DGMX_GPU=cuda -DGMX_MPI=on -DGMX_BUILD_OWN_FFTW=on -DGMX_CUDA_TARGET_COMPUTE="80" -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" ..
# make -j8 && make install
#
# You will also need to install python libraries AFTER loading BasePy/2025:
# pip3 install pyedr cython
# pip3 install mdtraj pymbar==3.0.3

job_prelude() {
    cd $PBS_O_WORKDIR
    source /etc/profile.d/modules.sh
    module purge
    module load BasePy/2025
    module load BaseGCC/2025 cuda/12.6
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
    # Each SQUID node has 8 A100 GPUs. Clamming 4+ jobs on each GPU is fine.
    HW_CPN=76
    HW_GPN=8
    CPN=76
    GPN=32

    if (( MULTI == 1 )); then
        # try to submit into shared (SQUID-S)
        CPN=38
        HW_CPN=38
        GPN=16
        HW_GPN=4
    fi

    if (( PPM > 1 )) ; then
        echo "$0 is designed to run with GPU, and for GPU, using more than one process per -multidir is not optimal. Set PARA=1 for GPU." 1>&2
        exit 1
    fi

    if [[ -z $CPU_ONLY_STAGE ]]; then
        GPP=1
    fi
}

job_submit() {
    # Export variables
    local -a exports
    exports=(PROCS TPP PPN GPP ID STEPNO)
    local -a key_vars
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
    #extra+=(-M nobody@example.com)
    # Assumes standard or HPCI use with GPU nodes.
    # HPCI user can't use SQUID-S.
    if (( GPP <= 4 )) && (( PROCS <= 38 )) && (( JOB_NODES == 1 )) && [[ $(id -g -n) != hpci ]]; then
        queue=SQUID-S
    elif (( JOB_NODES <= 32 )); then
        queue=SQUID
    else
        echo "Too large number of nodes (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
        exit 1
    fi

    # Set group name
    if [[ -z $GROUP ]]; then
        GROUP=$(id -G -n)
        GROUP=${GROUP##* }
    fi
    extra+=(--group=$GROUP)

    # use OpenMPI
    extra+=(-T openmpi)
    extra+=(-v NQSV_MPI_MODULE=BaseGCC)

    # Run actual code and get jobid
    local -a cmd

    # All preps done, now submitting the job
    # Required: elapsetim_req, cpunum_job,gpunum_job,memsz_job. 80GB should be enough and it should fit SQUID-S
    cmd=(qsub -q $queue -v ${(j:,:)key_vars} -N $JOB_NAME $dependency -b $JOB_NODES -l cpunum_job=$CPN -l gpunum_job=$JOB_GPN -l memsz_job=80GB -l elapstim_req=$timelimit $extra $BASEFILE --cancel-after)
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
    mpirun -x OMP_NUM_THREADS $NQSV_MPIOPTS_A -np $N -map-by ppr:$PPN:node --bind-to core $@
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


