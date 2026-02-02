# JCAHPC Miyabi
# Combination of PBS Pro + OpenMPI
#
# Installation notes for users:
# (1) Prepare python side.
# % module unload nvidia nv-hpcx     # need to unload because CC is overwritten to be nvc
# % module load gcc/12.4.0
# % module load cuda/12.6
# % module load ompi-cuda/4.1.6-12.6 
# % pip3 install numpy cython==0.29.37 six
# % pip3 install mdtraj pymbar==3.0.3 pyedr
# % python3 -c "import mdtraj"       # No error should be emit here
# (2) Apply patches if necessary.
# % cd PATH_TO_GROMACS             # For FEP/REST, do not forget patching here
# % mkdir build; cd build
# % cmake -DCMAKE_INSTALL_PREFIX=WHERE_TO_INSTALL_GROMACS -DGMX_GPU=CUDA -DGMX_CUDA_TARGET_COMPUTE="90" -DGMX_MPI=ON ..            # GH200 has compute capabaility of 9.0
# % make -j4 && make install
# (3) Update run.zsh, add group id, correct path and run.
# GROUP=hp000000 <- your group ID. Note it should end with number, not "g".
#

job_prelude() {
    cd $PBS_O_WORKDIR
    source /etc/profile.d/modules.sh
    module purge
    module load gcc/12.4.0
    module load cuda/12.6
    module load ompi-cuda/4.1.6-12.6 
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
    LIVE_ID_LIST=($(qstat | cut -d' ' -f1 | tail -n +4 | sed 's/$/.opbs/'))
}

job_set_preferred_resource() {
    HW_CPN=72
    HW_GPN=1
    CPN=72
    GPN=8    # overbook up to 8 procs
    GPP=1
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
        echo "DEBUG \"$jid\" \"$LIVE_ID_LIST\"" 1>&2
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
    queue=regular-g

    if [[ -z $GROUP ]]; then
        echo "GROUP is not set. Run \"show_token\" first, then set GROUP in run.zsh according to the appropriate group ID." 1>&2
        exit 1
    fi

    # Run actual code and get jobid
    local -a cmd
    # select: num nodes
    # mpiprocs: mpi procs per node
    # ompthreads: threads per mpi proc
    cmd=(qsub -q $queue -v ${(j:,:)key_vars} -N $JOB_NAME $dependency -l select=$JOB_NODES:mpiprocs=$JOB_PPN:ompthreads=$TPP -l walltime=$timelimit -W group_list=$GROUP $extra $BASEFILE)
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


