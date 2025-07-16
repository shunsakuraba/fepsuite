# The University of Tokyo Wisteria/BDEC-01
# Combination of PJM (Fujitsu job manager) + Intel MPI
#
# Installation notes for users:
# (0) (recommended) Check the Supercomputer User Guide and search "chhome"
# to change your home directory to /work/(groupid)/(uid), otherwise it's painful to run.
# (1) Prepare python side.
# % module load python/3.10.13          # Wisteria's default python 3.6.8 is way too old and is now deprecated.
# % pip3 install numpy cython==0.29.37 six
# % pip3 install mdtraj pymbar==3.0.3 pyedr
# (2) Load necessary modules. Apply patches if necessary.
# % module load gcc/8.3.1  cuda/12.2 impi/2021.8.0
# % cd PATH_TO_GROMACS             # For FEP/REST, do not forget patching
# % mkdir build; cd build
# % cmake -DCMAKE_INSTALL_PREFIX=WHERE_TO_INSTALL_GROMACS -DGMX_GPU=CUDA -DGMX_CUDA_TARGET_COMPUTE="80" ..            # A100 has compute capabaility of 8.0
# % make -j4 && make install
# (3) Update run.zsh to correct path and run.
#

job_prelude() {
    cd $PJM_O_WORKDIR
    source /etc/profile.d/modules.sh
    module load gcc/8.3.1  cuda/12.2 impi/2021.8.0
    module load python/3.10.13
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
    HW_CPN=72 # 2 sockets, 36 cores/socket
    HW_GPN=8
    CPN=72
    GPN=32

    if (( MULTI == 1 )); then
        # try to submit into share
        CPN=9 # 1/8
        HW_CPN=9
        GPN=1
        HW_GPN=1
    fi

    if (( PPM > 1 )); then
        echo "tokyo-wisteria.zsh is designed to run with GPU, and for GPU using more than one process per -multidir is not welcomed." 1>&2
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
    if [[ -n $TIMELIMIT ]]; then
        timelimit=$TIMELIMIT
    fi

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
    if (( GPP <= 4 )) && (( PROCS <= 9 * GPP )) ; then
        queue=share
        JOB_PROCS=$PROCS
        JOB_PPN=$PROCS  # fit into a single node...
        extra=(-L gpu=$GPP)
    elif (( JOB_NODES <= 8 )); then
        queue=regular-a
        extra=(-L node=$JOB_NODES)
    else
        echo "Unimplemented node number (JOB_PROCS=$JOB_PROCS JOB_NODES=$JOB_NODES)" 1>&2
        exit 1
    fi

    # Set gname. Use default group if not set
    if [[ -z $GROUP ]]; then
        GROUP=$(id -gn)
    fi

    # Run actual code and get jobid
    local -a cmd
    cmd=(pjsub -L rscgrp=$queue -x ${(j:,:)key_vars} -N $JOB_NAME $dependency -g $GROUP --mpi proc=$JOB_PROCS -L elapse=$timelimit $extra $BASEFILE)
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

    mpiexec.hydra -ppn $PPN -n $N $@
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


