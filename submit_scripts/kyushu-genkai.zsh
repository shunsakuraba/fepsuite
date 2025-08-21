# Kyushu University Genkai Node B (GPU)
# Combination of PJM (Fujitsu job manager) + GCC 12 + Intel MPI
#
# Installation notes for users:
# (1) Prepare python side.
# % pip3.11 install numpy cython==0.29.37 six --user
# % pip3.11 install mdtraj pymbar==3.0.3 pyedr --user
# (2) Load necessary modules. Apply patches if necessary.
# % module load cuda/12.2.2 gcc-toolset/12
# % module load impi/2021.12
# % ((( download and untargz GROMACS. For FEP/REST, also apply the patch )))
# % cd PATH_TO_GROMACS
# % mkdir build; cd build
# % cmake -DCMAKE_INSTALL_PREFIX=WHERE_TO_INSTALL_GROMACS -DGMX_GPU=CUDA -DGMX_CUDA_TARGET_COMPUTE="90" -DGMX_BUILD_OWN_FFTW -DGMX_MPI=on  ..
# ^- H100 has compute capabaility of 9.0
# % make -j4 && make install
# (3) Update run.zsh to correct path and run.
#

job_prelude() {
    cd $PJM_O_WORKDIR
    source /etc/profile.d/modules.sh
    module load cuda/12.2.2 gcc-toolset/12
    module load impi/2021.12
    PYTHON3=python3.11
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
    HW_CPN=120 # 2 sockets, 60 cores/socket, HT disabled
    HW_GPN=4
    CPN=120
    GPN=32

    if (( MULTI == 1 )); then
        # try to submit into share
        CPN=30 # 1/4
        HW_CPN=30
        GPN=1
        HW_GPN=1
    fi

    if (( PPM > 1 )); then
        echo "kyushu-genkai.zsh is designed to run with GPU, and for GPU using more than one process per -multidir is not welcomed." 1>&2
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
    queue=b-batch
    if (( GPP < 4 )) && (( PROCS <= 30 * GPP )) ; then
        # shared resource group
        extra=(-L gpu=$GPP)
        JOB_PROCS=$PROCS
        JOB_PPN=$PROCS
    else
        # request node instead of GPUs
        extra=(-L node=$JOB_NODES) 
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

    mpiexec -ppn $PPN -n $N $@
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


