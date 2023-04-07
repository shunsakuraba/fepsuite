#!/usr/bin/zsh

BASEFILE=$1
if [[ -z $BASEFILE ]]; then
    echo "This script is not directly callable. Read the manual for the proper use of the program." 2>&1
    exit 1
fi
shift

typeset -a ARGS
ARGS=($@)

if [[ -z $JOBTYPE ]] || [[ -z $FEPSUITE_ROOT ]] || [[ -z $JOBSYSTEM ]]; then
    echo "This script is not directly callable. Read the manual for the proper use of the program." 2>&1
    exit 1
fi

# Check GROMACS
if [[ -z $GROMACS_DIR ]]; then
    echo "You need to specify the location of Gromacs to GROMACS_DIR" 2>&1
    exit 1
fi
if [[ ! -d $GROMACS_DIR ]]; then
    echo "You specified GROMACS_DIR=\"$GROMACS_DIR\", but the directory does not exist or accessible" 2>&1
    exit 1
fi

FEPSUITE_ROOT=${FEPSUITE_ROOT%%/} # remove trailing "/" if exists
if [[ ! -d $FEPSUITE_ROOT/submit_scripts ]]; then
    echo "There is no \"submit_scripts\" directory under FEPSUITE_ROOT" 2>&1
    echo "(Perhaps you set the wrong directory to FEPSUITE_ROOT)" 2>&1
    exit 1
fi

JOB_SCRIPT=$FEPSUITE_ROOT/submit_scripts/$JOBSYSTEM.zsh
if [[ ! -e $JOB_SCRIPT ]]; then
    echo "You specified JOBSYSTEM as \"$JOBSYSTEM\" but no such job system exists under $FEPSUITE_ROOT/submit_scripts" 2>&1
    exit 1
fi

PIPELINE_SCRIPT=${PIPELINE_SCRIPT:-$FEPSUITE_ROOT/$JOBTYPE/pipeline.zsh}
if [[ ! -e $PIPELINE_SCRIPT ]]; then
    echo "You specified JOBTYPE as \"$JOBTYPE\" but no such work type exists under $FEPSUITE_ROOT" 2>&1
    exit 1
fi
# set ABFE_ROOT, FEPREST_ROOT etc.
typeset ${JOBTYPE:u}_ROOT=${PIPELINE_SCRIPT:h}

controller_get_jobid () {
    grep "^$1\t" $ID/jobid.txt | tail -n 1 | cut -f2
}

# Generic submission process. Most of them are handled inside $JOBSYSTEM.zsh,
# but tedious calculations and resource handling should fit here
controller_submit() {
    typeset -a STAGES
    STAGES=($@)

    job_init_queue_stat

    if [[ all = $STAGES[1] ]]; then
        STAGES=($($PIPELINE_SCRIPT query all))
    fi

    # ceil(a/b)
    zmath_ceildiv() {
        typeset -i a
        typeset -i b
        a=$1
        b=$2
        typeset -i ret
        (( ret = (a + b - 1) / b ))
        return ret
    }
    functions -M ceildiv 2 2 zmath_ceildiv

    for STEPNO in $STAGES; do
        # set default variables, these variables are overwritten in next "eval" statement
        MULTI=1
        DEPENDS=()

        unset PPM
        # get necessary resource information
        eval $($PIPELINE_SCRIPT query $STEPNO)
        if [[ -z $PPM ]]; then 
            echo "Error: pipeline script returned \$PPM as empty"
            exit 1
        fi

        # Calculate necessary information to submit
        # Note job_set_preferred_resource "default" resource may depend on requested variables
        unset JOB_PPN || true
        unset JOB_NODES || true
        job_set_preferred_resource

        # PROCS: Total MPI ranks (processes)
        # PPM: MPI ranks (processes) per multidir runs (for non-replica run PPM = PROCS)
        # MULTI: number of multidir runs, i.e., number of replicas in replex calculations (otherwise 1)
        (( PROCS = MULTI * PPM ))

        if (( PROCS == 0 )); then
            echo "Aborting because \$PROCS = 0 (\$MULTI = $MULTI, \$PPM = $PPM); this is probably because the parameters in either run.zsh or para_conf.zsh is not properly set" 2>&1
            exit 1
        fi

        # thread parallelism cannot cross the node boundary, so the intention of this division is floor(CPN/TPP)
        (( CPU_PPN = CPN / TPP ))
        if (( HW_GPN > 0 )) && (( GPP > 0 )) && [[ ! $CPU_ONLY_STAGE = yes ]]; then
            (( GPU_PPN = GPN / GPP ))
        else
            (( GPU_PPN = CPU_PPN ))
        fi
        # PPN is smaller of the two
        if (( CPU_PPN < GPU_PPN )); then
            (( PPN = CPU_PPN ))
        else
            (( PPN = GPU_PPN ))
        fi

        # Required number of nodes
        (( NODES = ceildiv(PROCS, PPN) ))

        # Number of resources the program requests to job system is different from the resource the program is actually using.
        # This is because oversubscribing / undersubscribing GPU resource is not "typical" operation in many job queueing systems.
        # While oversubscribing is unwelcomed on CPU, running 2-4 mdruns per 1 GPU is beneficial.

        # Total REQUESTED ranks (not CPN * NODES); REAL ranks are determined in the "job_mpirun" part.
        JOB_PPN=${JOB_PPN:-$PPN}
        JOB_NODES=${JOB_NODES:-$NODES}
        (( JOB_PROCS = JOB_NODES * JOB_PPN ))
        # Total CPU cores to be used
        (( JOB_CPU = JOB_PROCS * TPP ))

        # Total GPU used in the job; JOB_GPN is per node and JOB_GPU is total GPUs
        (( JOB_GPN = GPP * CPN / TPP ))
        if (( JOB_GPN > HW_GPN )); then
            # Pass HW_GPN to queue system (declaring full node use), but inside the pipeline we still use requested number of GPUs
            JOB_GPN=$HW_GPN
        fi
        (( JOB_GPU = JOB_GPN * JOB_NODES ))

        # Set default job name
        JOB_NAME="$ID.$STEPNO"
        # In almost all job systems job name starting from a number is not allowed
        case $JOB_NAME in
            [0-9]*)
                JOB_NAME="R$JOB_NAME"
                ;;
        esac

        if [[ -n $DEBUG_JOB_SUBMIT ]]; then
            typeset > $ID/debug.submit.$STEPNO.txt
            set -x
        fi

        echo "Submitting job with following settings:"
        echo " PROCS=$PROCS  # Number of requested processes (ranks)"
        echo " JOB_PPN=$JOB_PPN  # Process / node"
        echo " JOB_NODES=$JOB_NODES  # Number of requested nodes (rounding up)"
        echo " JOB_PROCS=$JOB_PROCS  # Number of requested processes in the job system (rounding up)"
        echo " TPP=$TPP  # Number of CPU threads / process"
        echo " JOB_CPU=$JOB_CPU  # Total CPU cores to be used"
        echo " GPU=$((GPP*CPN/TPP*JOB_NODES))  # Total logical GPUs used in the job"
        echo " JOB_GPU=$JOB_GPU  # Total physical GPUs used in the job"

        # run job_submit() and get the job ID
        unset JOBID
        job_submit
        if (( $? != 0 )); then
            echo "Job submission failed. Typical reason is your job in the queue exceeded the maximum allowed number." 2>&1
            exit 1
        fi
        if [[ -z $JOBID ]]; then
            echo "JOBID is not set; probably it means a failure to submit a job" 2>&1
            exit 1
        fi
        # keep record of submitted jobs, may use in some batch system
        echo "$STEPNO\t$JOBID" >> $ID/jobid.txt
    done
}

#################
# Main process starts from here

# Read Job-system specific settings and functions
source $JOB_SCRIPT

case $(job_get_mode) in
    run)
        # Execute job prelude (load modules, env vars etc)
        job_prelude $BASEFILE

        if [[ -z $ID ]] || [[ -z $STEPNO ]]; then
            echo "\$ID(=\"$ID\") or \$STEPNO(=\"$STEPNO\") is not set (likely an error in submit_scripts/$JOBSYSTEM.zsh)" 2>&1
            exit 1
        fi

        # Read para_conf after the prelude because $PWD may not be the same directory as that when submitted
        source para_conf.zsh
        [[ -e $ID/para_conf.zsh ]] && source $ID/para_conf.zsh

        # Note GROMACS_DIR is set inside GMXRC. If GMXRC is already sourced, this should overwrite GROMACS_DIR with the same variable again.
        source $GROMACS_DIR/bin/GMXRC.zsh

        # Currently all pipeline uses replica exchange (single precision), so we assume gmx_mpi exists.
        GMX=${GMX:-$(which gmx_mpi)}
        GMX_MPI=${GMX_MPI:-$(which gmx_mpi)}

        # typically you need nothing to do in this part
        job_prelude_after_gmx $BASEFILE

        if [[ -n $DEBUG_JOB_SUBMIT ]]; then
            typeset > $ID/debug.run.$STEPNO.txt
            env > $ID/debug.env.$STEPNO.txt
        fi

        # Run the real pipeline
        setopt ERR_EXIT
        source $PIPELINE_SCRIPT run $STEPNO
    ;;
    submit)
        ID=${ARGS[1]}
        if [[ -z $ID ]] || (( ${#ARGS} <= 1 )); then
            echo "Usage: $BASEFILE (subdirectory name) (stages)"
            echo "Example: $BASEFILE A31V all"
            echo "Example: $BASEFILE MYMOL 3 4 5 6"
            exit 1 
        fi
        shift ARGS
        source para_conf.zsh
        [[ -e $ID/para_conf.zsh ]] && source $ID/para_conf.zsh
        controller_submit $ARGS
    ;;
    *)
        echo "Unknown mode returned from job_get_mode(), likely an error in $JOB_SCRIPT"
        exit 1
    ;;
esac


