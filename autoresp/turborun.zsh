#!/usr/bin/zsh
# -*- Shell-script -*-
# submit turbomole job to RCCS

if [[ -z $2 ]]; then
    echo "Usage: $0 (dirname) FOO=BAR BAR=BAZ (turbomole script)" 1>&2
    exit 1
fi

rundir=$1
shift

ENVS=()

while true; do
    EQS=$1
    case $EQS in
	*=*)
	    ENVS+=$EQS
	    shift
	;;
	*)
	    input=$EQS
	    break
	;;
    esac
done

remotedirbase=save/turbomolequeue
queue=PG
tsubmit=/home/users/rq0/$remotedirbase/run.csh
dir=$remotedirbase/$rundir

script=${input:t}

thost=automation-ccfep.ims.ac.jp
ZSH=/home/users/rq0/cx/lcl/bin/zsh
mo='~/opt/pg/mo'

perform() {
    cmd=$1
    while true; do
        ssh $thost $ZSH -i -l <<< "set -x; mkdir -p $dir; cd $dir; $cmd"
        ecode=$?
        if (( ecode != 255 )); then
            echo "DEBUG: perform returned $ecode" 1>&2
            return $ecode
        fi
        sleep 60
    done
}

wait_completion() {
    jobid=$1
    trapoff
    while true; do
        perform "[ ! -e $jobid ]"
        ecode=$?
        echo "DEBUG: wait_completion returned $ecode" 1>&2
        if (( ecode == 1 )); then 
	    trapon
            return 0
        fi
        sleep 300
    done
    trapon
}

trapon() {
    trap "echo Error \$? returned while sending command 1>&2 ; exit 1" ZERR    
}

trapoff() {
    trap - ZERR
}



ENVF=$(tempfile -p envir)

for e in $ENVS; do
    K=${e%%=*}
    V=${e#*=}
    if [[ $K = CONSTRAINTS ]]; then
	CONSTRAINTFILE=$V
	continue
    fi
    echo "$K=\"${V}\""
done > $ENVF

trapon
set -x

ssh $thost 'mv -f $dir '${dir}'.bak; mkdir -p '$dir
scp $input $thost:$dir/coord.in
scp $ENVF $thost:$dir/env.txt
if [[ ! -z $CONSTRAINTFILE ]]; then
    scp $CONSTRAINTFILE $thost:$dir/constraints
fi

perform "{ set -x; cd $dir; rm -f GEO_OPT_CONVERGED; $tsubmit } |& tee run.log"
if [[ -z $NOWAIT ]]; then
    wait_completion GEO_OPT_CONVERGED
fi

rm $ENVF

