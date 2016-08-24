#!/usr/bin/zsh
# -*- Shell-script -*-
# submit turbomole job to magnolia

if [[ -z $1 ]]; then
    echo "Usage: $0 (dirname) [(dirname2)...]" 1>&2
    exit 1
fi

rundirs=($@)


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
        ssh $thost $ZSH -i -l <<< "set -x; $cmd"
        ecode=$?
        if (( ecode != 255 )); then
            echo "DEBUG: perform returned $ecode" 1>&2
            return $ecode
        fi
        sleep 60
    done
}

wait_completion() {
    trapoff
    while true; do
        perform "for d in $rundirs; do [[ -e $jobid/GEO_OPT_CONVERGED ]] || exit 1; done"
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

trapon
set -x

wait_completion

