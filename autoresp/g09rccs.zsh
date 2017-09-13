#!/usr/bin/zsh

# submit g09 job to magnolia

if [[ -z $2 ]]; then
    echo "Usage: $0 (dirname) (gaussian script)" 1>&2
    exit 1
fi

rundir=$1
input=$2

queue=PF
remotedirbase=save/g09queue
script=${input:t}
dir=$remotedirbase/$rundir/${script:r}

thost=ccfep.ims.ac.jp

GMAJOR=${GMAJOR:-09}
subg09=${subg09:-g${GMAJOR}sub}
if [[ -n $GMINOR ]]; then
    REV=(-rev g${GMAJOR}${GMINOR})
else
    REV=()
fi

perform() {
    cmd=$1
    while true; do
        ssh $thost /home/users/rq0/cx/lcl/bin/zsh -i -l <<< "set -x; mkdir -p $dir; cd $dir; $cmd"
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
        perform "while true; do sleep 300; qjobs -c | grep -q \" $jobid \" && continue; exit 1; done"
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

ssh $thost 'mkdir -p '$dir'; cd '$dir'; cat > '$script < $input

if [[ -n $EXTRA ]]; then
    FILES=($(tr ',' ' ' <<< $EXTRA))
    scp $FILES $thost:$dir
fi

if [[ -n $LONG ]]; then
    subopt=($subopt -walltime 72:00:00)
else
    subopt=($subopt -walltime 12:00:00)
fi

JOBID=$(perform "{ set -x; $CMD; $subg09 -q $queue $REV -np 1 $script $subopt |& tee jobid.txt } |& tee run.log | tail -1")
if [[ -n $NOWAIT ]]; then
    exit 0
fi
JOBID=$(tr '<>' '  ' <<< $JOBID | cut -f3 -d' ')

case $JOBID in
*.*)
JOBID=${JOBID%%.*}
;;
esac

if [[ -z $NOWAIT ]]; then
    wait_completion $JOBID
fi


