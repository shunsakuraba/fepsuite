#!/usr/bin/zsh

# submit g09 job to laurel

if [[ -z $2 ]]; then
    echo "Usage: $0 (dirname) (gaussian script)" 1>&2
    exit 1
fi

rundir=$1
input=$2

remotedirbase=LARGE0.shun/g09queue
queue=gr10341b
script=${input:t}
dir=$remotedirbase/$rundir/${script:r}

thost=automation-laurel01.kudpc.kyoto-u.ac.jp
subg09='~/'$remotedirbase/subg09
subopt="-W 168:00"

perform() {
    cmd=$1
    while true; do
        ssh $thost /bin/zsh -i -l <<< "set -x; mkdir -p $dir; cd $dir; $cmd"
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
        perform "while true; do sleep 300; qjobs | grep -q \"^$jobid \" && continue; exit 1; done"
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

ssh $thost 'mkdir -p '$dir'; cd '$dir'; awk "!x {print \"%nprocShared=32\n%mem=80GB\"; x=1} 1" > '$script < $input

if [[ -n $EXTRA ]]; then
    FILES=($(tr ',' ' ' <<< $EXTRA))
    scp $FILES $thost:$dir
fi

if [[ -n $COPY_ONLY ]]; then
    exit 0
fi

JOBID=$(perform "unsetopt share_history; { set -x; $CMD; module load gaussian09/e01; $subg09 $queue $script $subopt |& tee jobid.txt } |& tee run.log | grep \"[0-9]\\+\\.jb\"")
if [[ -n $NOWAIT ]]; then
    exit 0
fi
#JOBID=$(tr '<>' '  ' <<< $JOBID | cut -f3 -d' ')

wait_completion $JOBID


