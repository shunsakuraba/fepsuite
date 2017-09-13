#!/usr/bin/zsh
# -*- Shell-script -*-

if [[ -z $2 ]]; then
    echo "Usage: $0 (rundir) (output tmol file)" 1>&2
    exit 1
fi

rundir=$1
output=$2

thost=automation-ccfep.ims.ac.jp
remotedirbase=save/turbomolequeue
dir=$remotedirbase/$rundir

if [[ $ENERGYONLY = Y ]]; then
    scp $thost:$dir/'turbo.log' $output
else
    ssh $thost "[ -e $dir/GEO_OPT_CONVERGED ]" || { echo "tmole not converged" && exit 1 }
    scp $thost:$dir/coord $output
fi

