#!/usr/bin/zsh
# -*- Shell-script -*-

if [[ -z $2 ]]; then
    echo "Usage: $0 (rundir) (output tmol file) [(rundir) (output tmol file) [...]]" 1>&2
    exit 1
fi

thost=automation-ccfep.ims.ac.jp
remotedirbase=save/turbomolequeue

while true; do
    rundir=$1
    output=$2
    if [[ -z $rundir ]]; then
	break
    fi
    shift
    shift
    
    dir=$remotedirbase/$rundir
    echo get "$dir/coord" "$output"
done | sftp $thost

