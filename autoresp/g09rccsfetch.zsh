#!/usr/bin/zsh

# fetch files from magnolia

if [[ -z $2 ]]; then
    echo "Usage: $0 (dirname) (remote file) (local file to save)" 1>&2
    exit 1
fi

thost=ccfep.ims.ac.jp
remotedirbase=save/g09queue

while true; do
    rundir=$1
    input=$2
    output=$3

    if [[ -z $output ]]; then
        if [[ -z $input ]]; then
	    break
        else
            output=$input
        fi
    else
        shift || break
    fi

    shift || break
    shift || break

    script=${log:-$input}
    script=${script:t}
    script=${script%%.out}
    dir=$remotedirbase/$rundir/${script:r}

    echo get $dir/$input ./$output
done | sftp $thost 
