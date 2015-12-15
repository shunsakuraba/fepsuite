#!/usr/bin/zsh

# submit g09 job to magnolia

if [[ -z $2 ]]; then
    echo "Usage: $0 (dirname) (gaussian script)" 1>&2
    exit 1
fi

rundir=$1
input=$2

remotedirbase=LARGE2.asai/g09queue
queue=gr10341d
dir=$remotedirbase/$rundir

script=${input:t}

subg09='~/'$remotedirbase/subg09

trap "echo Error \$? returned while sending command 1>&2 ; exit 1" ZERR
set -x

ssh magnolia.kudpc.kyoto-u.ac.jp 'mkdir -p '$dir'; cd '$dir'; awk "/^%/ && !x {print \"%nprocShared=28\"; x=1} 1" > '$script < $input
ssh -t magnolia.kudpc.kyoto-u.ac.jp << EOF
cd $dir
{
  source /etc/zshrc
  set -x
  export MODULEPATH=/opt/app/modulefiles/app_ISV:\$MODULEPATH
  module load gaussian09/d01
  $subg09 $queue $script -K
} |& tee run.log
EOF


