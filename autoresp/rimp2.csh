#!/home/users/rq0/cx/lcl/bin/zsh -f
#PBS -l select=1:ncpus=4:mpiprocs=4:ompthreads=1:jobtype=small
#PBS -l walltime=72:00:00
# -*- mode: shell-script; sh-shell-file: zsh -*-

# this is a hack to run zsh instead of csh
set is_csh = 123
test X$is_csh = X123 && exec zsh $0

if [[ -z ${PBS_O_WORKDIR} ]]; then
  # submit itself
  ID=$(jsub -q PF -v FOO=BAR $0)
  echo $ID
  exit 0
fi

cd ${PBS_O_WORKDIR}

MO=$HOME/opt/pg/mo

NP=$(wc -l < $PBS_NODEFILE)

export PARA_ARCH=MPI
export OMP_NUM_THREADS=1
tmoldir=/local/apl/pg/turbomole
source $tmoldir/Config_turbo_env

source ./env.txt

SOURCEDIR=$PWD
WORKDIR=/work/users/rq0/${PWD##/*/users/rq0/}
TEMPLATE=$HOME/save/turbomolequeue/rimp2-sp.in.template

if [[ -e $WORKDIR ]]; then
  rm -rf ${WORKDIR}.bak
  mv ${WORKDIR} ${WORKDIR}.bak
fi
mkdir -p ${WORKDIR}

cp $SOURCEDIR/{env.txt,coord.in} $WORKDIR
if [[ -e $SOURCEDIR/constraints ]]; then
  cp $SOURCEDIR/constraints $WORKDIR
  CONSTRAINTS="$(cat constraints)"
  export CONSTRAINTS
else
  export CONSTRAINTS=dis
fi
cd $WORKDIR

export COORD=coord.in

if [[ -z $NEWBASIS ]]; then
  # Uses default basis set
  export LIBNUM=1
else
  # Uses new basis sets
  export LIBNUM=3
fi

# Default values
if [[ -z $ENERGY_CONV ]]; then
  ENERGY_CONV=6
fi

if [[ -z $COORD_CONV ]]; then
  COORD_CONV=3
fi

if [[ -z $SCF_CONV ]]; then
  # turbomole default is 6, but 7 is preferred almost always
  SCF_CONV=7
fi

if [[ -z $CYCLE ]]; then
  CYCLE=100
fi

# generate define.in
$MO --source=./env.txt $TEMPLATE | grep -v '^#' > define.in
define < define.in

if [[ ! -e mos ]] || [[ ! -e control ]] || [[ ! -e basis ]]; then
  echo "Define Failed!" && exit 1
fi

# modify control
if [[ ! -z $COSMO ]]; then
  sed -i -e '/^\$end$/i$cosmo\
epsilon='$COSMO control
fi

if [[ ! -z $DISP3 ]]; then
  sed -i -e '/^\$end$/i$disp3 '$DISP3 control
fi

# Done with preparation, running now!
sed -i '/^\$end/i $denconv='$ENERGY_CONV control

# rimp2prep part is tricky due to interactive input.
# save here for the data
natomtype=$(grep -v '^\$' $COORD | wc -l)
(( natomtype = natomtype - 1 ))

dscf 
ECODE=$?
if [[ $ECODE != 0 ]]; then
  echo "DSCF failed!" && exit 1
fi

# rimp2prep part is tricky due to interactive input.
{
  # energy
  echo e
  # just finish
  echo '*'
  # auxbasis (all default)
  for i in {1..$natomtype}; do
    echo
  done
  # finish
  echo '*'
} | rimp2prep

# rimp2 is bugged. Use ricc2 instead. Happy working-around!
sed -i '/^\$end/i $ricc2\
 mp2 energy only
' control
ricc2 > turbo.log
ECODE=$?

cp $WORKDIR/* $SOURCEDIR
echo $ECODE > $SOURCEDIR/done.txt
return $?

# vim:filetype=zsh
