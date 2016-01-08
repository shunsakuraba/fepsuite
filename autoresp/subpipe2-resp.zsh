#!/usr/bin/zsh

# Resp charge fitting.
# Assumes subpipe2-opt is done
if [[ -z $2 ]]; then
    echo "Usage: $0 (structure) (charge restraint file)" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure=$1
baserestraint=$2

RESP=${AMBERHOME}/bin/resp
ESPGEN=${AMBERHOME}/bin/espgen
RESPGEN=${AMBERHOME}/bin/respgen

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

basestructuretype=${basestructure:e}
if [[ $basestructuretype = "gau" ]]; then
    echo "Base structure file name must not be gau" 1>&2
    exit 1
fi
basestructurename=${basestructure:r}
basename=${basestructure:t:r}

respgau=${basename}.resp.gau
respcheck=${basename}.resp.chk
opt4check=${basename}.opt4.chk

# RESP gaussian
cat >> $respgau <<EOF
--Link1--
%oldchk=$opt4check
%chk=$respcheck
#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) sp
# iop(6/50=1) GEOM=Checkpoint

RESP calculation

$CHARGE 1

$basename.gesp

$basename.gesp


EOF

#zsh $basedir/g09run.zsh $basename $respgau
#zsh $basedir/g09fetch.zsh $basename $basename.gesp $basestructurename.gesp

ACFILE=$basestructurename.resp.ac

$ANTECHAMBER -i $basestructure -fi pdb -o $ACFILE -fo ac -at amber

ESP=$basestructurename.esp

# run espgen to extract electron cloud distribution
$ESPGEN -i $basestructurename.gesp -o $ESP

# generate respgen input (input of input)
RESPADDIN=$basestructurename.respadd
python $basedir/gen_respgen.py $basestructure $baserestraint > $RESPADDIN

RESPIN1=$basestructurename.respin1
RESPIN2=$basestructurename.respin2
RESPOUT1=$basestructurename.respout1
RESPOUT2=$basestructurename.respout2
QOUT1=$basestructurename.qout1
QOUT2=$basestructurename.qout2
# generate resp input
$RESPGEN -i $ACFILE -o $RESPIN1 -a $RESPADDIN -f resp1
$RESPGEN -i $ACFILE -o $RESPIN2 -a $RESPADDIN -f resp2

# run resp
$RESP -O -i $RESPIN1 -o $RESPOUT1 -e $ESP -t $QOUT1
$RESP -O -i $RESPIN2 -o $RESPOUT2 -e $ESP -t $QOUT2 -q $QOUT1

MOLRES=$basestructurename.resp.mol2
$ANTECHAMBER -i $basestructure -fi pdb -o $MOLRES -fo mol2 -at amber -c rc -cf $QOUT2
