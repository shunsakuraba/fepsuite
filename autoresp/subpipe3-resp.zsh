#!/usr/bin/zsh

# Resp charge fitting.
# Assumes subpipe2-opt is done
if [[ -z $0 ]]; then
    echo "Usage: $0" 1>&2
    echo "Example: $0" 1>&2
    exit 1
fi

choke me turn this script to derive both RNA / DNA version

source vars.zsh
basedir=${0:h}

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

RESP=${AMBERHOME}/bin/resp
ESPGEN=${AMBERHOME}/bin/espgen
RESPGEN=${AMBERHOME}/bin/respgen

basestructure_save=$basestructurename
basestructure=$basestructurename.pdb

basestructuretype=${basestructure:e}
if [[ $basestructuretype = "gau" ]]; then
    echo "Base structure file name must not be gau" 1>&2
    exit 1
fi

for iter in 1 2; do
    case $iter in
    1)
        basestructure=$basestructure_save.pdb
        basestructurename=${basestructure_save:r}
        basename=${basestructure_save:t:r}
        baserestraint=$basedir/resp-rna.txt
    ;;
    2)
        basestructure=$basestructure_save.DNA.pdb
        basestructurename=${basestructure_save:r}.DNA
        basename=${basestructure_save:t:r}.DNA
        baserestraint=$basedir/resp-dna.txt
    ;;
    esac

    optfintmol=${basestructurename}.opt2fin.tmol
    respgau=${basestructurename}.resp.gau
    respcheck=${basename}.resp.chk

    $OPENBABEL $optfintmol -o gau -xk "%chk=$respcheck
#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) sp" $respgau

    sed "4c Comment line" -i $respgau

    sed "6c $CHARGE 1" -i $respgau

    respgauout=${basestructurename}.resp.gau.out
    zsh $basedir/g09rccs.zsh $basename $respgau
    zsh $basedir/g09rccsfetch.zsh $basename $respgauout

    ACFILE=$basestructurename.resp.ac

    $ANTECHAMBER -i $basestructure -fi pdb -o $ACFILE -fo ac -at amber

    ESP=$basestructurename.esp

    # run espgen to extract electron cloud distribution
    $ESPGEN -i $respgauout -o $ESP

    # generate respgen input (input of input)
    RESPADDIN=$basestructurename.respadd
# FIXME TODO: respgen output discards precision. This problem is (miraclously) avoided in current condition, but we need to fix...
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

done

# generate monomer without phos group
# This part is now performed in genac
#$ANTECHAMBER -i $basestructurename.gesp -fi gesp -o $basestructurename.monomer.pre.ac -fo ac -at amber -c resp -gv 1 -nc $CHARGE
#$ANTECHAMBER -i $basestructurename.monomer.pre.ac -fi ac -a $basestructure -fa pdb -ao name -o $basestructurename.monomer.ac -fo ac -at amber -rn $RESNAME

