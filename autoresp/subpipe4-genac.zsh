#!/usr/bin/zsh
echo "FIXME THIS SCRIPT IS INCOMPLETE (file names should be replaced correctly)"
exit 1

# Resp charge fitting final phase
# Assumes subpipe2-opt is done
if [[ -z $4 ]]; then
    echo "Usage: $0 (base-only structure) (backbone) (charge overwrite file) (charge monomer file) (mainchain specification)" 1>&2
    echo "Example: $0 inosine.pdb RNA.pdb rna_charges.txt rna_monomer_charges.txt rna.mainchain" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure=$1
backbone=$2
chargefile=$3
mainchain=$4

basestructuretype=${basestructure:e}
if [[ $basestructuretype = "gau" ]]; then
    echo "Base structure file name must not be gau" 1>&2
    exit 1
fi
basestructurename=${basestructure:r}
basename=${basestructure:t:r}

PREPGEN=${AMBERHOME}/bin/prepgen
ANTECHAMBER=${AMBERHOME}/bin/antechamber

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

ACFILE=$basestructurename.final.ac
python $basedir/merger.py $basestructurename.base.pdb $backbone $basestructurename.full.pdb 

python $basedir/copycharge.py $basestructurename.full.pdb $basestructurename.resp.mol2 $chargefile $basestructurename.final.charge

$ANTECHAMBER -i $basestructurename.full.pdb -fi pdb -o $ACFILE -fo ac -at amber -c rc -cf $basestructurename.final.charge

# special treatment for RNA, what to do in general??
sed -i "/O3\'/s/ O$/OS/" $ACFILE || true

RESNAME=$(grep '^ATOM' $ACFILE | head -1 | cut -c18-20)

$PREPGEN -i $ACFILE -o $basestructurename.final.prep -m rna.mainchain -f int -rn $RESNAME && echo

#---------------------------------------------------
# Generate dot file for visual inspection
python $basedir/prep2dot.py < $basestructurename.final.prep > $basestructurename.final.dot
neato $basestructurename.final.dot -Tps -o $basestructurename.final.ps


RESSHORT=${RESNAME%?}

#---------------------------------------------------
# Generate 5' and 3' terminus version (5'-OH, 3'-OH)

for terminus in 3 5; do
    ACFILE=$basestructurename.$terminus.ac
    python $basedir/merger.py $basestructurename.base.pdb ${backbone:r}$terminus.${backbone:e} $basestructurename.$terminus.pdb 

    python $basedir/copycharge.py -nonint $basestructurename.$terminus.pdb $basestructurename.resp.mol2 ${chargefile:r}$terminus.${backbone:e} $basestructurename.$terminus.charge

    $ANTECHAMBER -i $basestructurename.$terminus.pdb -fi pdb -o $ACFILE -fo ac -at amber -c rc -cf $basestructurename.$terminus.charge

    # hack
    sed -i "/O3\'/s/ O$/OS/" $ACFILE || true

    $PREPGEN -i $ACFILE -o $basestructurename.$terminus.prep -m rna$terminus.mainchain -f int -rn ${RESSHORT}$terminus && echo

done


unset -x

echo "*********************************"
echo "Subpipe 4 finished"
echo "Check $basestructurename.final.ps, and (if necessary) manually update atomtypes of $basestructurename.final.prep / $basestructure.monomer.prep."
echo "*********************************"

