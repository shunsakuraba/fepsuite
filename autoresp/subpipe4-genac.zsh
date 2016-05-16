#!/usr/bin/zsh
echo "FIXME THIS SCRIPT IS INCOMPLETE (file names should be replaced correctly)"
exit 1
# Resp charge fitting final phase
# Assumes subpipe2-opt is done
if [[ -z $4 ]]; then
    echo "Usage: $0 (base-only structure) (backbone) (charge overwrite file) (mainchain specification)" 1>&2
    echo "Example: $0 inosine.pdb RNA.pdb rna_charges.txt rna.mainchain" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure=$1
basecharge=$2
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
python $basedir/merger.py $basestructure RNA.pdb $basestructurename.full.pdb 

$ANTECHAMBER -i $basestructurename.full.pdb -fi pdb -o $ACFILE -fo ac -at amber -c rc -cf $basestructurename.final.charge

python copycharge.py inosine.full.pdb inosine.resp.mol2 $basecharge inosine.final.charge

# special treatment for RNA, what to do in general??
sed -i "/O3\'/s/ O$/OS/" $ACFILE || true
RESNAME=$(grep '^ATOM' $ACFILE | head -1 | cut -c18-20)

$PREPGEN -i $ACFILE -o $basestructurename.final.prep -m rna.mainchain -f int -rn $RESNAME && echo
$PREPGEN -i $ACFILE -o $basestructurename.monomer.prep -f int -rn $RESNAME && echo

python $basedir/prep2dot.py < $basestructurename.final.prep > $basestructurename.final.dot
neato $basestructurename.final.dot -Tps -o $basestructurename.final.ps

echo "*********************************"
echo "Subpipe 4 finished"
echo "Check $basestructurename.final.ps, and (if necessary) manually update atomtypes of $basestructurename.final.prep / $basestructure.monomer.prep."
echo "*********************************"