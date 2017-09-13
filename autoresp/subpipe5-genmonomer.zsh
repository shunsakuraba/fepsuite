#!/usr/bin/zsh

if [[ -z $1 ]]; then
    echo "Usage: $0 (base-only structure)" 1>&2
    echo "Example: $0 inosine.pdb" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure=$1

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

source $basedir/defaults.zsh

basestructuretype=${basestructure:e}
if [[ $basestructuretype = "gau" ]]; then
    echo "Base structure file name must not be gau" 1>&2
    exit 1
fi

basestructure_save=$basestructure

PREPGEN=${AMBERHOME}/bin/prepgen
ANTECHAMBER=${AMBERHOME}/bin/antechamber

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

#for state in RNA DNA; do
for state in DNA; do
    case $state in
	RNA)
	    basestructure=$basestructure_save
	    backbone=$basedir/RNA-capped.pdb
	    chargemonomerfile=$basedir/rna_monomer_charges.txt
	    baseconstraint=$basedir/constraint-rna.txt
	    ;;
	DNA)
	    basestructure=${basestructure_save:r}.DNA.pdb
	    basestructuretype=pdb
	    backbone=$basedir/DNA-capped.pdb
	    chargemonomerfile=$basedir/dna_monomer_charges.txt
	    baseconstraint=$basedir/constraint-dna.txt
	    ;;
    esac

    basestructurename=${basestructure:r}
    basename=${basestructure:t:r}

    RESNAME=$(grep '^ATOM  ' $basestructure | head -1 | cut -c18-20)
    echo RES:$RESNAME

    # ANTECHAMBER overwrites atomtypes, prevent it with -ao type or -j 0
    #$ANTECHAMBER -i $basestructurename.monomer.ac -fi ac -o $basestructurename.monomer.prep -fo prepi -a $basestructurename.monomer.ac -fa ac -ao type

    # generate monomoer structure. monomer has CH3 terminus at 5'- region
    python $basedir/merger.py ${basestructure_save:r}.base.pdb $backbone $basestructurename.monomer.pdb

    #zsh $basedir/g09fetch.zsh $basename $basename.resp.log $basename.resp.log

    #$ANTECHAMBER -i $basestructurename.final.prep -fi prepi -o $basestructurename.final.mol2 -fo mol2 -j 0

    python $basedir/copycharge.py $basestructurename.monomer.pdb $basestructurename.resp.mol2 $chargemonomerfile $basestructurename.monomer.charge

    MONOMERAC=$basestructurename.monomer.ac

    $ANTECHAMBER -i $basestructurename.monomer.pdb -fi pdb -o $MONOMERAC -fo ac -at amber -c rc -cf $basestructurename.monomer.charge -rn $RESNAME

    # doing it by prepgen
    $PREPGEN -i $basestructurename.monomer.ac -o $basestructurename.monomer.prep -f int -rn $RESNAME && echo

    optmonomerprepmm=${basestructurename}.monomer.opt.xyz
    optmonomeroptlog=${basestructurename}.monomer.opt.log
    optmonomertmol=${basestructurename}.monomer.opt.tmol
    optmonomerdihf=${basestructurename}.monomer.opt.dih
    optmonomerfintmol=${basestructurename}.monomer.optfin.tmol
    optmonomerfinmol2=${basestructurename}.monomer.optfin.mol2
    rundir=${basename}/optmonomer

    # optimize monomer structure
    $OBMINIMIZE -ff UFF $basestructurename.monomer.pdb > $optmonomerprepmm 2> $optmonomeroptlog || true
    if grep -q "CONVERGED" $optmonomeroptlog; then
	true
    else
	echo "Monomer structure optimization failed to converge"
	false
    fi

    python $basedir/constrain-tm.py $optmonomerprepmm $basestructure $baseconstraint $optmonomertmol $optmonomerdihf

    # PBE/TZVPP (as in Zgarbova et al. 2011, for first structure optimization)
    zsh $basedir/turborun.zsh $rundir TITLE="Optimization" CONSTRAINTS=$optmonomerdihf BASIS="TZVPP" CHARGE=$CHARGE FUNCTIONAL=pbe GRID=m4 COSMO=78.4 ENERGY_CONV=6 COORD_CONV=3 SCF_CONV=7 CYCLE=1000 $optmonomertmol

    zsh $basedir/turbofetch.zsh $rundir $optmonomerfintmol

    $OPENBABEL $optmonomerfintmol $optmonomerfinmol2

    # generate correctly named, correct structure, and charged mol2
    optmonomer=$basestructurename.monomer.opt.mol2
    optmonomerprep=$basestructurename.monomer.opt.prep
    optmonomerpdb=$basestructurename.monomer.opt.pdb
    $ANTECHAMBER -i $MONOMERAC -fi ac -o $optmonomer -fo mol2 -j 0 -a $optmonomerfinmol2 -fa mol2 -ao crd
    $ANTECHAMBER -i $MONOMERAC -fi ac -o $optmonomerpdb -fo pdb -j 0 -a $optmonomerfinmol2 -fa mol2 -ao crd
    $ANTECHAMBER -i $MONOMERAC -fi ac -o $optmonomerprep -fo prepi -j 0 -a $optmonomerfinmol2 -fa mol2 -ao crd
done

#$ANTECHAMBER -i $MONOMERAC -fi ac -j 0 -o $basestructurename.optstructure.pdb -fo pdb

unsetopt -x
