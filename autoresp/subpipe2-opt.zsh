#!/usr/bin/zsh

# Constraint optimization.
# This optimization has two objectives:
# 1. Determine the structure of the molecule for the charge assignment.
# 2. Determine the angle for constraining in later calculations

source vars.zsh || exit 1

if [[ -z $0 ]]; then
    echo "Usage: $0" 1>&2
    echo "Example: $0" 1>&2
    exit 1
fi

basedir=${0:h}
if [[ -z $basestructurename ]]; then
    echo "Unable to find basestructure"
    exit 1
fi
basestructure=$basestructurename.pdb
baseconstraint=${basedir}/constraint-rna.txt
baseconstraint_D=${basedir}/constraint-dna.txt

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

for iter in 1 2; do
    case $iter in
    1)
    runstructure=$basestructure
    runconstraint=$baseconstraint
    ;;
    2)
    runstructure=$basestructurename.DNA.pdb
    runconstraint=$baseconstraint_D
    ;;
    esac
    basestructuretype=${runstructure:e}
    if [[ $basestructuretype = "gau" ]]; then
        echo "Base structure file name must not be gau" 1>&2
        exit 1
    fi
    basestructurename=${runstructure:r}
    basename=${runstructure:t:r}

    initialgau=${basestructurename}.init.gau

    opt1gau=${basestructurename}.opt1.gau
    opt1check=${basename}.opt1.chk
    opt1log=${basestructurename}.opt1.gau.out

    # Use babel, not antechamber (antechamber may swap atoms and emits no information for recovering atom names!)
    $OPENBABEL $runstructure -ogzmat -xk "%chk=$opt1check
# HF/6-31G POPT(Z-matrix,maxcycle=5000)" $initialgau

    # Fix 'Holmium' problem and broken multiplicity
    sed -i "s/Ho/H/;/^0 /c $CHARGE 1" $initialgau

    python $basedir/mod-zmatrix.py $initialgau $runstructure $runconstraint > $opt1gau

    # Optimize the structure with constraint
    # HF/6-31G level first to remove large crashes
    GMAJOR=09 GMINOR=e01 zsh $basedir/g09rccs.zsh $basename $opt1gau
    zsh $basedir/g09rccsfetch.zsh $basename $opt1log $opt1log

    opt2tmol=${basestructurename}.opt2.tmol
    opt2dihf=${basestructurename}.opt2.dih
    optfintmol=${basestructurename}.opt2fin.tmol
    rundir=${basename}/opt2

    # Switch to Turbomole
    # Define constraints
    python $basedir/constrain-tm.py $opt1log $runstructure $runconstraint $opt2tmol $opt2dihf

    # PBE/TZVPP (as in Zgarbova et al. 2011)
    # BJ3 dumping is not applied during the optimization (see page 2891)
    zsh $basedir/turborun.zsh $rundir TITLE="Optimization" CONSTRAINTS=$opt2dihf BASIS="TZVPP" CHARGE=$CHARGE FUNCTIONAL=pbe GRID=m4 COSMO=78.4 ENERGY_CONV=6 COORD_CONV=3 SCF_CONV=7 CYCLE=1000 $opt2tmol 

    zsh $basedir/turbofetch.zsh $rundir $optfintmol
done
