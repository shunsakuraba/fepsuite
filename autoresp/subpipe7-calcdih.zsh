#!/usr/bin/zsh

if [[ -z $1 ]]; then
    echo "Usage: $0 (base-only structure)" 1>&2
    echo "Example: $0 inosine.pdb" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure_save=$1

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

declare -A methmap basismap
declare -A multiplicity

#NOWAIT=y zsh $basedir/turborun.zsh $rundir TITLE="Optimization" CONSTRAINTS=$optdihf BASIS="cc-pVTZ" CHARGE=$CHARGE FUNCTIONAL=pbe GRID=m4 COSMO=78.4 ENERGY_CONV=6 COORD_CONV=3 SCF_CONV=7 CYCLE=1000 NEWBASIS=y $opttmol



#methmap=(mp2pvdz "MP2(SemiDirect)/cc-pVDZ"  mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ") # MP2/cc-pVDZ is no longer necessary since CCSD(T) includes it
#methmap=(mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ")
methmap=(mp2pvtz mp2  mp2pvqz mp2  pbelp pbe)
basismap=(mp2pvtz "cc-pVTZ"  mp2pvqz "cc-pVQZ"  pbelp "6-311++G(3df,3pd)")
#multiplicity=(mp2pvdz 2  mp2pvtz 4  mp2pvqz 3  ccsdpvdz 12)

DIRS=()
RECLAIM=()
touch dummy.constraint
for state in RNA DNA; do
    case $state in
	RNA)
	    basestructure=$basestructure_save
	    # See subpipe6.
	    range=({0..3} {15..35})
	    ;;
	DNA)
	    basestructure=${basestructure_save:r}.DNA.pdb
	    range=({3..15})
	    ;;
    esac
    
    basestructurename=${basestructure:r}
    basename=${basestructure:t:r}

    for i in $range; do
	for k in ${(k)methmap}; do
	    meth=${methmap[$k]}
	    basis=${basismap[$k]}

	    tmol_prev=$basestructurename.dihopt$i.done.tmol
	    tmol=$basestructurename.dihopt$i.done_clear.tmol
	    # This looks stupid, but turbomole may fail to define redundant coodinate.
	    $OPENBABEL $tmol_prev $tmol
	    case $k in
		pbelp)
		    NEWBASIS=y
		    SCRIPT=run.csh
		    COSMO=78.4
		    # For SP calculation BOTH COSMO and BJ dumping are needed
		    DISP3=bj
		    ;;
		*)
		    NEWBASIS=
		    SCRIPT=rimp2.csh
		    COSMO=
		    DISP3=
		    ;;
	    esac
	    rundir=$basestructurename/sp.$k.dih$i
# SCF_CONV=8 is required for some mp2 calculations.
	    NOWAIT=y zsh $basedir/turborun.zsh $rundir TITLE="$k" CONSTRAINTS=dummy.constraint BASIS=$basis CHARGE=$CHARGE FUNCTIONAL=pbe GRID=m4 COSMO=$COSMO ENERGY_CONV=6 COORD_CONV=3 SCF_CONV=8 CYCLE=0 NEWBASIS=$NEWBASIS SCRIPT=$SCRIPT DISP3=$DISP3 $tmol
	    DIRS+=$rundir

	    RECLAIM+=($rundir,$basestructurename.dihcalc$i.$k.log)
	done
	sleep 60
    done
done

ENERGYONLY=Y zsh $basedir/turbowait.zsh $DIRS

for f in $RECLAIM; do
    ENERGYONLY=Y zsh $basedir/turbofetch.zsh ${f%%,*} ${f##*,}
    sleep 60
done


