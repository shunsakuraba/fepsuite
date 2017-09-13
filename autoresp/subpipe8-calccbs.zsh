#!/usr/bin/zsh

# generate prep file for 3'- and 5'- termini 
if [[ -z $1 ]]; then
    echo "Usage: $0 (base-only structure) " 1>&2
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

declare -A methmap

#methmap=(mp2pvdz "MP2(SemiDirect)/cc-pVDZ"  mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ") # MP2/cc-pVDZ is no longer necessary since CCSD(T) already includes it
methmap=(mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ")

#fetches=()
#for k in ${(k)methmap}; do
#    for i in {0..35}; do
#	runlog=$basestructurename.dihcalc.$k.$i.log
#	locallog=$runlog
#	fetches+=($basename $runlog $locallog)
#    done
#done
#zsh $basedir/g09fetch.zsh $fetches

OUTF=${basestructure_save:r}.dihcalc.mp2cbs.log
mv -f $OUTF $OUTF.bak || true

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

    dihsummary=$basestructurename.dihsummary.txt
    for i in $range; do
	echo -n "$state $i " >> $OUTF
	python $basedir/calcCBS_mp2_tm.py $basestructurename.dihcalc$i.%s.log >> $OUTF
    done
done
