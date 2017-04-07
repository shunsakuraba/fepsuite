#!/usr/bin/zsh

# generate prep file for 3'- and 5'- termini 
if [[ -z $2 ]]; then
    echo "Usage: $0 (base-only structure) (backbone)" 1>&2
    echo "Example: $0 inosine.pdb RNA.pdb" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure=$1

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

basestructurename=${basestructure:r}
basename=${basestructure:t:r}

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

declare -A methmap

#methmap=(mp2pvdz "MP2(SemiDirect)/cc-pVDZ"  mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ") # MP2/cc-pVDZ is no longer necessary since CCSD(T) already includes it
methmap=(mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ")

fetches=()
for k in ${(k)methmap}; do
    for i in {0..35}; do
	runlog=$basestructurename.dihcalc.$k.$i.log
	locallog=$runlog
	fetches+=($basename $runlog $locallog)
    done
done
zsh $basedir/g09fetch.zsh $fetches

dihsummary=$basestructurename.dihsummary.txt
for i in {0..35}; do
    python $basedir/calcCBS.py $basestructurename.dihcalc.%s.$i.log
done > $basestructurename.dihcalc.ccsdtcbs.log

