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

basestructurename=${basestructure:r}
basename=${basestructure:t:r}

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

declare -A methmap
declare -A multiplicity

#methmap=(mp2pvdz "MP2(SemiDirect)/cc-pVDZ"  mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ") # MP2/cc-pVDZ is no longer necessary since CCSD(T) includes it
methmap=(mp2pvtz "MP2(SemiDirect)/cc-pVTZ"  mp2pvqz "MP2(SemiDirect)/cc-pVQZ"  ccsdpvdz "CCSD(T)/cc-pVDZ")
multiplicity=(mp2pvdz 2  mp2pvtz 4  mp2pvqz 3  ccsdpvdz 12)

N=36
for k in ${(k)methmap}; do
    meth=${methmap[$k]}
    sep=${multiplicity[$k]}
    for iter in {1..$sep}; do
	(( beg = (iter - 1) * N / sep )) || true
	(( end = iter * N / sep - 1 )) || true
	for i in {$beg..$end}; do
	# for i in 1; do
	    tmol=$basestructurename.dihopt$i.done.tmol
	    rungau=$basestructurename.dihcalc.$k.$i.gau
	    chk=$basename.dihcalc.$k.$i.chk
	    $OPENBABEL $tmol -xk "#SP $meth SCF=VeryTight" $rungau
	    sed -i "1i%chk=$chk\\

3c \ Remark

5c \ $CHARGE 1" $rungau
            # Laurel is unstable. Copy only
	    COPY_ONLY=y zsh $basedir/g09run.zsh $basename $rungau
	done &
	sleep 60
    done
done
wait


