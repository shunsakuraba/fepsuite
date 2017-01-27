#!/bin/zsh

# Assumes subpipe9-mmdihoptcalc is done
if [[ -z $4 ]]; then
    echo "Usage: $0 (base-only structure) (base constraint) (dihedral) (atomtype) (resname)" 1>&2
    echo "Example: $0 inosine.pdb constraint-rna-opt2.txt O4\'-C1\'-N9-C8 CP INO" 1>&2
    echo "The final atom in dihedral angle specification will be given the new atom type"
    exit 1
fi

basedir=${0:h}
basestructure=$1
baseconstraint=$basedir/$2
dihedral=$3
newatomtype=$4
resname=$5

basestructurename=${basestructure:r}
basename=${basestructure:t:r}

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

for i in {0..35}; do
    runlog=$basestructurename.optmm.$i.log
    locallog=$runlog
    fetches+=($basename $runlog $locallog)
done
#zsh $basedir/g09fetch.zsh $fetches

qme=$basestructurename.dihcalc.ccsdtcbs.log
mme=$basestructurename.mm.txt
rm -f $mme
for i in {0..35}; do
    MME=$(grep '^\s\+Energy=' $basestructurename.optmm.$i.log | tail -1 | awk '{print $2}')
    if (( $? != 0 )); then
	echo "Error while extracting MM-energy" 1>&2
	exit 1
    fi
    echo $MME >> $basestructurename.mm.txt
done

weight=$basestructurename.weight.txt
for i in {0..35}; do
    # double weight from 180-220 and 240-280.
    # note that chi angle definition in OL3 differes from angle definition by 180 degree 
    if (( (i >= 0 && i <= 4) || (i >= 6 && i <= 10) )); then
	echo 2
    else
	echo 1
    fi
done > $weight


result=$basestructurename.funopt.txt
python $basedir/funopt.py $qme $mme $weight $result --hartree
frcmodbase=$basestructurename.dihopt.frcmod
frcmodopt=$basestructurename.frcmod
python $basedir/replace_frcmod_opt_target.py $result $frcmodbase > $frcmodopt

for suf in final 3 5; do
    python $basedir/replace_prep.py $basestructurename.$suf.prep $dihedral $newatomtype > $basestructurename.$suf.tmp.prep
    # bsc0
    K=.$suf
    if [[ $suf = final ]]; then
	K=
    fi
    python $basedir/replace_prep.py $basestructurename.$suf.tmp.prep C2\'-C3\'-C4\'-C5\' CI > $basestructurename.opted$K.prep
done

