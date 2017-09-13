#!/bin/zsh

# Assumes subpipe9-mmdihoptcalc is done
if [[ -z $3 ]]; then
    echo "Usage: $0 (base-only structure) (dihedral) (atomtype)" 1>&2
    echo "Example: $0 inosine.pdb O4\'-C1\'-N9-C8 CP" 1>&2
    echo "The final atom in dihedral angle specification will be given the new atom type"
    exit 1
fi

basedir=${0:h}
basestructure_save=$1
dihedral=$2
newatomtype=$3

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

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
	runlog=$basestructurename.optmm.$i.log
	locallog=$runlog
	fetches+=($basename $runlog $locallog)
    done
done
zsh $basedir/g09rccsfetch.zsh $fetches

basestructurename=${basestructure_save:r}
qme=$basestructurename.dihcalc.mp2cbs.log
mme=$basestructurename.mm.txt
rm -f $mme

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
	MME=$(grep '^\s\+Energy=' $basestructurename.optmm.$i.log | tail -1 | awk '{print $2}')
	if (( $? != 0 )); then
	    echo "Error while extracting MM-energy" 1>&2
	    exit 1
	fi
	echo $state $i $MME >> $mme
    done
done

python $basedir/plot-energy.py $mme mm-energy.png
python $basedir/plot-energy.py $qme qm-energy.png

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

# Optimization done, generate RNA params

basestructure=$basestructure_save
basestructurename=${basestructure:r}
basename=${basestructure:t:r}


frcmodbase=$basestructurename.dihopt.frcmod
frcmodopt=$basestructurename.frcmod
python $basedir/replace_frcmod_opt_target.py $result $frcmodbase > $frcmodopt

for suf in final 3 5; do
    python $basedir/replace_prep.py $basestructurename.$suf.prep $dihedral $newatomtype > $basestructurename.$suf.tmp.prep
    # apply bsc0
    K=.$suf
    if [[ $suf = final ]]; then
	K=
    fi
    python $basedir/replace_prep.py $basestructurename.$suf.tmp.prep C2\'-C3\'-C4\'-C5\' CI > $basestructurename.opted$K.prep
done

