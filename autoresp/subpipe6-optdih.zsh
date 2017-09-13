#!/usr/bin/zsh

# Make a structure with 5'-OH replaced with H (H53)

# Assumes subpipe4-genac is done
if [[ -z $3 ]]; then
    echo "Usage: $0 (base-only structure) (dihedral) (atomtype)" 1>&2
    echo "Example: $0 inosine.pdb O4\'-C1\'-N9-C8 CP" 1>&2
    echo "The final atom in dihedral angle specification (e.g. C8) will be assigned a new atom type (e.g. CP)"
    exit 1
fi

basedir=${0:h}
basestructure_save=$1
dihedral=$2
newatomtype=$3

#resname=$(grep '^ATOM' $basestructure_save | head -1 | cut -c18-20)

set -x
trap "{ echo 'Error detected with code =' $? 2>&1 ; exit 1 }" ZERR

diheds=($(tr '-' ' ' <<< $dihedral))
repatm=${diheds[4]}

atomtype() {
    targ=$1
    prep=$2
    awk 'substr($0, 7, 4) == sprintf("%-4s", "'$targ'") { a=substr($0, 13, 4); gsub(/ /, "", a); print a }' $prep
}

RECLAIM=()
for state in RNA DNA; do
    case $state in
	RNA)
	    basestructure=$basestructure_save
	    baseconstraint=$basedir/constraint-rna-opt2.txt
	    # See below. We intentionally compute the overlap (3, 15) so that we can align them easily.
	    range=({0..3} {15..35})
	    ;;
	DNA)
	    basestructure=${basestructure_save:r}.DNA.pdb
	    baseconstraint=$basedir/constraint-dna-opt2.txt
	    # DNA is used "between 210 and 330" (Zgarbova 2013)
            # Our chi angle definition is chi + 180, so:
	    # DNA is used "between 30 and 150".
	    range=({3..15})
	    ;;
    esac

    basestructurename=${basestructure:r}
    basename=${basestructure:t:r}

    monomerprep=$basestructurename.monomer.opt.prep
    monomermol2=$basestructurename.monomer.opt.mol2
    monomerpdb=$basestructurename.monomer.opt.pdb


    if [[ ! -e $monomerprep ]] || [[ ! -e $monomermol2 ]] || [[ ! -e $monomerpdb ]]; then
	echo "Error: monomer prep file is not available" 1>&2
	exit 1
    fi

    prevatomtype=$(atomtype $repatm $monomerprep)

    # generate modified prep file
    dihoptprep=$basestructurename.dihopt.prep
    dihoptfrcmod=$basestructurename.dihopt.frcmod
    baseparm=$AMBERHOME/dat/leap/parm/parm10.dat

    python $basedir/copydih.py $monomerprep $dihedral $newatomtype $baseparm > $dihoptfrcmod
    python $basedir/replace_prep.py $monomerprep $dihedral $newatomtype > $dihoptprep

    if [[ ! -e $basestructurename.extra.frcmod ]]; then
	echo "$basestructurename.extra.frcmod does not exist, making dummy file" 2>&1
	echo "(Wating for 10 secs)" 2>&1
	sleep 10
	touch $basestructurename.extra.frcmod
    fi

    # generate amber topology
    leap_templ=$basedir/modrna.leap.in
    dihoptleap=$basestructurename.dihopt.leap

    sed "s/%base%/$basestructurename/g" < $leap_templ > $dihoptleap

    TLEAP=$AMBERHOME/bin/tleap

    $TLEAP -f $dihoptleap
    if [[ ! -e $basestructurename.dihopt.ambtop ]]; then
	echo "Failed to generate ambtop" 1>&2
	false
    fi

    # dirty hack: revert exchanged atomtype back only in ambtop
    # This is necessary to set correct PB radii in pbsa
    # Say f%!k to the hard-coded PB radii
    newatomtype=$(atomtype $repatm $dihoptprep)
    sed -i '/%FLAG AMBER_ATOM_TYPE/,/%FLAG/ {s/'$newatomtype'  /'$prevatomtype'  /}' $basestructurename.dihopt.ambtop
    # remaining MM thingies are in subpipe9

    # Run TMOL optimization
    for i in $range; do
	opttmol=$basestructurename.dihopt$i.tmol
	optdihf=$basestructurename.dihopt$i.txt
	{
	    cat $baseconstraint
	    echo "dihedral $diheds "$((i * 10))
	} > $basestructurename.cons$i.txt
        # monomermol2 should be properly named, no need to hassle
	python $basedir/constrain-tm.py $monomerpdb $monomermol2 $basestructurename.cons$i.txt $opttmol $optdihf
    done

    DIRS=()
    for i in $range; do
	opttmol=$basestructurename.dihopt$i.tmol
	optdihf=$basestructurename.dihopt$i.txt
	rundir=$basestructurename/optdih$i
	DIRS+=$rundir
        # Optimization is performed with pbe/LP/COSMO. No BJ3 dumping (see p.2891)
	NOWAIT=y zsh $basedir/turborun.zsh $rundir TITLE="Optimization" CONSTRAINTS=$optdihf BASIS="6-311++G(3df,3pd)" CHARGE=$CHARGE FUNCTIONAL=pbe GRID=m4 COSMO=78.4 ENERGY_CONV=6 COORD_CONV=3 SCF_CONV=7 CYCLE=1000 NEWBASIS=y $opttmol
	sleep 60
	RECLAIM+=($basestructurename/optdih$i,$basestructurename.dihopt$i.done.tmol)
    done
done

zsh $basedir/turbowait.zsh $DIRS

for f in $RECLAIM; do
    zsh $basedir/turbofetch.zsh ${f%%,*} ${f##*,}
    sleep 60
done

