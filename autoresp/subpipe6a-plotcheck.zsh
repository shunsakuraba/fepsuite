#!/usr/bin/zsh

# Make a structure with 5'-OH replaced with H (H53)

# Assumes subpipe4-genac is done
if [[ -z $1 ]]; then
    echo "Usage: $0 (base-only structure)"  1>&2
    echo "Example: $0 inosine.pdb" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure_save=$1
basestructurename_save=${basestructure_save:r}

#resname=$(grep '^ATOM' $basestructure_save | head -1 | cut -c18-20)

trap "{ echo 'Error detected with code =' $? 2>&1 ; exit 1 }" ZERR

declare -A dihedrals
dihedrals=(c2c3 "C1'-C2'-C3'-C4'"  c1c2 "O4'-C1'-C2'-C3'"  o4c1 "C4'-O4'-C1'-C2'"  c4o4 "C3'-C4'-O4'-C1'")

for d in ${(k)dihedrals}; do
    for state in RNA DNA; do
	case $state in
	    RNA)
		basestructure=$basestructure_save
                # See below. We intentionally compute the overlap (3, 15) so that we can align them easily.
		range=({0..3} {15..35})
		;;
	    DNA)
		basestructure=${basestructure_save:r}.DNA.pdb
	        # DNA is used "between 210 and 330" (Zgarbova 2013)
                # Our chi angle definition is chi + 180, so:
	        # DNA is used "between 30 and 150".
		range=({3..15})
		;;
	esac

	basestructurename=${basestructure:r}
	monomerpdb=$basestructurename.monomer.opt.pdb

	for i in $range; do
	    opttmol=$basestructurename.dihopt$i.done.tmol
	    
	    echo -n "$state $i "
	    python $basedir/dump-dihedral.py $opttmol $monomerpdb ${dihedrals[$d]}
	done
    done > $basestructurename_save.dih.$d.txt
done

python $basedir/plot-dihedral.py $basestructurename_save
