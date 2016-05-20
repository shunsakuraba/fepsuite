#!/usr/bin/zsh

# Assumes subpipe4-genac is done
if [[ -z $4 ]]; then
    echo "Usage: $0 (base-only structure) (base constraint) (dihedral) (atomtype) (resname)" 1>&2
    echo "Example: $0 inosine.pdb constraint-rna-opt2.txt O4\'-C1\'-N9-C8 CP INO" 1>&2
    echo "The final atom in dihedral angle specification will be given the new atom type"
    exit 1
fi

basedir=${0:h}
basestructure=$1
baseconstraint=$2
dihedral=$3
newatomtype=$4
resname=$5

basestructurename=${basestructure:r}
basename=${basestructure:t:r}

monomerprep=$basestructurename.monomer.prep

if [[ ! -e $monomerprep ]]; then
    echo "Error: monomer prep file is not available" 1>&2
    exit 1
fi

set -x
trap "{ echo 'Error detected with code =' $? 2>&1 ; exit 1 }" ZERR

diheds=($(tr '-' ' ' <<< $dihedral))
repatm=${diheds[4]}

atomtype() {
    targ=$1
    prep=$2
    awk 'substr($0, 7, 4) == sprintf("%-4s", "'$targ'") { a=substr($0, 13, 4); gsub(/ /, "", a); print a }' $prep
}

prevatomtype=$(atomtype $repatm $monomerprep)

# generate modified prep file
dihoptprep=$basestructurename.dihopt.prep
dihoptfrcpre=$basestructurename.dihopt.frcpre
dihoptfrcmod=$basestructurename.dihopt.frcmod
baseparm=$AMBERHOME/dat/leap/parm/parm10.dat

python copydih.py $monomerprep $dihedral $newatomtype $baseparm > $dihoptfrcpre
python replace_prep.py $monomerprep $dihedral $newatomtype > $dihoptprep

if [[ ! -e $basestructurename.extra.frcmod ]]; then
    echo "$basestructurename.extra.frcmod does not exist, making dummy file" 2>&1
    touch $basestructurename.extra.frcmod
fi

# generate amber topology
leap_templ=$basedir/modrna.leap.in
dihoptleap=$basestructurename.dihopt.leap

sed "s/%base%/$basestructurename/g" < $leap_templ > $dihoptleap

TLEAP=$AMBERHOME/bin/tleap

$TLEAP -f $dihoptleap

# dirty hack: revert exchanged atomtype back only in ambtop
# This is necessary to set correct PB radii in pbsa
# Say f%!k to the hard-coded PB radii
newatomtype=$(atomtype $repatm $dihoptprep)
sed -i '/%FLAG AMBER_ATOM_TYPE/,/%FLAG/ {s/'$newatomtype'  /'$prevatomtype'  /}' $basestructurename.dihopt.ambtop


sed -n '/DUMM/,/^$/p' $basestructurename.dihopt.prep | sed '$d' |grep -v 'DUMM' | cut -c 7-10 | tr -d ' ' > $basestructurename.mm.atoms

grep '^ATOM' $basestructure  | cut -c 13-16 | tr -d ' ' > $basestructurename.gaussian.atoms

natom=$(wc -l < $basestructurename.gaussian.atoms)
{
    echo $basestructurename.dihopt.ambtop
    echo $natom
    cat $basestructurename.mm.atoms
    cat $basestructurename.gaussian.atoms
} > $basestructurename.ext.parameter


# XXX: dihedral assignment may be in incorrect order
zsh $basedir/g09fetch.zsh $basename $basename.resp.log $basename.resp.log

template=$basestructurename.dihopt.templ.gau
babel -i g09 $basename.resp.log -o gzmat -xk "%chk=@check@
#PBEPBE/6-311++G(3df,3pd) Int=UltraFine
# EmpiricalDispersion=GD3 SCRF(CPCM,Solvent=Water,Read)
# popt(Z-matrix,maxcycle=9999)" $template

for i in {0..35}; do
    {
	cat $baseconstraint
	echo "dihedral $diheds "$((i * 10))
    } > $basestructurename.cons$i.txt
    optgau=$basestructurename.dihopt.gau$i.gau
    python $basedir/mod-zmatrix.py $template $basestructure $basestructurename.cons$i.txt > $optgau
    echo "QConv=VeryTight\n" >> $optgau
    sed -i "s/@check@/$basename.dih$i.chk/" $optgau
done

