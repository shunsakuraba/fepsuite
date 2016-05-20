#!/usr/bin/zsh

# Constraint optimization.
if [[ -z $2 ]]; then
    echo "Usage: $0 (structure) (constraint file)" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure=$1
baseconstraint=$2

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

basestructuretype=${basestructure:e}
if [[ $basestructuretype = "gau" ]]; then
    echo "Base structure file name must not be gau" 1>&2
    exit 1
fi
basestructurename=${basestructure:r}
basename=${basestructure:t:r}

initialgau=${basestructurename}.init.gau

opt1gau=${basestructurename}.opt1.gau
opt1check=${basename}.opt1.chk

# Use babel, not antechamber (antechamber may swap atoms and emits no information for recovering atom names!)
$OPENBABEL $basestructure -ogzmat -xk "%chk=$opt1check
# HF/6-31G POPT(Z-matrix,maxcycle=5000)" $initialgau

# Fix 'Holmium' problem and broken multiplicity
sed -i "s/Ho/H/;/^0 /c $CHARGE 1" $initialgau

python $basedir/mod-zmatrix.py $initialgau $basestructure $baseconstraint > $opt1gau

# Optimize the structure with constraint
# HF/6-31G level first to remove large crashes
zsh $basedir/g09run.zsh $basename $opt1gau

opt2gau=${basestructurename}.opt2.gau
opt2gauin=${basestructurename}.opt2.gau.in
opt2check=${basename}.opt2.chk

# PBE/Def2TZVP/COSMO for optimization
cat > $opt2gauin << EOF
--Link1--
%oldchk=$opt1check
%chk=$opt2check
#PBEPBE/6-31G* EmpiricalDispersion=GD3 SCRF(CPCM,Solvent=Water)
# Geom=Checkpoint popt(Z-matrix,maxcycle=5000)

Optimization phase 2

$CHARGE 1


EOF

zsh $basedir/g09run.zsh $basename $opt2gau

opt3gau=${basestructurename}.opt3.gau
opt3gauin=${basestructurename}.opt3.gau.in
opt3check=${basename}.opt3.chk

# PBE/Def2TZVPP/COSMO for optimization 2nd phase
sed "/%oldchk=/c%oldchk=$opt2check
/%chk=/c%chk=$opt3check
s/phase 2/phase 3/;s/6-31G\*/Def2TZVPP Int=UltraFine/" $opt2gau > $opt3gau
zsh $basedir/g09run.zsh $basename $opt3gau


opt4gau=${basestructurename}.opt4.gau
opt4gauin=${basestructurename}.opt4.gau.in
opt4check=${basename}.opt4.chk

# PBE/6-311++G(3df,3pd)/COSMO for final phase
# Tight optimization & tight integration
# Unfortunately G09 does not support Def2TZVPPD (as of rev. d01)
# thus we use LP instead.
sed "/%oldchk=/c%oldchk=$opt3check
/%chk=/c%chk=$opt4check
s/phase 3/phase 4/;s/Def2TZVPP/6-311++G(3df,3pd)/;s/SCRF(/SCRF(Read,/" $opt3gau > $opt4gau
echo "QConv=VeryTight\n" >> $opt4gau
zsh $basedir/g09run.zsh $basename $opt4gau


# work around Gaussian bug (?)
logfile=${opt4gau:r}.log

while true; do
    zsh $basedir/g09fetch.zsh $basename $logfile ${logfile:t}
    has_error=n
    tail -4 ${opt4gau:r}.log | grep -q 'link 9999' && has_error=y || true
    if [[ $has_error = n ]]; then
	break
    fi
    opt4gauretry=${basestructurename}.opt4r1.gau
    sed '/oldchk/d;s/popt(/popt(restart,/' $opt4gau > $opt4gauretry
    zsh $basedir/g09run.zsh $basename $opt4gauretry
    logfile=${opt4gauretry:r}.log
done
