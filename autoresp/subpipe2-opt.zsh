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
$ANTECHAMBER \
    -fi $basestructuretype -i $basestructure \
    -fo gcrt -o $initialgau \
    -ch $opt1check -gk "# HF/6-31G POPT"

python $basedir/add-constraint.py $initialgau $basestructure $baseconstraint > $opt1gau

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
#PBEPBE/Def2TZVP EmpiricalDispersion=GD3 SCRF(CPCM,Solvent=Water)
# Geom=Checkpoint Opt

Optimization phase 2

$CHARGE 1


EOF

python $basedir/add-constraint.py $opt2gauin $basestructure $baseconstraint > $opt2gau

zsh $basedir/g09run.zsh $basename $opt2gau

opt3gau=${basestructurename}.opt3.gau
opt3gauin=${basestructurename}.opt3.gau.in
opt3check=${basename}.opt3.chk

# PBE/Def2TZVPP/COSMO for optimization 2nd phase
sed "/%oldchk=/c%oldchk=$opt2check
/%chk=/c%chk=$opt3check
s/phase 2/phase 3/;s/Def2TZVP/Def2TZVPP/" $opt2gau > $opt3gau
zsh $basedir/g09run.zsh $basename $opt3gau


opt4gau=${basestructurename}.opt4.gau
opt4gauin=${basestructurename}.opt4.gau.in
opt4check=${basename}.opt4.chk

# PBE/6-311++G(3df,3pd)/COSMO for final phase
# Tight optimization & tight integration
# Unfortunately G09 does not support Def2TZVPPD (as of rev. d01)
sed "/%oldchk=/c%oldchk=$opt3check
/%chk=/c%chk=$opt4check
s/phase 3/phase 4/;s/Def2TZVPP/6-311++G(3df,3pd) Int=UltraFine/;s/\\bpopt\\b/popt=VeryTight/" $opt3gau > $opt4gau
zsh $basedir/g09run.zsh $basename $opt4gau


