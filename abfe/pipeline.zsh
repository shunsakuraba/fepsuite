#!/bin/zsh

reqstate=$1
stateno=$2
if [[ -z $stateno ]]; then
    echo "This file should be called from template.zsh" 2>&1
    exit 1
fi
shift
shift

TEMPLATE=$ABFE_ROOT/template-pullrestr

#-------- subroutines for production run
mdrun_find_possible_np() {
    least_unit=$1
    shift
    args=($@)
    log_basename="" || true
    dir_basename="."
    #set +e
    is_log=0
    is_prefix=0
    for arg in $args; do
        if [[ $is_log = 1 ]]; then
            log_basename=$arg
        elif [[ $is_prefix = 1 ]]; then
            dir_basename=$arg
            is_prefix=0
        fi
        case $arg in
            -deffnm|-l)
                is_log=1
                ;;
            -multidir)
                is_prefix=1
                ;;
            *)
                is_log=0
                ;;
        esac
    done
    log_basename=$dir_basename/$log_basename
    if [[ -n $PROCS_SAVED ]]; then
        NP=$PROCS_SAVED
    else
        NP=$PROCS
    fi
    while true; do
        echo "Trying with NP=$NP"
        unsetopt ERR_EXIT
        job_mpirun $NP $GMX_MPI mdrun $args
        if [[ $? != 0 ]]; then
            # fail to run. Check log file to see whether it is domain decomposition problem
            if [[ ! -e $log_basename.log ]]; then
                echo "Error: failed during mpirun-mdrun startup. This typically happens when there are a problem with mpirun command."
                exit 1
            fi
            if tail -20 $log_basename.log | grep -q -i "\\(domain\\|prime\\)"; then
                true
            else
                echo "Error: mdrun stopped with errors unrelated to domain size"
                exit 1
            fi 
            if tail -20 $log_basename.log | grep -q -i "There are perturbed non-bonded pair interactions beyond the pair-list cutoff"; then
                echo "Error: mdrun stopped because ligand atom-atom distance exceeded automatically determined rlist value. Try specifying LIGAND_DIAMETER in para_conf.zsh manually."
                exit 1
            fi 
            PREVNP=$NP
            # round up
            (( NP = (NP - 1) / least_unit * least_unit ))
            if (( NP == 0 )) || (( PREVNP == NP )); then
                echo "Error: unable to find proper parallel processes"
                exit 1
            fi
        else
            setopt ERR_EXIT
            break
        fi
        setopt ERR_EXIT
    done
    PROCS_SAVED=$NP
    echo "Normal termination of mdrun"
}

check_rlimit ()
{
    stderrfile=$1
    if grep -q 'This is likely either a 1,4 interaction,' $stderrfile; then
        echo "Error: 1,4- interaction distance exceeds the predefined table distance. Set table-extension= in mdps to larger values. (Must be the possible diameter of the ligand minus initial ligand radius)"
        exit $ERRORCODE
    fi
}

check_replica_probs ()
{
    logfile=$1
    threshold=0.03
    repl_probs=($(grep -A2 "average probabilities:" $logfile | tail -n 1))
    for v in $repl_probs; do
        if [[ $v = "Repl" ]]; then
            continue
        fi
        if (( $v < threshold )); then
            echo "Error: replica exchange probabilities are too low. You need to increase \$NRESTR, \$NCHARGE or \$NANNIH."
            exit $ERRORCODE
        fi
    done
}

set_and_add_rlist ()
{
    MDP=$1

    sync

    integrator=""

    # set rlist, required because we decouple via couple-intramol=no, where exclusions are used within the ligand
    if (( LIGAND_DIAMETER == 0 )); then
        rlist=0.
        while read line; do
            case $line in
            rvdw|rcoulomb)
                d=$(cut -f2 -d '=' <<< $line | tr -d ' \t') || true
                if (( rlist < d )); then
                    (( rlist = d * 1.2 )) # 1.2 for buffering, should be enough for any kernel
                fi
                ;;
            esac
        done < $MDP
        safe_rlist=$(grep safe $ID/diameter.txt | cut -f2 -d ' ')
        if (( safe_rlist > rlist )); then
            (( rlist = safe_rlist ))
        fi
    else
        rlist=$LIGAND_DIAMETER
    fi
    echo "rlist=$rlist" >> $MDP
    echo "verlet-buffer-tolerance=-1" >> $MDP
    if [[ -n $NSTLIST ]] && \
        ! grep -q "integrator[[:space:]]*=[[:space:]]*cg" $MDP && \
        ! grep -q "integrator[[:space:]]*=[[:space:]]*steep" $MDP ; then
        echo "nstlist=$NSTLIST" >> $MDP    
    fi
    
    sync
}

do_run() {
    topol=$1
    conf=$2
    outprefix=$3
    phase=$4
    restr=$5
    ndx=$6
    mdpadd=$7
    cont=$8
    prep=$9
    maxwarn=$10
    if [[ -z $restr ]]; then
        restrcmd=()
    else
        restrcmd=(-r $ID/$restr.pdb)
    fi
    if [[ -n $ndx ]]; then
        index=(-n $ID/$ndx)
    else
        index=()
    fi
    if [[ -n $mdpadd ]]; then
        newmdp=$ID/$conf.mdp
        cat mdp/$phase.mdp $ID/$mdpadd > $newmdp
        mdp=(-f $newmdp)
    else
        mdp=(-f mdp/$phase.mdp)
    fi
    case $cont in
        cont)
            contcmd=(-t $ID/$conf.cpt)
            ;;
        nocont)
            contcmd=()
            ;;
        *)
            echo "do_run $1 $2 $3 $4 $5 $6 $7 $8 $9 $10: \$cont invalid"
            exit $ERRORCODE
            ;;
    esac
    case $prep in
        prep)
            # currently nothing to do
            ;;
        product)
            ;;
        *)
            echo "do_run $1 $2 $3 $4 $5 $6 $7 $8 $9 $10: \$prep invalid"
            exit $ERRORCODE
            ;;
    esac
    job_singlerun $GMX grompp $mdp $index -p $ID/$topol.top -c $ID/$conf.pdb $contcmd -o $ID/$outprefix.$phase.tpr -po $ID/$outprefix.$phase.out.mdp $restrcmd -maxwarn $maxwarn
    stdout=$ID/$outprefix.$phase.stdout
    stderr=$ID/$outprefix.$phase.stderr
    mdrun_find_possible_np 1 -deffnm $ID/$outprefix.$phase -c $ID/$outprefix.$phase.pdb > >(tee $stdout) 2> >(tee $stderr >&2)
    check_rlimit $stderr
}

do_prep_runs() {
    topol=$1
    conf=$2
    output=$3
    initialrestr=$4
    ndx=$5
    pullmdp=$6
    EXPECTED_WARN=0
    case $conf in 
        ligand-*)
        # expect warning from -DPOSRES. This is hacky but other solutions are equally dirty...
        EXPECTED_WARN=1
        ;;
    esac
    do_run $topol $conf $output steep "$initialrestr" "$ndx" "$pullmdp" nocont prep $((EXPECTED_WARN + 1)) # maxwarn 1 for switching
    do_run $topol ${output}.steep $output cg "$initialrestr" "$ndx" "$pullmdp" nocont prep $((EXPECTED_WARN + 1)) # 1 for switching
    do_run $topol ${output}.cg $output nvt "$initialrestr" "$ndx" "$pullmdp" nocont prep $((EXPECTED_WARN + 0)) # 
    do_run $topol ${output}.nvt $output npt "$initialrestr" "$ndx" "$pullmdp" cont prep $((EXPECTED_WARN + 1)) # 1 for gen-vel
}

do_product_runs() {
    topol=$1
    phase=$2
    prev=$3
    output=$4
    ndx=$5
    pullmdp=$6
    nrepl=$7
    maxwarn=$8
    enable_traj_only=$9
    mkdir -p $ID/mdp_addenda || true
    mkdir -p $ID/mdp || true
    mdprestr=()
    if [[ -n $pullmdp ]]; then
        mdprestr=(--additional-mdp $ID/$pullmdp)
    fi
    nprerun=0
    # In both annihilation phases we optimize separately, because solvent and protein are totally different environments.
    # In annihilation-complex phase we restart from optimized lambda values in annihilation-lig
    if [[ $phase = "annihilation-lig" ]] || [[ $phase = "annihilation-complex" ]]; then
        nprerun=$ANNIH_LAMBDA_OPT
    fi
    prevcrd=$ID/$prev
    update=()
    for iprerun in {0..$nprerun}; do
        if (( iprerun < nprerun )); then
            infix=".pre$iprerun"
        else
            infix=""
        fi
        if (( iprerun != 0 )); then
            update=(--update $ID/$phase.0/$phase.pre$((iprerun-1)).log --update-nth $iprerun)
        elif (( iprerun == 0 )) && [[ $phase = "annihilation-complex" ]]; then
            # use annihilation-lig's final prerun output to update the lambda values
            update=(--update $ID/annihilation-lig.0/annihilation-lig.pre$(($ANNIH_LAMBDA_OPT-1)).log --update-nth $ANNIH_LAMBDA_OPT)
        fi
        $PYTHON3 $ABFE_ROOT/generate_decoupling.py --template-dir $TEMPLATE $mdprestr --output-mdp $ID/mdp_addenda --mode $phase -N $nrepl --mol-name $LIG_GMX $update
        DIRS=()
        for i in {0..$((nrepl-1))}; do
            CODE=$phase.$i
            RUNDIR=$ID/$CODE
            if (( iprerun != 0 )); then
                prevcrd=$RUNDIR/$phase.pre$((iprerun-1))
            fi
            mkdir -p $RUNDIR || true
            MDP=$ID/mdp/$phase-$i$infix.mdp
            postprocess=(cat)
            if [[ -n $enable_traj_only ]] && (( i != $enable_traj_only )); then
                postprocess=(sed "1,/LRCONLY_BEGIN/{/compressed/d};/LRCONLY_BEGIN/,/LRCONLY_END/d")
            elif [[ -n $enable_traj_only ]]; then
                postprocess=(sed "1,/LRCONLY_BEGIN/{/compressed/d}")
            fi
            cat ./mdp/run.mdp $ID/mdp_addenda/${phase}-$i.mdp | $postprocess  > $MDP

            # modifying rlist is needed only for charging and (due to GROMACS internal setup) annihilation
            case $phase in 
                charging-complex|charging-lig|annihilation-complex|annihilation-lig)
                set_and_add_rlist $MDP
                ;;
            esac

            if (( iprerun < nprerun )); then
                dt=$(grep "^\\s*dt\\s*=" $MDP | cut -d '=' -f2 | cut -d ';' -f1)
                declare -i nsteps # fix to integer
                (( nsteps = ANNIH_LAMBDA_OPT_LENGTH / dt ))
                sed -i "/nsteps/c nsteps = $nsteps" $MDP # in the future thie may be sent to postprocess
            fi
            ndx_grompp=()
            if [[ -n $ndx ]]; then
                ndx_grompp=(-n $ID/$ndx)
            fi
            job_singlerun $GMX grompp -f $MDP -p $ID/$topol -c $prevcrd.pdb -t $prevcrd.cpt -o $RUNDIR/$phase$infix.tpr -po $ID/$output.$phase$infix.$i.mdp $ndx_grompp -maxwarn $maxwarn
            DIRS+=$RUNDIR
        done

        stdout=$ID/$phase.stdout
        stderr=$ID/$phase.stderr
        mdrun_find_possible_np $nrepl -multidir $DIRS -deffnm $phase$infix -c $phase$infix.pdb -replex 500 > >(tee $stdout) 2> >(tee $stderr >&2)
        if (( iprerun == nprerun )); then
            check_replica_probs $ID/$phase.0/$phase.log
        fi
    done
}

do_eval_run() {
    topol=$1
    phase=$2
    prev=$3
    output=$4
    ndx=$5
    maxwarn=$6

    { grep -v "dispcorr" ./mdp/run.mdp | grep -v 'nstxtcout' ; cat $TEMPLATE/$phase.mdp } | sed "s/{group_mol}/$LIG_GMX/g" > $ID/$phase.mdp
    ndx_grompp=()
    if [[ -n $ndx ]]; then
        ndx_grompp=(-n $ID/$ndx)
    fi
    case $phase in 
        charging-complex|charging-lig|annihilation-complex|annihilation-lig|lr-*)
        set_and_add_rlist $ID/$phase.mdp
    ;;
    esac
    
    job_singlerun $GMX grompp -f $ID/$phase.mdp -p $ID/$topol -c $ID/$prev.pdb -t $ID/$prev.cpt -o $ID/$phase.tpr -po $ID/$output.$phase.mdp $ndx_grompp -maxwarn $maxwarn
    mdrun_find_possible_np 1 -deffnm $ID/$phase -rerun $ID/$prev.xtc
}

do_bar() {
    mode=$1
    nrepl=$2
    job_singlerun $GMX bar -b $RUN_PROD -f $ID/$mode.{0..$((nrepl-1))}/$mode.xvg -o $ID/$mode.bar_diff.xvg > $ID/$mode.bar.log
}

do_exp() {
    lr_output=$1
    non_lr_output=$2
    non_lr_repl=$3
    TEMP=$4
    $PYTHON3 $ABFE_ROOT/lr_exp.py --long $ID/$lr_output.edr --short $ID/$non_lr_output.$non_lr_repl/$non_lr_output.edr --temp $TEMP --time-begin $RUN_PROD --output $ID/$lr_output.lrc.txt
}

do_solvate() {
    solute=$1
    topol=$2
    output_trunk=$3
    topol_sol=${topol%.top}-sol.top
    solute_sol=${solute%.*}-sol.pdb

    cp $topol $topol_sol
    job_singlerun $GMX solvate -cp $solute -p $topol_sol -cs $WATER_STRUCTURE -o $solute_sol
    touch ${output_trunk}-dummy.mdp
    job_singlerun $GMX grompp -f ${output_trunk}-dummy.mdp -p $topol_sol -c $solute_sol -po ${output_trunk}-dummy-out.mdp -o ${output_trunk}-sol.tpr
    topol_ion=${output_trunk}-ion.top
    cp $topol_sol $topol_ion
    echo SOL | job_singlerun $GMX genion -s ${output_trunk}-sol.tpr -o ${output_trunk}-ion.pdb -p $topol_ion -pname $ION_POSITIVE -nname $ION_NEGATIVE -conc $IONIC_STRENGTH -neutral
}

generate_ndx() {
    conf=$1
    pptop=$2
    output=$3
    is_complex_or_ligand=$4

    receptor=()
    if [[ $is_complex_or_ligand = complex ]]; then
        receptor=(--receptor $RECEPTOR_MDTRAJ)
    fi

    if [[ -e $conf ]]; then
        ; # do nothing
    elif [[ -e $conf.pdb ]]; then
        conf=$conf.pdb
    elif [[ -e $conf.gro ]]; then
        conf=$conf.gro
    else
        echo "Conformation not found (expected $conf.{pdb or gro}" 1>&2
        return 1
    fi
    $PYTHON3 $ABFE_ROOT/make_ndx.py --structure $conf --topology $pptop --ligand $LIG_GMX $receptor --output $output
}

run_charge_correction() {
    local top=$ID/pp_run.top
    local reftpr=$ID/prerun.run.final.tpr
    local complex_traj=$ID/prerun.run.recpbc.xtc
    local ndx=$ID/complex.ndx

    APBS=${APBS:-apbs}

    # Gromacs has problems shifting back the ligand/receptor to the center, so we reuse the structure generated at the pipeline 3
    local simlen
    simlen=$(job_singlerun $GMX check -f $complex_traj 2>&1 >/dev/null | grep '^Time' | awk '{ print ($2-1) * $3 }')
    local tt
    for i in {1..$CHARGE_CORRECTION_NSAMP}; do
        # first RUN_PROD ps removed
        (( tt = RUN_PROD + simlen * (i + 0.0) / CHARGE_CORRECTION_NSAMP ))
        charge_temp=$ID/charge_corr_in.pdb
        rm -f $charge_temp
        echo "System" | job_singlerun $GMX trjconv -s $reftpr -f $complex_traj -o $charge_temp -n $ndx -dump $tt
        pushd $ID # due to intermediate files it is better executed at each $ID
        #$PYTHON3 $ABFE_ROOT/charge_correction_generator.py --top ${top##$ID/} --pdb ${charge_temp##$ID/} --ndx ${ndx##$ID/} --summary charge_result$i.txt --ligand-group Ligand --receptor-group Receptor > charge_detail$i.txt
        popd
    done
    cat $ID/charge_detail{1..$CHARGE_CORRECTION_NSAMP}.txt > $ID/charge_correction.txt

    top=$ID/ligand-ion.top
    reftpr=$ID/charging-lig.0/charging-lig.tpr
    local lig_traj=$ID/charging-lig.0/charging-lig.xtc
    ndx=$ID/ligand.ndx
    simlen=$(job_singlerun $GMX check -f $lig_traj 2>&1 >/dev/null | grep '^Time' | awk '{ print ($2-1) * $3 }')
    for i in {1..$CHARGE_CORRECTION_NSAMP}; do
        # skip first RUN_PROD ps
        (( tt = RUN_PROD + (simlen - RUN_PROD) * (i + 0.0) / CHARGE_CORRECTION_NSAMP ))
        charge_temp=$ID/charge_corr_in.pdb
        rm -f $charge_temp
        echo "Ligand\nSystem" | job_singlerun $GMX trjconv -s $reftpr -f $lig_traj -o $charge_temp -n $ndx -dump $tt -pbc mol -center -ur compact
        pushd $ID # due to intermediate files it is better executed at each $ID
        #$PYTHON3 $ABFE_ROOT/charge_correction_generator.py --top ${top##$ID/} --pdb ${charge_temp##$ID/} --ndx ${ndx##$ID/} --summary charge_result_lig$i.txt --ligand-group Ligand > charge_detaillig$i.txt
        popd
    done
    cat $ID/charge_detaillig{1..$CHARGE_CORRECTION_NSAMP}.txt > $ID/charge_correction_lig.txt
}

charge_correction() {
    local chargefile=$1

    # FIXME: do we need ligand structures as well?
    local totcharge=$(cat $chargefile)
    local abstotcharge=$totcharge
    if (( abstotcharge < 0 )); then
        (( abstotcharge = - abstotcharge ))
    fi
    if (( abstotcharge < 1e-4 )); then
        echo "No charge correction needed"
        rm -f $ID/charge_correction.txt
        rm -f $ID/charge_correction_lig.txt
    else
        run_charge_correction
    fi
}

# ---- end of subroutines for production run

main() {
    case $reqstate,$stateno in
        query,all)
            echo {1..12}
            ;;
        query,1)
            echo "DEPENDS=(); PPM=\$COMPLEX_PARA ; MULTI=1"
            ;;
        run,1)
            # Check version. There are too many issues...
            version=$(job_singlerun $GMX --version | grep 'GROMACS version:' | tr -s ' ' | cut -f3 -d' ')
            case $version in
                202[01]*|201[6789]*|2022|2022.[1-4])
                    echo "GROMACS version $version is not supported because it cannot detect error cases with excluded interactions properly"
                    exit 1
                ;;
                2022.[56789]*)
                    echo "GROMACS version $version (supported version)"
                ;;
                *)
                    echo "GROMACS version $version. We are not sure the program works on this version; we will continue anyway."
                ;;
            esac
            touch mdp/dummy.mdp
            # TODO: check whether maxwarn 0 is sufficient
            job_singlerun $GMX grompp -f mdp/dummy.mdp -p $ID/topol_ionized.top -c $ID/conf_ionized -o $ID/pp.tpr -po $ID/pp.mdp -maxwarn 0 -pp $ID/pp.top
            if [[ ! -e mdp/dummy_flex.mdp ]]; then
                echo "define = -DFLEXIBLE" > mdp/dummy_flex.mdp
            fi
            job_singlerun $GMX grompp -f mdp/dummy_flex.mdp -p $ID/topol_ionized.top -c $ID/conf_ionized -o $ID/pp_flex.tpr -po $ID/pp_flex.mdp -maxwarn 0 -pp $ID/pp_flex.top
            generate_ndx $ID/conf_ionized $ID/pp.top $ID/complex.ndx complex
            do_prep_runs topol_ionized conf_ionized prep conf_ionized "" ""
            ;;
        query,2)
            echo "DEPENDS=(1); PPM=\$COMPLEX_PARA ; MULTI=1"
            ;;
        run,2)
            do_run topol_ionized prep.npt prerun run "" "" "" cont prep 0
            ;;
        query,3)
            echo "DEPENDS=(2); PPM=1; MULTI=1"
            ;;
        run,3)
            # Compute RMSD of ligands for thresholding and later uses.
            # This seeming strange procedure was set to support cases which do not work with gmx trjconv -pbc cluster.
            # Since this is time-evolved simulation without replica exchange, we can assume -pbc nojump works (assumes no >0.5pbc jump during prep.min/prep.nvt/prep.npt).
            echo "Ligand+Receptor\nSystem" | job_singlerun $GMX trjconv -s $ID/pp_flex.tpr -f $ID/prerun.run.xtc -o $ID/prerun.run.nojump.xtc -b $RUN_PROD -center -n $ID/complex.ndx -pbc nojump
            # Then wrap the system around using -pbc mol and -ur compact. This and above line prevents "split-ligand" and "split-receptor" artifacts.
            echo "System" | job_singlerun $GMX trjconv -s $ID/pp_flex.tpr -f $ID/prerun.run.nojump.xtc -o $ID/prerun.run.recpbc.xtc -b $RUN_PROD -pbc mol -n $ID/complex.ndx -ur compact
            # Get the final simulation time. Note that prerun.run.recpbc.xtc has been chopped off initial RUN_PROD ps, so using prerun.run.xtc instead.
            # gmx check outputs into the standard error
            LAST=$(job_singlerun $GMX check -f $ID/prerun.run.xtc 2>&1 >/dev/null | grep '^Time' | awk '{ print ($2-1) * $3 }')
            # Convert the final structure into PBC-fixed one
            echo "System" | job_singlerun $GMX trjconv -s $ID/pp_flex.tpr -f $ID/prerun.run.recpbc.xtc -o $ID/prerun.run.final.pdb -n $ID/complex.ndx -dump $LAST
            # Then calculate the rms
            # remake the tpr file because atoms may not have defined mass
            job_singlerun $GMX grompp -f mdp/dummy.mdp -p $ID/topol_ionized.top -c $ID/prerun.run.final.pdb -o $ID/prerun.run.final.tpr -po $ID/prerun_dummy.mdp
            echo "Receptor\nLigand" | job_singlerun $GMX rms -s $ID/prerun.run.final.tpr -f $ID/prerun.run.recpbc.xtc -o $ID/prerun.rms.fromfinal.xvg -n $ID/complex.ndx
            $PYTHON3 $ABFE_ROOT/rms_check.py --rms=$ID/prerun.rms.fromfinal.xvg --threshold=$EQ_RMSD_CUTOFF || { echo "RMS of ligands too large, aborting the calculation" 1>&2; false }
            $PYTHON3 $ABFE_ROOT/find_restr_from_md.py --lig-sel "Ligand" --prot-sel "Receptor" --index $ID/complex.ndx --topology $ID/conf_ionized.pdb --trajectory $ID/prerun.run.recpbc.xtc --output $ID/restrinfo
            # Restraint for annihilation and charging
            $PYTHON3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --mdp $ID/restr_pull.mdp --ndx $ID/restr_pull.ndx
            # Restraint decoupling, used in restraint phase
            $PYTHON3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --decouple-B --mdp $ID/restr_pull_decouple.mdp --ndx $ID/restr_pull_decouple.ndx
            diff -q $ID/restr_pull.ndx $ID/restr_pull_decouple.ndx || { echo "Two pull indices are not identical" 1>&2; false }

            # generate topologies and conformations:
            # stage 4, 5: ligand discharge calculation \
            # stage 6: ligand annihilation calculation +-> for these three, generate ligand only structure and topologies
            # stage 7: LRC ligand                      / 
            # stage 8: complex discharge calculation     \
            # stage 9: complex annihilation calculation  +-> for these four, complex structures as well as topologies are used
            # stage 10: complex restraint calculation    |
            # stage 11: LRC                              /
            
            # prepare pp because we calculate ligand's total charge here
            job_singlerun $GMX grompp -f mdp/run.mdp -p $ID/topol_ionized.top -c $ID/prerun.run.pdb -t $ID/prerun.run.cpt -o $ID/pp.tpr -po $ID/pp.mdp $restrcmd -maxwarn 0 -pp $ID/pp_run.top
            $PYTHON3 $ABFE_ROOT/generate_ligand_topology.py --mol $LIG_GMX --topology $ID/pp_run.top --structure $ID/prerun.run.final.pdb --index $ID/complex.ndx --output-ligand-structure $ID/ligand.pdb --output-ligand-topology $ID/ligand.top --total-charge $ID/totalcharge.txt

            # run ligand-only system to determine the cutoff distance, maxwarn = 1 for PME charge
            job_singlerun $GMX editconf -f $ID/ligand.pdb -d 3.0 -o $ID/ligand-bigbox.pdb
            job_singlerun $GMX grompp -f mdp/ligsample.mdp -c $ID/ligand-bigbox.pdb -p $ID/ligand.top -o $ID/ligand_only.tpr -maxwarn 1
            # run directly without do_run because I don't want to modify NSTLIST and SAVE_*
            job_mpirun 1 $GMX_MPI mdrun -deffnm $ID/ligand_only
            # get the diameter
            $PYTHON3 $ABFE_ROOT/ligand_diameter.py --traj=$ID/ligand_only.xtc --structure=$ID/ligand.pdb > $ID/diameter.txt

            # solvate ligand-only confs
            real_water_thickness=$WATER_THICKNESS
            diameter_safe=$(grep safe $ID/diameter.txt | cut -f2 -d ' ')
            if (( real_water_thickness < diameter_safe )); then
                (( real_water_thickness = diameter_safe + 0.1 ))
            fi
            job_singlerun $GMX editconf -f $ID/ligand.pdb -d $real_water_thickness -bt dodecahedron -o $ID/ligand-box.pdb
            do_solvate $ID/ligand-box.pdb $ID/ligand.top $ID/ligand
            generate_ndx $ID/ligand-ion $ID/ligand-ion.top $ID/ligand.ndx ligand

            # generate index files for each topology
            # complex: index can be reused
            cat $ID/complex.ndx $ID/restr_pull.ndx > $ID/complex_pull.ndx
            # ligand: not needed because never pulled
        
            # generate mdp files for input
            for d in $TEMPLATE/*-{lig,complex}.mdp; do
                cp $d $ID/
            done
            $PYTHON3 $ABFE_ROOT/resurrect_flexible.py --flexible $ID/pp_flex.top --top $ID/ligand-ion.top --output $ID/ligand-ion-flex.top
            cat $ID/complex.ndx $ID/restr_pull.ndx > $ID/complex_with_pull.ndx
            ;;
        query,4)
            echo "DEPENDS=(3); PPM=\$LIG_PARA; MULTI=1"
            ;;
        run,4)
            do_prep_runs ligand-ion-flex ligand-ion charging-lig ligand-ion ligand ""
            ;;
        query,5)
            echo "DEPENDS=(4); PPM=\$LIG_PARA; MULTI=\$NCHARGE"
            ;;
        run,5)
            # product run for charge-discharge ligand
            do_product_runs ligand-ion charging-lig charging-lig.npt charging-lig "" "" $NCHARGE 0 0
            ;;
        query,6)
            echo "DEPENDS=(5); PPM=\$LIG_PARA; MULTI=\$NANNIH"
            ;;
        run,6)
            # product run for ligand annihilation, starting from totally discharged conformation
            do_product_runs ligand-ion annihilation-lig charging-lig.$((NCHARGE - 1))/charging-lig annihilation-lig "" "" $NANNIH 0 $(( NANNIH - 1 ))
            ;;
        query,7)
            echo "DEPENDS=(5 6); PPM=\$LIG_PARA; MULTI=1"
            ;;
        run,7)
            # re-eval for LRC for charging #0
            do_eval_run ligand-ion lr-lig charging-lig.0/charging-lig lr-lig ligand 0
            # re-eval for LRC for annihilation last
            do_eval_run ligand-ion lr-annihilation-lig annihilation-lig.$((NANNIH-1))/annihilation-lig lr-annihilation-lig ligand $((NANNIH-1))
            ;;
        query,8)
            echo "DEPENDS=(3); PPM=\$COMPLEX_PARA; MULTI=\$NCHARGE"
            ;;
        run,8)
            # product run for complex system (discharging)
            do_product_runs topol_ionized charging-complex prerun.run charging-complex complex_with_pull restr_pull.mdp $NCHARGE 0
            ;;
        query,9)
            echo "DEPENDS=(6 8); PPM=\$COMPLEX_PARA; MULTI=\$NANNIH"
            ;;
        run,9)
            # product run for complex system (annihilation)
            do_product_runs topol_ionized annihilation-complex charging-complex.$((NCHARGE-1))/charging-complex annihilation-complex complex_with_pull restr_pull.mdp $NANNIH 0 $((NANNIH - 1))
            ;;
        query,10)
            echo "DEPENDS=(3); PPM=\$COMPLEX_PARA; MULTI=\$NRESTR"
            ;;
        run,10)
            # product run for restraint decopuling (100% restraint -> 0% restraint) FEP
            do_product_runs topol_ionized restraint prerun.run restraint complex_with_pull restr_pull_decouple.mdp $NRESTR 1 $((NRESTR - 1))
            ;;
        query,11)
            echo "DEPENDS=(9 10); PPM=\$COMPLEX_PARA; MULTI=1"
            ;;
        run,11)
            # eval run (LRC) for restraint final (unrestrained)
            do_eval_run topol_ionized lr-complex restraint.$((NRESTR - 1))/restraint lr-complex complex 0
            # re-eval for LRC annihilation final
            do_eval_run topol_ionized lr-annihilation-complex annihilation-complex.$((NANNIH - 1))/annihilation-complex lr-annihilation-complex complex $((NANNIH - 1))
            ;;
        query,12)
            echo "DEPENDS=(3 4 5 6 7 8 9 10 11); PPM=1; MULTI=1; CPU_ONLY_STAGE=yes"
            ;;
        run,12)
            TEMP=$(grep ref_t mdp/run.mdp | cut -d '=' -f2)
            # do gmx bar
            do_bar charging-lig $NCHARGE
            do_bar charging-complex $NCHARGE
            do_bar annihilation-lig $NANNIH
            do_bar annihilation-complex $NANNIH
            do_bar restraint $NRESTR
            # LRC
            do_exp lr-lig charging-lig 0 $TEMP
            do_exp lr-annihilation-lig annihilation-lig $((NANNIH - 1)) $TEMP
            do_exp lr-complex restraint $((NRESTR - 1)) $TEMP
            do_exp lr-annihilation-complex annihilation-complex $((NANNIH - 1)) $TEMP
            # Charge correction
            charge_correction $ID/totalcharge.txt
            # Sum evertyhing up
            $PYTHON3 $ABFE_ROOT/calc_bar_replex.py --basedir $ID --restrinfo $ID/restrinfo --temp $TEMP > $ID/result.txt
            ;;
        *)
            echo "Unexpected stage outputs"
            exit 1
            ;;
    esac
}

main

