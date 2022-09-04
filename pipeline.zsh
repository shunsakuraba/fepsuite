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
    #set +e
    is_log=0
    for arg in $args; do
        if [[ $is_log = 1 ]]; then
            log_basename=$arg
        fi
        case $arg in
            -deffnm|-l)
                is_log=1
                ;;
            *)
                is_log=0
                ;;
        esac
    done
    if [[ -n $PROCS_SAVED ]]; then
        NP=$PROCS_SAVED
    else
        NP=$PROCS
    fi
    while true; do
        echo "Trying with NP=$NP"
        mpirun_ $NP $GMX_MPI mdrun $args $NSTLIST_CMD
        if [[ $? != 0 ]]; then
            # fail to run. Check log file to see whether it is domain decomposition problem
            tail -20 $log_basename.log | grep -q -i "\\(domain\\|prime\\)" || { echo "Error: domain-unrelated error"; exit $ERRORCODE }
            PREVNP=$NP
            # round up
            (( NP = (NP - 1) / least_unit * least_unit ))
            if (( NP == 0 )) || (( PREVNP == NP )); then
                echo "Error: unable to find proper parallel processes"
                exit $ERRORCODE
            fi
        else
            break
        fi
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
            echo "Error: replica exchange probabilities are too low. You need to increase \$NCHARGE or \$NANNIH."
            exit $ERRORCODE
        fi
    done
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
    nstlist_cmd=()
    if [[ -n $NSTLIST ]] && [[ $phase != steep ]]; then
        nstlist_cmd=(-nstlist $NSTLIST)
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
    $SINGLERUN $GMX grompp $mdp $index -p $ID/$topol.top -c $ID/$conf.pdb $contcmd -o $ID/$outprefix.$phase.tpr -po $ID/$outprefix.$phase.out.mdp $restrcmd -maxwarn $maxwarn
    stdout=$ID/$outprefix.$phase.stdout
    stderr=$ID/$outprefix.$phase.stderr
    mdrun_find_possible_np 1 -deffnm $ID/$outprefix.$phase -c $ID/$outprefix.$phase.pdb $nstlist_cmd > >(tee $stdout) 2> >(tee $stderr >&2)
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
    if [[ $phase = "annihilation-lig" ]]; then
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
            if (( iprerun < nprerun )); then
                sync $MDP
                dt=$(grep "^\\s*dt\\s*=" $MDP | cut -d '=' -f2 | cut -d ';' -f1)
                declare -i nsteps # fix to integer
                (( nsteps = ANNIH_LAMBDA_OPT_LENGTH / dt ))
                sed -i "/nsteps/c nsteps = $nsteps" $MDP # in the future thie may be sent to postprocess
                sync $MDP
            fi
            ndx_grompp=()
            if [[ -n $ndx ]]; then
                ndx_grompp=(-n $ID/$ndx)
            fi
            $SINGLERUN $GMX grompp -f $MDP -p $ID/$topol -c $prevcrd.pdb -t $prevcrd.cpt -o $RUNDIR/$phase$infix.tpr -po $ID/$output.$phase$infix.$i.mdp $ndx_grompp -maxwarn $maxwarn
            DIRS+=$RUNDIR
        done
        nstlist_cmd=()
        if [[ -n $NSTLIST ]]; then
            nstlist_cmd=(-nstlist $NSTLIST)
        fi
        stdout=$ID/$phase.stdout
        stderr=$ID/$phase.stderr
        mdrun_find_possible_np $nrepl -multidir $DIRS -deffnm $phase$infix -c $phase$infix.pdb -replex 500 $nstlist_cmd > >(tee $stdout) 2> >(tee $stderr >&2)
        if (( iprerun = nprerun )); then
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
    $SINGLERUN $GMX grompp -f $ID/$phase.mdp -p $ID/$topol -c $ID/$prev.pdb -t $ID/$prev.cpt -o $ID/$phase.tpr -po $ID/$output.$phase.mdp $ndx_grompp -maxwarn $maxwarn
    mdrun_find_possible_np 1 -deffnm $ID/$phase -rerun $ID/$prev.xtc $nstlist_cmd
}

do_bar() {
    mode=$1
    nrepl=$2
    $SINGLERUN $GMX bar -b $RUN_PROD -f $ID/$mode.{0..$((nrepl-1))}/$mode.xvg -o $ID/$mode.bar_diff.xvg > $ID/$mode.bar.log
}

do_exp() {
    lr_output=$1
    non_lr_output=$2
    non_lr_repl=$3
    TEMP=$4
    python3 $ABFE_ROOT/lr_exp.py --long $ID/$lr_output.edr --short $ID/$non_lr_output.$non_lr_repl/$non_lr_output.edr --temp $TEMP --time-begin $RUN_PROD --output $ID/$lr_output.lrc.txt
}

do_solvate() {
    solute=$1
    topol=$2
    output_trunk=$3
    topol_sol=${topol%.top}-sol.top
    solute_sol=${solute%.*}-sol.pdb

    cp $topol $topol_sol
    $SINGLERUN $GMX solvate -cp $solute -p $topol_sol -cs $WATER_STRUCTURE -o $solute_sol
    touch ${output_trunk}-dummy.mdp
    $SINGLERUN $GMX grompp -f ${output_trunk}-dummy.mdp -p $topol_sol -c $solute_sol -po ${output_trunk}-dummy-out.mdp -o ${output_trunk}-sol.tpr
    topol_ion=${output_trunk}-ion.top
    cp $topol_sol $topol_ion
    echo SOL | $SINGLERUN $GMX genion -s ${output_trunk}-sol.tpr -o ${output_trunk}-ion.pdb -p $topol_ion -pname $ION_POSITIVE -nname $ION_NEGATIVE -conc $IONIC_STRENGTH -neutral
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
    python3 $ABFE_ROOT/make_ndx.py --structure $conf --topology $pptop --ligand $LIG_GMX $receptor --output $output
}

charge_correction() {
    CHARGEFILE=$1
    REFSTRUCTURE=$2

    # FIXME: do we need ligand structures as well?
    TOTCHARGE=$(cat $CHARGEFILE)
    ABSTOTCHARGE=$TOTCHARGE
    if (( ABSTOTCHARGE < 0 )); then
        (( ABSTOTCHARGE = - ABSTOTCHARGE ))
    fi
    if (( ABSTOTCHARGE < 1e-4 )); then
        echo "No charge correction needed"
        echo 0.0 > $ID/charge_correction.txt
    else
        echo "Error: charge correction is not implemented yet."
        exit 1
    fi
}

# ---- end of subroutines for production run

main() {
    if [[ $reqstate = run ]]; then
        set -x
    fi
    case $reqstate,$stateno in
        query,all)
            echo {1..12}
            ;;
        query,1)
            echo "DEPENDS=(); (( PROCS = COMPLEX_PARA ))"
            ;;
        run,1)
            touch mdp/dummy.mdp
            # TODO: check whether maxwarn 0 is sufficient
            $SINGLERUN $GMX grompp -f mdp/dummy.mdp -p $ID/topol_ionized.top -c $ID/conf_ionized -o $ID/pp.tpr -po $ID/pp.mdp -maxwarn 0 -pp $ID/pp.top
            if [[ ! -e mdp/dummy_flex.mdp ]]; then
                echo "define = -DFLEXIBLE" > mdp/dummy_flex.mdp
            fi
            $SINGLERUN $GMX grompp -f mdp/dummy_flex.mdp -p $ID/topol_ionized.top -c $ID/conf_ionized -o $ID/pp_flex.tpr -po $ID/pp_flex.mdp -maxwarn 0 -pp $ID/pp_flex.top
            generate_ndx $ID/conf_ionized $ID/pp.top $ID/complex.ndx complex
            do_prep_runs topol_ionized conf_ionized prep conf_ionized "" ""
            ;;
        query,2)
            echo "DEPENDS=(1); (( PROCS = COMPLEX_PARA ))"
            ;;
        run,2)
            do_run topol_ionized prep.npt prerun run "" "" "" cont prep 0
            ;;
        query,3)
            echo "DEPENDS=(2);PROCS=1;NOGPU_STAGE=yes"
            ;;
        run,3)
            # Compute RMSD of ligands for thresholding
            echo "Receptor\nSystem" | $SINGLERUN $GMX trjconv -s $ID/prerun.run.tpr -f $ID/prerun.run.xtc -o $ID/prerun.run.recpbc.xtc -b $RUN_PROD -center -pbc mol -n $ID/complex.ndx -ur compact
            echo "Receptor\nSystem" | $SINGLERUN $GMX trjconv -s $ID/prerun.run.tpr -f $ID/prerun.run.pdb -o $ID/prerun.run.final.pdb -center -pbc mol -n $ID/complex.ndx -ur compact
            echo "Receptor\nLigand" | $SINGLERUN $GMX rms -s $ID/prerun.run.final.pdb -f $ID/prerun.run.recpbc.xtc -o $ID/prerun.rms.fromfinal.xvg -n $ID/complex.ndx
            python3 $ABFE_ROOT/rms_check.py --rms=$ID/prerun.rms.fromfinal.xvg --threshold=$EQ_RMSD_CUTOFF || { echo "RMS of ligands too large, aborting the calculation" 1>&2; false }
            python3 $ABFE_ROOT/find_restr_from_md.py --lig-sel "Ligand" --prot-sel "Receptor" --index $ID/complex.ndx --topology $ID/conf_ionized.pdb --trajectory $ID/prerun.run.recpbc.xtc --output $ID/restrinfo
            # Restraint for annihilation and charging
            python3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --mdp $ID/restr_pull.mdp --ndx $ID/restr_pull.ndx
            # Restraint decoupling, used in restraint phase
            python3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --decouple-B --mdp $ID/restr_pull_decouple.mdp --ndx $ID/restr_pull_decouple.ndx
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
            $SINGLERUN $GMX grompp -f mdp/run.mdp -p $ID/topol_ionized.top -c $ID/prerun.run.pdb -t $ID/prerun.run.cpt -o $ID/pp.tpr -po $ID/pp.mdp $restrcmd -maxwarn 0 -pp $ID/pp_run.top
            python3 $ABFE_ROOT/generate_ligand_topology.py --mol $LIG_GMX --topology $ID/pp_run.top --structure $ID/prerun.run.final.pdb --index $ID/complex.ndx --output-ligand-structure $ID/ligand.pdb --output-ligand-topology $ID/ligand.top --total-charge $ID/totalcharge.txt

            # solvate ligand-only confs
            $SINGLERUN $GMX editconf -f $ID/ligand.pdb -d $WATER_THICKNESS -bt dodecahedron -o $ID/ligand-box.pdb
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
            python3 $ABFE_ROOT/resurrect_flexible.py --flexible $ID/pp_flex.top --top $ID/ligand-ion.top --output $ID/ligand-ion-flex.top
            cat $ID/complex.ndx $ID/restr_pull.ndx > $ID/complex_with_pull.ndx
            ;;
        query,4)
            echo "DEPENDS=(3); (( PROCS = LIG_PARA ))"
            ;;
        run,4)
            do_prep_runs ligand-ion-flex ligand-ion charging-lig ligand-ion ligand ""
            ;;
        query,5)
            echo "DEPENDS=(4); (( PROCS = NCHARGE * LIG_PARA ))"
            ;;
        run,5)
            # product run for charge-discharge ligand
            do_product_runs ligand-ion charging-lig charging-lig.npt charging-lig "" "" $NCHARGE 0 0
            ;;
        query,6)
            echo "DEPENDS=(5); (( PROCS = NANNIH * LIG_PARA ))"
            ;;
        run,6)
            # product run for ligand annihilation, starting from totally discharged conformation
            do_product_runs ligand-ion annihilation-lig charging-lig.$((NCHARGE - 1))/charging-lig annihilation-lig "" "" $NANNIH 0 $(( NANNIH - 1 ))
            ;;
        query,7)
            echo "DEPENDS=(5 6); (( PROCS = LIG_PARA ))"
            ;;
        run,7)
            # re-eval for LRC for charging #0
            do_eval_run ligand-ion lr-lig charging-lig.0/charging-lig lr-lig ligand 0
            # re-eval for LRC for annihilation last
            do_eval_run ligand-ion lr-annihilation-lig annihilation-lig.$((NANNIH-1))/annihilation-lig lr-annihilation-lig ligand $((NANNIH-1))
            ;;
        query,8)
            echo "DEPENDS=(3); (( PROCS = NCHARGE * COMPLEX_PARA ))"
            ;;
        run,8)
            # product run for complex system (discharging)
            do_product_runs topol_ionized charging-complex prerun.run charging-complex complex_with_pull restr_pull.mdp $NCHARGE 0
            ;;
        query,9)
            echo "DEPENDS=(6 8); (( PROCS = NANNIH * COMPLEX_PARA ))"
            ;;
        run,9)
            # product run for complex system (annihilation)
            do_product_runs topol_ionized annihilation-complex charging-complex.$((NCHARGE-1))/charging-complex annihilation-complex complex_with_pull restr_pull.mdp $NANNIH 0 $((NANNIH - 1))
            ;;
        query,10)
            echo "DEPENDS=(3); (( PROCS = NRESTR * COMPLEX_PARA ))"
            ;;
        run,10)
            # product run for restraint decopuling (100% restraint -> 0% restraint) FEP
            do_product_runs topol_ionized restrain prerun.run restrain complex_with_pull restr_pull_decouple.mdp $NRESTR 1 $((NRESTR - 1))
            ;;
        query,11)
            echo "DEPENDS=(9 10); (( PROCS = COMPLEX_PARA ))"
            ;;
        run,11)
            # eval run (LRC) for restraint final (unrestrained)
            do_eval_run topol_ionized lr-complex restrain.$((NRESTR - 1))/restrain lr-complex complex 0
            # re-eval for LRC annihilation final
            do_eval_run topol_ionized lr-annihilation-complex annihilation-complex.$((NANNIH - 1))/annihilation-complex lr-annihilation-complex complex $((NANNIH - 1))
            ;;
        query,12)
            echo "DEPENDS=(3 4 5 6 7 8 9 10 11); NOGPU_STAGE=yes"
            ;;
        run,12)
            TEMP=$(grep ref_t mdp/run.mdp | cut -d '=' -f2)
            # do gmx bar
            do_bar charging-lig $NCHARGE
            do_bar charging-complex $NCHARGE
            do_bar annihilation-lig $NANNIH
            do_bar annihilation-complex $NANNIH
            do_bar restrain $NRESTR
            # LRC
            do_exp lr-lig charging-lig 0 $TEMP
            do_exp lr-annihilation-lig annihilation-lig $((NANNIH - 1)) $TEMP
            do_exp lr-complex restrain $((NRESTR - 1)) $TEMP
            do_exp lr-annihilation-complex annihilation-complex $((NANNIH - 1)) $TEMP
            # Charge correction
            charge_correction $ID/totalcharge.txt restrain.$((NRESTR - 1))
            # Sum evertyhing up
            $PYTHON3 $ABFE_ROOT/calc_bar_replex.py --basedir $ID --restrinfo $ID/restrinfo --temp $TEMP > $ID/result.txt
            ;;
        *)
            echo "Unexpected stage outputs"
            exit 1
            ;;
    esac
    if [[ $reqstate == run ]]; then
        echo $STEPNO >> $ID/done_step.txt
    fi
}

main

