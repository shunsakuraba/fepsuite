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
    set +e
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
    set -e
    PROCS_SAVED=$NP
    echo "Normal termination of mdrun"
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
    trap '{ echo "Aborting job"; set +e; exit $ERRORCODE }' ZERR
    set -x
    $SINGLERUN $GMX grompp $mdp $index -p $ID/$topol.top -c $ID/$conf.pdb $contcmd -o $ID/$outprefix.$phase.tpr -po $ID/$outprefix.$phase.out.mdp $restrcmd -maxwarn $maxwarn
    mdrun_find_possible_np 1 -deffnm $ID/$outprefix.$phase -c $ID/$outprefix.$phase.pdb $nstlist_cmd
}

do_prep_runs() {
    topol=$1
    conf=$2
    output=$3
    initialrestr=$4
    ndx=$5
    pullmdp=$6
    trap '{ echo "Aborting job"; set +e; exit $ERRORCODE }' ZERR
    set -x
    do_run $topol $conf $output steep "$initialrestr" "$ndx" "$pullmdp" nocont prep 1 # 1 for switching
    do_run $topol ${output}.steep $output nvt "$initialrestr" "$ndx" "$pullmdp" nocont prep 0 # due to instability of SETTLE, and inability to use -DFLEXIBLE, only perform steepest descent.
    do_run $topol ${output}.nvt $output npt "$initialrestr" "$ndx" "$pullmdp" cont prep 1 # 1 for gen-vel
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
    DIRS=()
    mkdir -p $ID/mdp_addenda || true
    mkdir -p $ID/mdp || true
    mdprestr=()
    if [[ -n $pullmdp ]]; then
        mdprestr=(--additional-mdp $ID/$pullmdp)
    fi
    $PYTHON3 $ABFE_ROOT/generate_decoupling.py --template-dir $TEMPLATE --cominfo $ID/cominfo $mdprestr --output-mdp $ID/mdp_addenda --mode $phase -N $nrepl --mol-name $LIG_TOP || return 1
    for i in {0..$((nrepl-1))}; do
        CODE=$phase.$i
        RUNDIR=$ID/$CODE
        mkdir -p $RUNDIR || true
        MDP=$ID/mdp/$phase-$i.mdp
        cat ./mdp/run.mdp $ID/mdp_addenda/$phase-$i.mdp > $MDP
        ndx_grompp=()
        if [[ -n $ndx ]]; then
            ndx_grompp=(-n $ID/$ndx)
        fi
        $SINGLERUN $GMX grompp -f $MDP -p $ID/$topol -c $ID/$prev.pdb -t $ID/$prev.cpt -o $RUNDIR/$phase.tpr -po $ID/$output.$phase.$i.mdp $ndx_grompp -maxwarn $maxwarn || return 1
        DIRS+=$RUNDIR
    done
    nstlist_cmd=()
    if [[ -n $NSTLIST ]]; then
        nstlist_cmd=(-nstlist $NSTLIST)
    fi
    mdrun_find_possible_np $nrepl -multidir $DIRS -deffnm $phase -c $phase.pdb -replex 500 $nstlist_cmd || return 1
}

do_eval_run() {
    topol=$1
    phase=$2
    prev=$3
    output=$4
    ndx=$5
    pullmdp=$6
    maxwarn=$7

    { grep -v "dispcorr" ./mdp/run.mdp | grep -v 'nstxtcout'; cat $TEMPLATE/$phase.mdp } > $ID/$phase.mdp
    $SINGLERUN $GMX grompp -f $ID/$phase.mdp -p $ID/$topol -c $ID/$prev.pdb -t $ID/$prev.cpt -o $ID/$phase.tpr -po $ID/$output.$phase.mdp $ndx_grompp -maxwarn $maxwarn || return 1
    mdrun_find_possible_np 1 -deffnm $ID/$phase -rerun $ID/$prev.xtc $nstlist_cmd
}

do_bar() {
    mode=$1
    nrepl=$2
    $SINGLERUN $GMX bar -b $RUN_PROD -f $ID/$mode.{0..$((nrepl-1))}/$mode.xvg -o $ID/$mode.bar_diff.xvg > $ID/$mode.bar.log
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
    $SINGLERUN $GMX grompp -f ${output_trunk}-dummy.mdp -p $topol_sol -c $solute_sol -po ${output_trunk}-dummy-out.mdp -o ${output_trunk}-sol.tpr || return 1
    topol_ion=${output_trunk}-ion.top
    cp $topol_sol $topol_ion
    echo SOL | $SINGLERUN $GMX genion -s ${output_trunk}-sol.tpr -o ${output_trunk}-ion.pdb -p $topol_ion -pname $ION_POSITIVE -nname $ION_NEGATIVE -conc $IONIC_STRENGTH -neutral || return 1
}

generate_ndx() {
    conf=$1
    index=$2
    output=$3

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
    python3 $ABFE_ROOT/make_ndx.py --structure $conf --index $index --ligand $LIG_GMX --receptor $RECEPTOR_MDTRAJ --output $output
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
    fi
}

# ---- end of subroutines for production run

main() {
    if [[ $reqstate = run ]]; then
        trap '{ echo "Aborting job"; set +e; exit $ERRORCODE }' ZERR
        set -x
    fi
    case $reqstate,$stateno in
        query,1)
            echo "DEPENDS=(); (( PROCS = NPRERUN_PARA ))"
            ;;
        run,1)
            touch mdp/dummy.mdp
            # TODO: check whether maxwarn 0 is sufficient
            $SINGLERUN $GMX grompp -f mdp/dummy.mdp -p $ID/topol_ionized.top -c $ID/conf_ionized -o $ID/pp.tpr -po $ID/pp.mdp -maxwarn 0 -pp $ID/pp.top
            echo "q" | $SINGLERUN $GMX make_ndx -f $ID/pp.tpr -o $ID/index_vanilla.ndx
            generate_ndx $ID/conf_ionized $ID/index_vanilla.ndx $ID/complex.ndx
            do_prep_runs topol_ionized conf_ionized prep conf_ionized "" ""
            ;;
        query,2)
            echo "DEPENDS=(1); (( PROCS = NPRERUN_PARA ))"
            ;;
        run,2)
            do_run topol_ionized prep.npt prerun run "" "" "" cont prep 0
            ;;
        query,3)
            echo "DEPENDS=(2);PROCS=1"
            ;;
        run,3)
            # Compute RMSD of ligands for thresholding
            echo "Receptor\nSystem" | $SINGLERUN $GMX trjconv -s $ID/prerun.run.tpr -f $ID/prerun.run.xtc -o $ID/prerun.run.recpbc.xtc -b $RUN_PROD -center -pbc mol -n $ID/complex.ndx -ur compact
            echo "Receptor\nLigand" | $SINGLERUN $GMX rms -s $ID/prerun.run.pdb -f $ID/prerun.run.recpbc.xtc -o $ID/prerun.rms.fromfinal.xvg -n $ID/complex.ndx
            python3 $ABFE_ROOT/rms_check.py --rms=$ID/prerun.rms.fromfinal.xvg --threshold=$EQ_RMSD_CUTOFF || { echo "RMS of ligands too large, aborting the calculation" 1>&2; false }
            python3 $ABFE_ROOT/find_restr_from_md.py --lig-sel "Ligand" --prot-sel "Receptor" --index $ID/complex.ndx --topology $ID/conf_ionized.pdb --trajectory $ID/prerun.run.recpbc.xtc --output $ID/restrinfo
            # Restraint for annihilation and charging
            python3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --mdp $ID/restr_pull.mdp --ndx $ID/restr_pull.ndx
            # Restraint decoupling, used in restraint phase
            python3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --decouple-B --mdp $ID/restr_pull_decouple.mdp --ndx $ID/restr_pull_decouple.ndx
            diff -q $ID/restr_pull.ndx $ID/restr_pull_decouple.ndx || { echo "Two pull indices are not identical" 1>&2; false }

            # generate topologies and conformations:
            # stage 4, 5, 6: ligand discharge calculation (+LRC) \
            # stage 7: ligand annihilation calculation           +-> for these two, generate ligand only structure and topologies
            # stage 8, 9: complex discharge calculation          \
            # stage 10: complex annihilation calculation         +-> for these three, complex structures as well as topologies are used
            # stage 11, 12: restrain calculation (+LRC)          /
            
            # prepare pp because we calculate ligand's total charge here
            $SINGLERUN $GMX grompp -f mdp/run.mdp -p $ID/topol_ionized.top -c $ID/prerun.run.pdb -t $ID/prerun.run.cpt -o $ID/pp.tpr -po $ID/pp.mdp $restrcmd -maxwarn 0 -pp $ID/pp_run.top
            python3 $ABFE_ROOT/generate_ligand_topology.py --mol $LIG_GMX --topology $ID/topol_ionized.top --structure $ID/prerun.pbc.pdb --index $ID/complex.ndx --output-ligand-structure $ID/ligand.pdb --output-ligand-topology $ID/ligand.top --total-charge $ID/totalcharge.txt

            # solvate ligand-only confs
            $SINGLERUN $GMX editconf -f $ID/lig.pdb -d $WATER_THICKNESS -bt dodecahedron -o $ID/lig-box.pdb
            do_solvate $ID/lig-box.pdb $ID/ligand.top $ID/ligand

            # generate index files for each topology
            # complex: index can be reused
            cat $ID/complex.ndx $ID/restr_pull.ndx > $ID/complex_pull.ndx
            # ligand: not needed because never pulled
        
            # generate mdp files for input
            for d in $TEMPLATE/*-lig.mdp; do
                cp $d $ID/
            done
            for d in $TEMPLATE/*-complex.mdp; do
                if [[ ${d##**/} = lr-complex.mdp ]]; then
                    cp $d $ID/
                else
                    cat $d $ID/restr_pull.mdp > $ID/${d:t}
                fi
            done
            ;;
        query,4)
            echo "DEPENDS=(3)"
            ;;
        run,4)
            # prep run for charge-discharge ligand
            do_prep_runs ligand ligand ligand.charging ligand ligand charging-lig
            ;;
        query,5)
            echo "DEPENDS=(4); (( PROCS = NCHARGE * LIG_PARA ))"
            ;;
        run,5)
            # product run for charge-discharge ligand
            do_product_runs ligand charging-lig ligand.npt ligand.charging "" charge-lig $NCHARGE 0 0 # FIXME: what was the number I wanted to add to do_product_runs functionality?? -> Saving trajectory for LRC, enabled on 0 (unmodded ligand)
            ;;
        query,6)
            echo "DEPENDS=(5)"
            ;;
        run,6)
            # re-eval for LRC. Sample 0th trajectory
            do_eval_run ligand ligand.charging.0 lr-lig ligand 0 # FIXME TODO: do_eval_run must also remove dispcorr
            ;;
        query,7)
            echo "DEPENDS=(5); (( PROCS = NANNIH * LIG_PARA ))"
            ;;
        run,7)
            # product run for ligand annihilation, starting from totally discharged conformation
            do_product_runs ligand ligand ligand.charging.$((NCHARGE - 1))/ligand.charging annihilation-lig ligand.annihilation ligand annihilation-lig $NCHARGE 0 
            ;;
        query,8)
            echo "DEPENDS=(3)"
            ;;
        run,8)
            # prep run for complex system (discharging)
            do_prep_runs topol_ionized prerun.run charging-complex prerun.run complex charging-complex
            ;;
        query,9)
            echo "DEPENDS=(8); (( PROCS = NCHARGE * COMPLEX_PARA ))"
            ;;
        run,9)
            # product run for complex system (discharging)
            do_product_runs topol_ionized charging-complex charging-complex.npt charging-complex complex charging-complex $NCHARGE 0
            ;;
        query,10)
            echo "DEPENDS=(9); (( PROCS = NANNIH * COMPLEX_PARA ))"
            ;;
        run,10)
            # product run for complex system (annihilation)
            do_product_runs topol_ionized annihilation-complex charging-complex.$((NCHARGE-1))/charging-complex complex restr_pull.mdp $NCHARGE 0
            ;;
        query,11)
            echo "DEPENDS=(3); (( PROCS = NCPLX_RESTR * NCPLX_RESTR_PARA ))"
            ;;
        run,11)
            # product run for restraint decopuling (100% restraint -> 0% restraint) FEP
            do_product_runs topol_ionized restrain prerun.run restrain complex restr_pull_decouple.mdp $NCPLX_RESTR 1 $((NCPLX_RESTR - 1))
            ;;
        query,12)
            echo "DEPENDS=(11)"
            ;;
        run,12)
            # eval run (LRC) for restraint final (unrestrained)
            do_eval_run topol_ionized restrain.$((NCPLX_RESTR - 1))/restrain complex-lr "" "" 0
            ;;
        query,13)
            echo "DEPENDS=(3 5 7 9 10)"
            ;;
        run,13)
            # do gmx bar
            do_bar charging-lig $NCHARGE
            do_bar charging-complex $NCHARGE
            do_bar annihilation-lig $NANNH
            do_bar annihilation-complex $NANNH
            do_bar restrain $NCPLX_RESTR
            # LRC
            do_exp lr-lig charging-lig.0 
            do_exp lr-complex restrain.$((NCPLX_RESTR - 1))
            # Charge correction
            charge_correction $ID/totalcharge.txt restrain.$((NCPLX_RESTR - 1))
            # Sum evertyhing up
            TEMP=$(grep ref_t mdp/run.mdp | cut -d '=' -f2)
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

