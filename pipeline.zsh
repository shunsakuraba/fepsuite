#!/bin/zsh

reqstate=$1
stateno=$2
if [[ -z $stateno ]]; then
    echo "This file should be called from template.zsh" 2>&1
    exit 1
fi
shift
shift

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
        mpirun_ $NP $GMX_MPI mdrun $args
        if [[ $? != 0 ]]; then
            # fail to run. Check log file to see whether it is domain decomposition problem
            tail -20 $log_basename.log | grep -q -i "\\(domain\\|prime\\)" || { echo "Error: domain-unrelated error"; exit 1 }
            PREVNP=$NP
            # round up
            (( NP = (NP - 1) / least_unit * least_unit ))
            if (( NP == 0 )) || (( PREVNP == NP )); then
                echo "Error: unable to find proper parallel processes"
                exit 1
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
    case $cont in
        cont)
            contcmd=(-t $ID/$conf.cpt)
            ;;
        nocont)
            contcmd=()
            ;;
        *)
            echo "do_run $1 $2 $3 $4 $5 $6 $7 $8 $9 $10: \$cont invalid"
            exit 1
            ;;
    esac
    case $prep in
        prep)
            # currently nothing to do
            ;;
        product)
            ;;
        *)
            echo "do_run $1 $2 $3 $4 $5 $6 $7 $8 $9 $10: \$cont invalid"
            exit 1
            ;;
    esac
    $SINGLERUN $GMX grompp $mdp $index -p $ID/$topol.top -c $ID/$conf.pdb $contcmd -o $ID/$outprefix.$phase.tpr -po $ID/$outprefix.$phase.out.mdp $restrcmd -maxwarn $maxwarn || exit 1
    mdrun_find_possible_np 1 -deffnm $ID/$outprefix.$phase -c $ID/$outprefix.$phase.pdb
}

do_prep_runs() {
    topol=$1
    conf=$2
    output=$3
    initialrestr=$4
    ndx=$5
    pullmdp=$6
    do_run $topol $conf $output steep "$initialrestr" "$ndx" "$pullmdp" nocont prep 1 # 1 for switching
    #do_run $topol ${output}.steep $output cg "$initialrestr" "$ndx" "$pullmdp" nocont prep 0 # 1 for switching
    #do_run $topol ${output}.cg $output nvt "$initialrestr" "$ndx" "$pullmdp" nocont prep 0
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
    $PYTHON3 $ABFE_ROOT/generate_decoupling.py --template-dir $ABFE_ROOT/template-pullrestr --cominfo $ID/cominfo $mdprestr --output-mdp $ID/mdp_addenda --mode $phase -N $nrepl --mol-name $MOL
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
        $SINGLERUN $GMX grompp -f $MDP -p $ID/$topol -c $ID/$prev.pdb -t $ID/$prev.cpt -o $RUNDIR/$phase.tpr -po $ID/$output.$phase.$i.mdp $ndx_grompp -maxwarn $maxwarn
        DIRS+=$RUNDIR
    done
    mdrun_find_possible_np $nrepl -multidir $DIRS -deffnm $phase -c $phase.pdb -replex 500
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
    $SINGLERUN $GMX grompp -f ${output_trunk}-dummy.mdp -p $topol_sol -c $solute_sol -po ${output_trunk}-dummy-out.mdp -o ${output_trunk}-sol.tpr
    topol_ion=${output_trunk}-ion.top
    cp $topol_sol $topol_ion
    echo SOL | $SINGLERUN $GMX genion -s ${output_trunk}-sol.tpr -o ${output_trunk}-ion.pdb -p $topol_ion -pname $ION_POSITIVE -nname $ION_NEGATIVE -conc $IONIC_STRENGTH -neutral
}


# ---- end of subroutines for production run

main() {
    if [[ $reqstate = run ]]; then
        set -ex
    fi
    case $reqstate,$stateno in
        query,1)
            echo "DEPENDS=(); (( PROCS = NPRERUN_PARA ))"
            ;;
        run,1)
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
            $SINGLERUN $GMX grompp -f mdp/run.mdp -p $ID/topol_ionized.top -c $ID/prerun.run.pdb -t $ID/prerun.run.cpt -o $ID/pp.tpr -po $ID/pp.mdp $restrcmd -maxwarn 0 -pp $ID/pp.top
            echo "System" | $SINGLERUN $GMX trjconv -f $ID/prerun.run.xtc -o $ID/prerun.run.cut.xtc -b $RUN_PROD
            if typeset -f generate_ndx > /dev/null; then
                generate_ndx $ID/topol_ionized.top $ID/prerun.run.pdb $ID/pp.tpr $ID/index_receptor.ndx
                ndxcmd=(-n $ID/index_receptor.ndx)
            else
                echo "q" | $SINGLERUN $GMX make_ndx -f $ID/pp.tpr -o $ID/index_receptor.ndx
                ndxcmd=(-n $ID/index_receptor.ndx)
            fi
            echo "$RECEPTOR_NDX\n$RECEPTOR_NDX\nSystem" | $SINGLERUN $GMX trjconv -s $ID/prep.steep.tpr $ndxcmd -f $ID/prerun.run.pdb -pbc cluster -center -ur compact -o $ID/prerun.pbc.pdb
            python3 $ABFE_ROOT/find_restr_from_md.py --lig-sel "$LIG_MDTRAJ" --prot-sel "$RECEPTOR_MDTRAJ" --topology $ID/conf_ionized.pdb --trajectory $ID/prerun.run.cut.xtc --output $ID/restrinfo
            # Restraint for annihilation and charging
            python3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --mdp $ID/restr_pull.mdp --ndx $ID/restr_pull.ndx
            # Restraint decoupling, used in restraint phase
            python3 $ABFE_ROOT/generate_restr.py --restrinfo $ID/restrinfo --decouple-B --mdp $ID/restr_pull_decouple.mdp --ndx $ID/restr_pull_decouple.ndx
            diff -q $ID/restr_pull.ndx $ID/restr_pull_decouple.ndx || { echo "Two pull indices are not identical"; exit 1 }
            python3 $ABFE_ROOT/generate_warpdrive.py --mol $LIG_TOP --topology $ID/pp.top --structure $ID/prerun.pbc.pdb --output-structure $ID/wd.pdb --output-charging $ID/charging.top --output-ligand-q0 $ID/ligand-q0.top --output-complex-q0 $ID/complex-q0.top --output-com-info $ID/cominfo --output-index $ID/for_pull.ndx --distance $WATER_THICKNESS_CHARGING
            do_solvate $ID/wd.pdb $ID/charging.top $ID/charging
            # slice out the ligand
            python3 -c "import mdtraj; c = mdtraj.load(\"$ID/prerun.pbc.pdb\"); s = c.topology.select(\"$LIG_MDTRAJ\"); c.atom_slice(s).save_pdb(\"$ID/lig.pdb\")"
            $SINGLERUN $GMX editconf -f $ID/lig.pdb -d $WATER_THICKNESS -bt dodecahedron -o $ID/lig-box.pdb
            do_solvate $ID/lig-box.pdb $ID/ligand-q0.top $ID/ligand-q0

            # generate index files for each topology
            # 1. complex
            cat $ID/index_receptor.ndx $ID/restr_pull.ndx > $ID/complex.ndx
            # 2. warp-drive (complex + ligand), reuse complex restr_pull.ndx, and add for_pull.ndx
            echo "q" | $SINGLERUN $GMX make_ndx -f $ID/charging-ion.pdb -o $ID/charging.ndx
            cat $ID/restr_pull.ndx $ID/for_pull.ndx >> $ID/charging.ndx
            ;;
        query,4)
            echo "DEPENDS=(3)"
            ;;
        run,4)
            # prep run for charge-discharge ligand/complex warp-drive 
            do_prep_runs charging-ion charging-ion charging charging-ion charging restr_pull.mdp
            ;;
        query,5)
            echo "DEPENDS=(4); (( PROCS = NCHARGE * NCHARGE_PARA ))"
            ;;
        run,5)
            # product run for charge-discharge warp-drive system, production run
            do_product_runs charging-ion charging charging.npt charging charging restr_pull.mdp $NCHARGE 2 # 2 for copying states and absolute com pulling
            ;;
        query,6)
            echo "DEPENDS=(3)"
            ;;
        run,6)
            # prep run for ligand system (annihilation)
            do_prep_runs ligand-q0-ion ligand-q0-ion annihilation-lig ligand-q0-ion "" ""
            ;;
        query,7)
            echo "DEPENDS=(6); (( PROCS = NANN * NANN_PARA ))"
            ;;
        run,7)
            # product run for ligand system (annihilation)
            do_product_runs ligand-q0-ion annihilation-lig annihilation-lig.npt annihilation-lig "" "" $NANN 0
            ;;
        query,8)
            echo "DEPENDS=(3)"
            ;;
        run,8)
            # prep run for complex system (annihilation)
            do_prep_runs complex-q0 prerun.run annihilation-complex prerun.run complex restr_pull.mdp
            ;;
        query,9)
            echo "DEPENDS=(8); (( PROCS = NANN_PRO * NANN_PRO_PARA ))"
            ;;
        run,9)
            # product run for complex system (annihilation)
            do_product_runs complex-q0 annihilation-complex annihilation-complex.npt annihilation-complex complex restr_pull.mdp $NANN_PRO 1 # maxwarn for copying states
            ;;
        query,10)
            echo "DEPENDS=(3); (( PROCS = NCPLX_RESTR * NCPLX_RESTR_PARA ))"
            ;;
        run,10)
            # product run for restraint decopuling (100% restraint -> 0% restraint) FEP
            do_product_runs pp restrain prerun.run restrain complex restr_pull_decouple.mdp $NCPLX_RESTR 1
            ;;
        query,11)
            echo "DEPENDS=(3 5 7 9 10)"
            ;;
        run,11)
            # do gmx bar
            do_bar charging $NCHARGE
            do_bar annihilation-lig $NANN
            do_bar annihilation-complex $NANN_PRO
            do_bar restrain $NCPLX_RESTR
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

