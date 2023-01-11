#!/bin/zsh

reqstate=$1
stateno=$2
if [[ -z $stateno ]]; then
    echo "This file should be called from template.zsh" 2>&1
    exit 1
fi
shift
shift

if [[ -z $ABFE_ROOT ]]; then
    echo "ABFE_ROOT is not set"
    exit 1
fi
if [[ -z $REST2PY_ROOT ]]; then
    echo "REST2PY_ROOT is not set"
    exit 1
fi

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


# ---- end of subroutines for production run

main() {
    if [[ $reqstate = run ]]; then
        #set -e
        #trap '{ echo "Aborting job"; set +e; exit $ERRORCODE }' ZERR
        set -x
    fi
    case $reqstate,$stateno in
        query,1)
            echo "DEPENDS=(); (( PROCS = COMPLEX_PARA ))"
            ;;
        run,1)
            # min for state A
            $SINGLERUN $GMX grompp -f mdp/steep.mdp -c $ID/conf_ionized -p $ID/topol_ionized.top -r $ID/conf_ionized -o $ID/steep -po $ID/steep.mdout -maxwarn 1 # 1 for switching
            mdrun_find_possible_np 1 -deffnm $ID/steep -c $ID/steep.pdb
            $SINGLERUN $GMX grompp -f mdp/cg.mdp -c $ID/steep.pdb -p $ID/topol_ionized.top -r $ID/conf_ionized -o $ID/min -po $ID/min.mdout -maxwarn 1
            mdrun_find_possible_np 1 -deffnm $ID/min -c $ID/min.pdb
            ;;
        query,2)
            echo "DEPENDS=(1); (( PROCS = COMPLEX_PARA ))"
            ;;
        run,2)
            sed '/dt/c dt=0.0005' mdp/nvt.mdp > $ID/nvt05.mdp
            $SINGLERUN $GMX grompp -f $ID/nvt05.mdp -c $ID/min.pdb -p $ID/topol_ionized.top -r $ID/conf_ionized -o $ID/nvt05 -po $ID/nvt05.mdout -maxwarn 1
            mdrun_find_possible_np 1 -deffnm $ID/nvt05 -c $ID/nvt05.pdb
            $SINGLERUN $GMX grompp -f mdp/nvt.mdp -c $ID/min.pdb -p $ID/topol_ionized.top -r $ID/conf_ionized -o $ID/nvt -po $ID/nvt.mdout -t $ID/nvt05.cpt -maxwarn 1
            mdrun_find_possible_np 1 -deffnm $ID/nvt -c $ID/nvt.pdb
            ;;
        query,3)
            echo "DEPENDS=(2); (( PROCS = COMPLEX_PARA ))"
            ;;
        run,3)
            $SINGLERUN $GMX grompp -f mdp/npt.mdp -c $ID/nvt.pdb -p $ID/topol_ionized.top -t $ID/nvt.cpt -r $ID/conf_ionized -o $ID/npt -po $ID/npt.mdout -maxwarn 1
            mdrun_find_possible_np 1 -deffnm $ID/npt -c $ID/npt.pdb
            ;;
        query,4)
            echo "DEPENDS=(3); (( PROCS = 1 ))"
            ;;
        run,4)
            $SINGLERUN $GMX grompp -f mdp/run.mdp -c $ID/npt.pdb -t $ID/npt.cpt -p $ID/topol_ionized.top -o $ID/pp -po $ID/runpp.mdout -pp $ID/topol_pp.top -maxwarn 1
            python3 $REST2PY_ROOT/canonicalize_top.py $ID/topol_pp.top $ID/topol_pp_ca.top
            python3 $ABFE_ROOT/tools/initial-grest/rest_region.py --structure $ID/steep.pdb --target-gmx "$LIG_GMX" --receptor "$RECEPTOR_MDTRAJ" --range $GREST_DISTANCE --input $ID/topol_pp_ca.top --output $ID/topol.underlined.top
            TEMP=$(grep "^\\s*ref[_-]t\\s*=" mdp/run.mdp | cut -d '=' -f2 | cut -d ';' -f1 | tr -d ' ')
            if [[ -z $TEMP ]]; then 
                TEMP=300.0
            fi
            python3 $REST2PY_ROOT/grestdih_nohrex.py --temp0 $TEMP --temp $TEMP_HIGH $ID/topol.underlined.top $ID/topol.fep.top
            LAMBDAS=()
            for i in {0..$((NGREST-1))}; do
                # equispaced.
                (( LAMBDA = i / (NGREST - 1.0) )) || true
                LAMBDAS+=$LAMBDA
            done
            for i in {0..$((NGREST-1))}; do
                wdir=$ID/rest$i
                mkdir $wdir || true
                sed -e "s/%LAMBDAS%/$LAMBDAS/;s/%LAMBDASTATE%/$i/" mdp/run_fep.mdp > $wdir/run.mdp
                $SINGLERUN $GMX grompp -f $wdir/run.mdp -c $ID/npt.pdb -t $ID/npt.cpt -p $ID/topol.fep.top -o $wdir/run -po $wdir/run.po.mdp -maxwarn 1
            done
            ;;
        query,5)
            echo "DEPENDS=(4); (( PROCS = NGREST * COMPLEX_PARA ))"
            ;;
        run,5)
            reps=()
            for i in {0..$((NGREST-1))}; do
                wdir=$ID/rest$i
                reps+=$wdir
            done
            mdrun_find_possible_np $NGREST -deffnm prodrun -multidir $reps -s run.tpr -replex 500
            ;;
        query,6)
            echo "DEPENDS=(5); (( PROCS = 1 ))"
            ;;
        run,6) 
            base_structure=$ID/conf_ionized.pdb
            [[ -e $base_structure ]] || base_structure=$ID/conf_ionized.gro
            python3 $ABFE_ROOT/make_ndx.py --structure $base_structure --topology $ID/topol_pp_ca.top --ligand $LIG_GMX --receptor $RECEPTOR_MDTRAJ --output $ID/analysis.ndx
            echo "Receptor\nSystem\n" | $SINGLERUN $GMX trjconv -s $ID/rest0/run.tpr -f $ID/rest0/prodrun.xtc -o $ID/rest0/prodrun.pbc.xtc -n $ID/analysis.ndx -center -pbc mol -ur compact
            echo "Receptor\nReceptor\nSystem\n" | $SINGLERUN $GMX trjconv -s $ID/rest0/run.tpr -f $ID/rest0/prodrun.pbc.xtc -o $ID/rest0/prodrun.fit.xtc -n $ID/analysis.ndx -center -fit rot+trans
            echo "Ligand\nSystem\n" | $SINGLERUN $GMX cluster -s $ID/rest0/run.tpr -f $ID/rest0/prodrun.fit.xtc -b $RUN_PROD -method jarvis-patrick -nofit -g $ID/cluster.log -n $ID/analysis.ndx -o $ID/rmsd-clust -dist $ID/rmsd-dist
            python3 $ABFE_ROOT/tools/initial-grest/select_cluster.py --base-structure $base_structure --traj $ID/rest0/prodrun.pbc.xtc --cluster-log $ID/cluster.log --threshold $CLUSTER_THRESHOLD --output-prefix $ID/result
            N=$(head -n 1 $ID/result.txt) || true
            if (( N > 0 )); then
                for i in {0..$((N-1))}; do
                    echo "Receptor\nLigand\n" | $SINGLERUN $GMX rms -s $ID/steep.tpr -f $ID/result$i.pdb -o $ID/rms$i.xvg -n $ID/analysis.ndx
                done
            fi
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
