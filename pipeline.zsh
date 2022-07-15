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
    #set +e
    is_log=0
    multidir=()
    for arg in $args; do
        if [[ $is_log = 1 ]]; then
            log_basename=$arg
        elif [[ $is_multidir = 1 ]]; then
            multidir+=$arg
        fi
        case $arg in
            -deffnm|-l)
                is_log=1
                ;;
            -multidir)
                is_multidir=1
                ;;
            -*)
                is_multidir=0
                is_log=0
                ;;
            *)
                ;;
        esac
    done
    case $log_basename in
        *.log)
            log_basename=${log_basename%%.log}
            ;;
    esac
    if [[ -z $multidir ]]; then
        multidir=(.)
    fi
    if [[ -n $PROCS_SAVED ]]; then
        NP=$PROCS_SAVED
    else
        NP=$PROCS
    fi
    ntomp=()
    if [[ -z $OMP_NUM_THREADS ]] || (( OMP_NUM_THREADS == 1 )); then
        ntomp=(-ntomp 1) # this is required for forcing GROMACS to run
    fi
    while true; do
        echo "Trying with NP=$NP"
        mpirun_ $NP $GMX_MPI mdrun $args $NSTLIST_CMD $ntomp
        if [[ $? != 0 ]]; then
            # fail to run. Check log file to see whether it is domain decomposition problem
            domain_error=0
            for d in $multidir; do
                if tail -20 $log_basename.log | grep -q -i "\\(domain\\|prime\\)"; then
                    domain_error=1
                fi
            done
            if (( domain_error = 0 )); then
                echo "Error: domain-unrelated error"
                exit $ERRORCODE 
            fi
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

do_bar ()
{
    ID=$1
    MINPART=$2
    MAXPART=$3
    TEMP=$4
    mkdir $ID/bar || true
    python3 $FEPREST_ROOT/bar_deltae.py --minpart $MINPART --maxpart $MAXPART --xvgs $ID/prodrun/rep%sim/deltae.part%part.xvg  --nsim $NREP --temp $TEMP --save-dir $ID/bar | tee $ID/bar1.log || echo "BAR failed due to bad convergence, please continue the run to get it fixed"
}

prep_gpucpu_parameters ()
{
    if (( GPP > 0 )); then
        NB_WHICH=(-nb gpu)
    else
        NB_WHICH=(-nb cpu)
    fi
}

initref() {
    REFCMDINIT=()
    if [[ -n $REFINIT ]]; then
        REFCMDINIT=(-r $ID/$REFINIT)
    fi
    REFCMD=()
    if [[ -n $REFCRD ]]; then
        REFCMD=(-r $ID/$REFCRD)
    fi
}

main() {
    if [[ $reqstate = run ]]; then
        #set -e
        #trap '{ echo "Aborting job"; set +e; exit $ERRORCODE }' ZERR
        set -x
    fi
    INITTIME=$(date +%s)
    initref

    case $reqstate,$stateno in
        query,all)
            echo {1..8}
            ;;
        query,1)
            echo "DEPENDS=(); (( PROCS = PARA ))"
            ;;
        run,1)
            # min for state A
            sed -e "/%LAMBDA%/d;/%VDWLAMBDA%/d;/%STATE%/d;/free-energy/c free-energy = no" mdp/cginit.mdp > $ID/minA.mdp
            sed -e "/integrator/c integrator = steep" $ID/minA.mdp > $ID/steepA.mdp
            # now that grompp emits warning for -Dfoo if #ifdef foo does not exist
            if ! grep -q POSRES $ID/$BASETOP ; then
                sed -e "s/-DPOSRES/ /" -i $ID/minA.mdp
                sed -e "s/-DPOSRES/ /" -i $ID/steepA.mdp
            fi
            $SINGLERUN $GMX grompp -f $ID/steepA.mdp -c $ID/$BASECONF -p $ID/$BASETOP -o $ID/steepA -po $ID/steepA.mdout -maxwarn $((BASEWARN+0)) $REFCMDINIT
            mdrun_find_possible_np 1 -deffnm $ID/steepA
            $SINGLERUN $GMX grompp -f $ID/minA.mdp -c $ID/steepA.gro -p $ID/$BASETOP -o $ID/minA -po $ID/minA.mdout -maxwarn $((BASEWARN+0)) $REFCMDINIT
            mdrun_find_possible_np 1 -deffnm $ID/minA
            ;;
        query,2)
            echo "DEPENDS=(1); (( PROCS = PARA ))"
            ;;
        run,2)
            # NVT run
            sed -e "/%LAMBDA%/d;/%VDWLAMBDA%/d;/%STATE%/d;/free-energy/c free-energy = no" mdp/nvtinit.mdp > $ID/nvtA.mdp
            if ! grep -q POSRES $ID/$BASETOP ; then
                sed -e "s/-DPOSRES/ /" -i $ID/nvtA.mdp
            fi
            sed '/^\[ atoms \]/,/^\[ \.\+ \]/s/1.0080\+e+00/8.00000e+00/g' $ID/$BASETOP > $ID/heavy.top
            $SINGLERUN $GMX grompp -f $ID/nvtA.mdp -c $ID/minA.gro -p $ID/heavy.top -o $ID/nvtA -po $ID/nvtA.mdout -maxwarn $((BASEWARN+1)) $REFCMDINIT
            mdrun_find_possible_np 1 -deffnm $ID/nvtA
            ;;
        query,3)
            echo "DEPENDS=(2); (( PROCS = PARA ))"
            ;;
        run,3)
            # NPT run (10 ns)
            cp mdp/nptinit.mdp $ID/nptA.mdp
            $SINGLERUN $GMX grompp -f $ID/nptA.mdp -c $ID/nvtA.gro -p $ID/heavy.top -o $ID/nptA -po $ID/nptA.mdout -maxwarn $((BASEWARN+1)) -pp $ID/fep_pp.top $REFCMD
            mdrun_find_possible_np 1 -deffnm $ID/nptA 
            ;;
        query,4)
            echo "DEPENDS=(3); (( PROCS = PARA ))"
            ;;
        run,4)
            # initialize for state A to B
            python3 $FEPREST_ROOT/add_underline.py $ID/$BASECONF $ID/fep_pp.top $ID/fep_underlined.top
            prev=$ID/nptA
            top=$ID/fep_underlined.top
            if [[ $CHARGED = yes ]]; then
                python3 $FEPREST_ROOT/neutralize.py $ID/fep_underlined.top $ID/nptA.gro $ID/fep_underlined_neut.top $ID/nptA_neut.gro $AT_POSITIVE $AT_NEGATIVE
                prev=$ID/nptA_neut
                top=$ID/fep_underlined_neut.top
            fi
            python3 $FEPREST_ROOT/underlined_group.py $top $ID/for_rest.ndx
            python3 $FEPREST_ROOT/rest2py/replica_optimizer.py init $NREP feprest --basedir $ID
            mkdir $ID/genmdps || true
            python3 $FEPREST_ROOT/rest2py/replica_optimizer.py update-mdp mdp/cg.mdp $ID/genmdps/cg%d.mdp --basedir $ID
            mkdir $ID/gentops || true
            python3 $FEPREST_ROOT/rest2py/replica_optimizer.py update-topology $top $ID/gentops/fep_%d.top --basedir $ID
            ln -s $FEPREST/itp_addenda/*.itp $ID || true
            # this is a hack to enable #include "foo.itp" or "../foo.itp" in the top file. FIXME: how to deal with this?
            ln -s $PWD/$ID/*.itp gentops || true 
            ln -s $PWD/*.itp $ID || true
            for i in {0..$((NREP - 1))}; do
                work=$ID/min$i
                mkdir $work || true
                sed -e "s/cg/steep/;/nsteps/s/5000/500/;" $ID/genmdps/cg$i.mdp > $work/steep$i.mdp
                python3 $FEPREST_ROOT/re-tip3p.py $ID/gentops/fep_$i.top $ID/gentops/fep_tip3p_$i.top
                # turn atoms heavy. This is a bit hacky and may break with custom force fields
                sed -i '/^\[ atoms \]/,/^\[ bonds \]/s/1.0080\+e+00/8.00000e+00/g' $ID/gentops/fep_tip3p_$i.top
                $SINGLERUN $GMX grompp -f $work/steep$i.mdp -c $prev -p $ID/gentops/fep_tip3p_$i.top -o $work/steep$i -po $work/steep.mdout.$i -maxwarn $((BASEWARN+1)) $REFCMD
                mdrun_find_possible_np 1 -deffnm $work/steep$i -rdd $DOMAIN_SHRINK
                $SINGLERUN $GMX grompp -f $ID/genmdps/cg$i.mdp -c $work/steep$i.gro -p $ID/gentops/fep_tip3p_$i.top -o $work/min$i -po $work/min.mdout.$i -maxwarn $((BASEWARN+1)) $REFCMD
                mdrun_find_possible_np 1 -deffnm $work/min$i -rdd $DOMAIN_SHRINK
                prev=$work/min$i
            done
            ;;
        query,5)
            echo "DEPENDS=(4); (( PROCS = PARA * NREP ))"
            ;;
        run,5)
            # tune replex
            top=$ID/fep_underlined.top
            if [[ $CHARGED = yes ]]; then
                top=$ID/fep_underlined_neut.top
            fi
            prevgro=()
            for i in {0..$((NREP-1))}; do
                prevgro+=$ID/min$i/min$i.gro
            done
            for p in {1..$NTUNE}; do
                work=$ID/nvt$p
                mkdir $work || true
                python3 $FEPREST_ROOT/rest2py/replica_optimizer.py update-mdp mdp/nvt.mdp $work/nvt${p}_%d.mdp --basedir $ID
                python3 $FEPREST_ROOT/rest2py/replica_optimizer.py update-topology $top $work/fep_%d.top --basedir $ID
                reps=()
                for i in {0..$((NREP - 1))}; do
                    mrundir=$work/rep$i
                    mkdir $mrundir || true
                    reps+=$mrundir
                    python3 $FEPREST_ROOT/re-tip3p.py $work/fep_$i.top $work/fep_tip3p_$i.top
                    sed -i '/^\[ atoms \]/,/^\[ bonds \]/s/1.0080\+e+00/8.00000e+00/g' $work/fep_tip3p_$i.top
                    { echo "energygrps = hot"; echo "userint1 = 1" } >> $work/nvt${p}_$i.mdp 
                    $SINGLERUN $GMX grompp -f $work/nvt${p}_$i.mdp -c ${prevgro[$((i+1))]} -p $work/fep_tip3p_$i.top -o $mrundir/nvt -po $mrundir/nvt.mdout -maxwarn $((BASEWARN+1)) $REFCMD -n $ID/for_rest.ndx
                done
                mdrun_find_possible_np $NREP -deffnm nvt -multidir $reps -rdd $DOMAIN_SHRINK # extend to 100 ps
                for i in {0..$((NREP - 1))}; do
                    mrundir=$work/rep$i
                    $SINGLERUN $GMX convert-tpr -s $mrundir/nvt -o $mrundir/nvt_c -extend 50
                done
                prep_gpucpu_parameters # sets NB_WHICH
                mdrun_find_possible_np $NREP -deffnm nvt -multidir $reps -s nvt_c -cpi nvt -hrex -replex 100 -rdd $DOMAIN_SHRINK $NB_WHICH -bonded cpu # extend 50 ps
                python3 $FEPREST_ROOT/rest2py/replica_optimizer.py optimize $work/rep0/nvt.log --basedir $ID --step $((i + 1))
                cat $ID/replica_states
                prev=$work
                prevgro=()
                for i in {0..$((NREP-1))}; do
                    prevgro+=${reps[$((i+1))]}/nvt.gro
                done
            done
            cp $prev/fep_tip3p_{0..$((NREP-1))}.top $ID/gentops/
            ;;
        query,6)
            echo "DEPENDS=(5); (( PROCS = PARA * NREP ))"
            ;;
        run,6)
            # NPT run
            work=$ID/npt
            mkdir $work || true
            python3 $FEPREST_ROOT/rest2py/replica_optimizer.py update-mdp mdp/npt.mdp $work/npt%d.mdp --basedir $ID
            reps=()
            for i in {0..$((NREP - 1))}; do
                mrundir=$work/rep$i
                mkdir $mrundir || true
                reps+=$mrundir
                #echo "energygrps = hot" >> $work/npt$i.mdp
                $SINGLERUN $GMX grompp -f $work/npt$i.mdp -c $ID/nvt${NTUNE}/rep$i/nvt.gro -t $ID/nvt${NTUNE}/rep$i/nvt.cpt -p $ID/gentops/fep_tip3p_$i.top -o $mrundir/npt -po $mrundir/npt.mdout -maxwarn $((BASEWARN+1)) $REFCMD
            done
            mdrun_find_possible_np $NREP -deffnm npt -multidir $reps -rdd $DOMAIN_SHRINK
            ;;
        query,7)
            echo "DEPENDS=(6); (( PROCS = PARA * NREP ))"
            ;;
        run,7)
            # 50 ps initialization
            mkdir $ID/run.mdout || true
            mkdir $ID/runmdps || true
            python3 $FEPREST_ROOT/rest2py/replica_optimizer.py update-mdp mdp/run.mdp $ID/runmdps/run%d.mdp --basedir $ID
            work=$ID/prodrun
            mkdir $work || true
            reps=()
            for i in {0..$((NREP - 1))}; do
                mrundir=$work/rep$i
                mkdir $mrundir || true
                reps+=$mrundir
                { echo "energygrps = hot"; echo "userint1 = 1" } >> $ID/runmdps/run$i.mdp
                $SINGLERUN $GMX grompp -f $ID/runmdps/run$i.mdp -c $ID/npt/rep$i/npt.gro -t $ID/npt/rep$i/npt.cpt -p $ID/gentops/fep_tip3p_$i.top -o $mrundir/run -po $mrundir/run.mdout -maxwarn $((BASEWARN+1)) $REFCMD -n $ID/for_rest.ndx
                {
                    sed -e "/nstxout/d" $ID/runmdps/run$i.mdp
                    echo "nstxout=0\nnstxtcout=0\nnstvout=0\nnstdhdl=0"
                } > $ID/runmdps/eval$i.mdp
                $SINGLERUN $GMX grompp -f $ID/runmdps/eval$i.mdp -c $ID/npt/rep$i/npt.gro -t $ID/npt/rep$i/npt.cpt -p $ID/gentops/fep_tip3p_$i.top -o $mrundir/eval -po $mrundir/eval.mdout -maxwarn $((BASEWARN+1)) $REFCMD -n $ID/for_rest.ndx
            done
            mdrun_find_possible_np $NREP -deffnm run -multidir $reps -rdd $DOMAIN_SHRINK
            mkdir $ID/checkpoint_7 || true
            for d in $reps; do
                mkdir -p $ID/checkpoint_7/$d || true
                cp $d/run.cpt $ID/checkpoint_7/$d
                
                cp $d/run.cpt $d/prodrun.cpt
                cp $d/run.tpr $d/run_ph0.tpr
            done 
            ;;
        # step 99: for analysis
        query,99)
            echo "DEPENDS=(); (( PROCS = 1 ))"
            ;;
        run,99)
            for state in A B; do
                case $state in
                A)
                    REPNO=0 || true
                    ;;
                B)
                    REPNO=$((NREP - 1))
                    ;;
                esac
                ipart=1
                ndxfile=$ID/prodrun/fepbase_$state.ndx
                if [[ ! -e $ndxfile ]]; then
                    python3 $FEPREST_ROOT/make_ndx_trjconv_analysis.py -i $ID/fepbase_$state.pdb -o $ndxfile
                fi
                while true; do
                    (( ipart += 1 ))
                    ipartstr=$(printf '%04d' $ipart)
                    sourcefile=$ID/prodrun/rep$REPNO/prodrun.part$ipartstr.trr
                    destfile=$ID/prodrun/state$state.part$ipartstr.xtc
                    if [[ ! -e $sourcefile ]]; then
                        break
                    fi
                    if [[ -e $destfile ]] && [[ $destfile -nt $sourcefile ]]; then
                        continue
                    fi
                    echo "centering\noutput" | $SINGLERUN $GMX trjconv -s $ID/prodrun/rep$REPNO/run.tpr -f $sourcefile -o $destfile -pbc atom -ur compact -center -n $ndxfile
                done
            done
            ;;
        query,*)
            (( PREV = stateno - 1 ))
            echo "DEPENDS=($PREV); (( PROCS = PARA * NREP ))"
            ;;
        run,*)
            # extend run
            (( PHASE = stateno - 7 )) || true
            TEMP=$(grep ref_t mdp/run.mdp | cut -d '=' -f2)
            work=$ID/prodrun
            reps=()
            for i in {0..$((NREP - 1))}; do
                mrundir=$work/rep$i
                mkdir $mrundir || true
                reps+=$mrundir
                $SINGLERUN $GMX convert-tpr -s $mrundir/run_ph$((PHASE-1)) -o $mrundir/run_ph$PHASE -extend $SIMLENGTH
            done
            mdrun_find_possible_np $NREP -deffnm prodrun -s run_ph$PHASE -cpi prodrun -cpt 60 -hrex -othersim deltae -othersiminterval $SAMPLING_INTERVAL -multidir $reps -replex $REPLICA_INTERVAL -noappend -rdd $DOMAIN_SHRINK $NB_WHICH -bonded cpu

            mkdir $ID/checkpoint_$stateno || true
            for d in $reps; do
                mkdir -p $ID/checkpoint_$stateno/$d || true
                cp $d/prodrun.cpt $ID/checkpoint_$stateno/$d
            done
            MINPART=2
            (( MAXPART = PHASE + 1 ))
            do_bar $ID $MINPART $MAXPART $TEMP
        ;;
    esac
    if [[ $reqstate == run ]]; then
        echo $STEPNO >> $ID/done_step.txt
        CURTIME=$(date +%s)
        echo "Pipeline stage $stateno finished in $(( CURTIME - INITTIME )) sec"
    fi
}

main

