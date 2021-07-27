## FEP-REST-PP: free-energy perturbation pipeline for protein mutational analysis

# Preparation 1: environment

Before running the script, prepare mdtraj and pymbar (ver 3.0.3) package for python3.
```sh
% pip3 install mdtraj pymber==3.0.3 --user
% python3 -c "import mdtraj"   # checks whether it works correctly
```
pymbar must be 3.0.3 because pymbar changes the API within minor version update. (We are planning to remove pymbar dependency)

You may also need to make sure these modules can be accessed through your job submission system (a.k.a. qsub).

Also, you need `zsh` to be installed. You can check by, for example:
```sh
% zsh --version
```

# Preparation 2: compiling GROMACS

To run this script you need to compile enhanced version of GROMACS. 

```
% tar xf gromacs-2020.x.tar.gz
% cd gromacs-2020.x
% patch -p1 < gromacs-2020-hrex.patch
% mkdir build; cd build
% cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs-2020-sphrex -DGMX_MPI=on ..
```

# Calculating the free energy
Then, follow the next order to run the FEP-REST-PP.

1. Prepare topology with `nucfepgen`.
2. Set topology file name as `topol_ionized.top` and coorinate file as `conf_ionized.pdb` or `conf_ionized.gro`. (You can use symlink for that).
3. Prepare a directory with two depths, e.g. `feprestrun/A123V`. Put topology and structure files under `feprestrun/A123V`. Similarly put reference state files under `feprestrun/A123V_ref`.
4. Clone this git repository to somewhere, e.g. `$HOME/repos/feprest`.
5. Copy `para_conf.zsh` from feprest repository to `feprestrun`. Find an appropriate submission script under `$HOME/repos/feprest/example_submit_script` to `feprestrun/`. Copy whole `mdp_template` directory to `feprestrun/`, and rename the directory name to `mdp`. If necessary, modify mdp files according to your calculation conditions.
6. Modify submission script's `FEPREST_ROOT=` part as follows:
```sh
export FEPREST_ROOT=$HOME/repos/feprest
export PIPELINE=$FEPREST_ROOT/pipeline.zsh
export GMX_DIR=$HOME/opt/gromacs-2020-sphrex
```
7. Modify `para_conf.zsh` according to your calculation environment.
8. Chdir to `feprest` directory, e.g. `cd $HOME/feprestrun`. 
9. Run the submission script with secondary directory name (in this case `A123V`) and `all`, e.g. `./titech_tsubame3.zsh A123V all`. The script will submit all necessary jobs.
10. Wait until all calculations finish.
11. If everything works fine, output should be generated at `feprestrun/A123V/bar1.log`. The final line of each log file represents the free-energy change. Typically you need to subtract the reference state from the calculation (such like `A123V/bar.log` - `A123V_ref/bar.log`)

