# FEP/REST pipeline: free-energy perturbation pipeline for relative free-energy difference analysis

# What's this?

FEP/REST pipeline can be used to calculate the relative free energy change upon changing one part of molecule to the other. Examples include:
* FE change upon inducing mutation into the protein. By calculating relative change of unfolding FE, we can find which mutation stabilizes or destablizes the protein, protein complex, or protein-ligand binding.
* FE change upon changing a part of a ligand. Calculating relative change of total binding FE, we can find which ligand can bind to the protein the most.

# How to use the pipeline

In the example below I show how to calculate the stability change in terms of the Gibbs free-energy of folding. 

## Preparation 1: compiling fepgen, FEP input file generator

If you have supercomputer, this part is recommended to be run locally. You will need to `git clone` this (fepsuite) repository both locally and remotely.

You first need [Eigen3|https://eigen.tuxfamily.org/] library to compile `fepgen`. If you are using Ubuntu you can install via 

```sh
sudo apt install libeigen3-dev
```

If you are using other OS, you can just download the archive and extract it. In that case, modify `EIGENFLAGS` in `Makefile` to `include/` directory of Eigen.

After that, compile the fepgen program:

```sh
cd feprest/fepgen
make
```

If you have successfully compiled the program, you will be able to run `fepgen`:

````
% ./fepgen
usage: ./fepgen --structureA=string --structureB=string --topologyA=string --topologyB=string --structureO=string --topologyO=string [options] ... 
options:
      --help                   Print this message
  -v, --verbose                Be more verbose
  -q, --quiet                  Suppress unnecessary information
      --maxdist                Maximum distances to fit (double [=1])
      --force-bond             Cluster bond force constants (double [=500000])
…
````

This `fepgen` program is used to generate input file. We assume that you will copy generated files by `fepgen` to remote computer center via `scp` or `rsync`; how to generate files will be described later.

## Preparation 2: preparing FASPR to mutate amino acid on protein

We are currently using [FASPR](https://github.com/tommyhuangthu/FASPR) to prepare protein structure with mutation.[^1] `FASPR` can be prepared as:

````sh
git clone https://github.com/tommyhuangthu/FASPR.git
cd FASPR
g++ --ffast-math -O3 -o FASPR src/*.cpp
````

If you successfully compiled `FASPR`, you should be able to run FASPR as a command:
````
% ./FASPR
###########################################################################
                    FASPR (Version 20200309)                 
  A method for fast and accurate protein side-chain packing, 
which is an important problem in protein structure prediction
and protein design.

Copyright (c) 2020 Xiaoqiang Huang
Yang Zhang Lab
Dept. of Computational Medicine and Bioinformatics
Medical School
University of Michigan
Email:tommyhuangthu@foxmail.com, xiaoqiah@umich.edu
###########################################################################
Usage: ./FASPR -i input.pdb -o output.pdb
[-s sequence.txt] to load a sequence file
````

[^1]: You can use - and devs are actually using - other software to generate mutated structure, e.g. [pymol](https://pymol.org/), [MODELLER](https://salilab.org/modeller/), [AlphaFold2](https://github.com/deepmind/alphafold), or [PDBfixer](https://github.com/openmm/pdbfixer), to name a few.

## Preparation 3: preparing python libraries on remote

If you are using supercomputer, I recommend to run this on on remote side.

Before running the pipeline, prepare mdtraj and pymbar (ver 3.0.3) package for python3.
```sh
pip3 install mdtraj pymbar==3.0.3 --user
python3 -c "import mdtraj"   # checks whether it works correctly
```
`pymbar` must be `3.0.3` because `pymbar` frequently changes the API within minor version update. (We are planning to remove `pymbar` dependency in future.)

## Preparation 4: compiling patched GROMACS on remote

To run `feprest` pipeline you need to compile function-enhanced version of GROMACS. Currently we are only supporting GROMACS 2020.

```sh
tar xf gromacs-2020.x.tar.gz
cd gromacs-2020.x
patch -p1 < /path/to/fepsuite/feprest/gmx_patch/gromacs-2020-hrex.patch
mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs-2020-hrex -DGMX_MPI=on ..
make -j8 && make install
```

If you have trouble compiling gromacs, consult your supercomputer center's support team.

## Preparing input structure

Here, we show how to generate FEP/REST input in the protein mutation case. For other case see the advanced guide at the bottom of this document.

### Preparing Force field (not really FEP/REST-specific)
First, we prepare AMBER14SB force field file for GROMACS. [Dr. Hong Viet did excellent work of porting](https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2015-December/102497.html), and [Olomuc group added OL14/OL15 nucleic FF](https://fch.upol.cz/ff_ol/gromacs.php). But we need a few fix to complete. Let's begin from downloading:

````sh
wget https://fch.upol.cz/ff_ol/amber14sb_OL15.ff_corrected-Na-cation-params.tar.gz
tar xf amber14sb_OL15.ff_corrected-Na-cation-params.tar.gz
````

Then, first we [fix improper angle parameter](http://zhenglz.blogspot.com/2017/05/fixing-bugs-in-ff14sb-port-for-gromacs.html) in `[ dihedraltypes ]` section of improper angles:

````
NA  CV  CC  CT       4      180.00     4.60240     2
NA  CW  CC  CT       4      180.00     4.60240     2
````

We also need to workaround [the bond atom type issue](https://gitlab.com/gromacs/gromacs/-/issues/4120) existing in old (<2022) GROMACS, since we are using patched one over GROMACS 2020. We need to replace `2C` and `3C` atom type with `A2C` and `A3C`.

````
cd amber14sb_OL15.ff
sed -i 's/2C/A2C/g;s/3C/A3C/g' *
````

### Preparing FEP input file with prep_mutation_fep.py

Make a working directoy and download PDB file of T4 Lysozyme.

````sh
mkdir ~/feprest-test
cd ~/feprest-test
cp -r /path/to/amber14sb_OL15.ff .
wget https://files.rcsb.org/download/2LZM.pdb.gz
gunzip 2LZM.pdb.gz
````

Fortunately, PDB [2LZM](https://www.rcsb.org/structure/2LZM) do not have any missing residues, so we can directly feed this structure into the mutant structure prepration. 

````sh
python3 /path/to/fepsuite/feprest/utils/prep_mutation_fep.py --feprest /path/to/fepsuite/feprest --faspr /path/to/FASPR/FASPR --pdb 2LZM.pdb --mutation L99A --ff amber14sb_OL15
````

If the program finishes without problem, you should find four directories: `wt`, `L99A`, `wt_L99A` and `wt_L99A_ref`.[^2] The first one `wt` corresponds to wild-type structure (without mutation), the second one `L99A` corresponds to Leu99Ala mutation. The third one `wt_L99A` is the "mixed" structure of wild-type and Leu99Ala mutant. If you visualize `wt_L99A/conf.pdb`, you will see at residue 99 that Leu sidechain and Ala hydrogens are overlapping. 

Likewise, if you look into `topol.top` you should see:
````
TODO
````

The fourth, `wt_L99A_ref` contains the "reference state", where (TODO)

[^2]: If the program fails and you are going to rerun, you may need to remove all output directories before running `prep_mutation_fep.py`, i.e. `rm -rf wt L99A wt_L99A`.

### Copy FEP input file to remote server

```sh
scp -r wt_L99A wt_L99A_ref amber14sb_OL15.ff your-supercomputer.somewhere:somedirectory
ssh your-supercomputer.somewhere:
cd somedirectory
```

### Running pipeline on remote server 

Then, copy template file from rundir_template to start the calculation.

````sh
cp -r /path/to/fepsuite/feprest/rundir_template/* .
````

Update `run.zsh` and `para_conf.zsh` according to your environment. If you are using GPU you typically may want to update the following part:

```
# Number of MPI processes per replica. For GPU, recommended = 1.
PARA=8

# Number of threads per process. Recommended: for CPU: 1, for GPU: 2-6.
TPP=1
```

After the update, run the pipeline with

```sh
./run.zsh wt_L99A all
./run.zsh wt_L99A_ref all
````

and the pipeline will run to the end. The pipeline consists of stages 1-8. Stage 1-3 runs with single replica (minimization, NVT, NPT equilibration run), and stage 4+ runs with multiple copies of the system (NVT, NPT, replica-exchanged production run). 

If the script stops due to an error, and if you want to rerun the stage, you can rerun the stage by giving appropriate stage number. For example

````
./run.zsh wt_L99A 3 4
````

runs stage 3, waits the completion of stage 3, then runs stage 4 of `wt_L99A`. You can see the current completed stages in `wt_L99A/done_step.txt`.

### Understanding output

If everything works fine, output should be generated at `wt_L99A/bar1.log` and `wt_L99A_ref/bar1.log`. The final line of each log file represents the free-energy change. By subtracting the reference state from the calculation (such like `wt_L99A/bar.log` - `wt_L99A_ref/bar.log`), you will get the free-energy change upon mutation. If the calculation goes without problem it should be something around +1-5 kcal/mol, i.e. L99A is predicted to destabilize the protein significantly. Note the experimental value for the free-energy difference of folding is (TODO) kcal/mol ().[^3]

[^3]: If you are going to look into the article, the sign is the opposite because  et al. considered the free-energy difference of *unfolding*.

## Advanced guide: preparing FEP/REST input manually


(TODO)