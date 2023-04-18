# Force fields for FEP pipelines

This directory contains modified force field files.
Ideally, the pipeline programs should work with any force field files.
However, some force field files are released with known problems, or some force field files are too exotic from others and needs a modification to be run in combination with fepsuite. Hence this directory.

## `amber14sb_OL15_fs1`

This directory contains AMBER14SB (= 99SB + dihed correction) + OL14/15 (bsc0 (DNA/RNA) + OL14 (RNA) + OL15 (DNA)) force field set.
The origin and the modification are as follows:

* The origin of AMBER14SB ff file is [Dr. Hong Viet's porting](https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2015-December/102497.html) to GROMACS.
* [Olomuc group added OL14/OL15 nucleic FF](https://fch.upol.cz/ff_ol/gromacs.php).
* We added [improper angle parameters](http://zhenglz.blogspot.com/2017/05/fixing-bugs-in-ff14sb-port-for-gromacs.html) in `[ dihedraltypes ]` section of improper torsion angles to absorb difference in AMBER Leap (alphabetical ordering) and GROMACS (verbatim ordering).
* We also needed to workaround [the bond atom type issue](https://gitlab.com/gromacs/gromacs/-/issues/4120) existing in old (<2022) GROMACS, since the current feprest pipeline is using patched GROMACS 2020. We needed to replace `2C` and `3C` atom type with `A2C` and `A3C`.

## `charmm36-jul2022fs1`

This directory contains CHARMM36m 2022 Jul version with the following modifications.

* The origin of this file is [Mackerell's CHARMM36 files](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs), CHARMM36m 2022 July version.
* "ions.itp" now contain only monoatomic ions, and polyatomic ions are moved to `"ions_molecular.itp"`. This was necessary due to feprest pipeline's restrictions.
