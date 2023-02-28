````
    ________________                   _ __     
   / ____/ ____/ __ \      _______  __(_) /____ 
  / /_  / __/ / /_/ /_____/ ___/ / / / / __/ _ \
 / __/ / /___/ ____/_____(__  ) /_/ / / /_/  __/
/_/   /_____/_/         /____/\__,_/_/\__/\___/ (beta)
````

# FEP-suite: scripts and pipelines for free-energy perturbation calculation works

FEP-suite is a collection of programs for calculating various free energy differences. The pipeline uses GROMACS as a backend. There are currently two pipelines:

* ABFE (absolute binding free energy) estimator, calculating binding free-energy betweeen receptor and ligand (typically protein and drug).
* FEP/REST (free energy perturbation / replica exchange solute tempering), a relative free-energy difference calculation pipeline with an enhanced sampling technique.

Under both directory named `abfe` and `feprest` you can see README.md corresponding to each pipeline and how to use. Have fun!

Note the program is currently "beta" version.

# License

GPL 3 or later (see COPYING for details)

feprest/fepgen/cmdline.h is licensed under BSD 3-clause. See the file for details.

# Authors

Main contributor:
Shun Sakuraba (National Institutes for Quantum Science and Technology, Japan)

# Acknowledgements

This software was supported by following fundings:
* a Grant-in-aid for Young Scientists (B) (grant no. 16K17778) from Japan Society for the Promotion of Science (JSPS).
* a Grant-in-aid for Scientific Research on Innovative Areas "Molecular Engine" (grant no. 19H05410) from JSPS.
* Platform Project for Supporting Drug Discovery and Life Science Research (Basis for Supporting Innovative Drug Discovery and Life Science Research (BINDS)) under Grant Number JP21am0101106, Agency for Medical Research and Development (AMED), Japan
* Project Focused on Developing Key Technology for Discovering and Manufacturing Drugs for Next-Generation Treatment and Diagnosis from AMED (JP21ae0121005s0101).
