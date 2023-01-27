To support new water model, provide an itp file with solvent included, a new itp files supporting three new moleculetypes, and coordinate generation information.

(a) *.water.itp
Water topology file. "#ifdef FLEXIBLE"~"#endif# must exist to allow stable steep/cg optimization.

(b) *.ion.itp
(1) Base atom information
"; SOLINFO N (posres) (negres)" where N is number of atoms in the water model. (posres)/(negres) are atom names & residue names of the ion.
Atom name, residue names, molecule names, are assumed to be the same.

Example:
; SOLINFO 3  NA CL

(2) Coordinate generation information
You can skip coordinate generation information if you don't use CHARGE=posonly . 
This value is used to generate water structure when randomly picked ion is perturbed into an water.
This starts from "; SOLCOORD" and 3N values follow (N is number of atoms in the solvent). These values are relative location of atoms to be generated, from the ion (unit is nm).

Example:
; SOLCOORD  0.000 0.000 0.000   0.100 0.000 0.000   -0.020 0.080 0.000
;           \---------------/   \---------------/   \---------------/
;             Relative OW pos         HW1 pos             HW2

(3) Three moleculetypes for the perturbation
SOL2pos    water -> positive ion
SOL2neg    water -> negative ion
pos2SOL    positive ion -> water (only required when CHARGE=posonly)

Note that you can't use [ settle ] because [ settle ] can only be applied to one molecule type (= original "SOL" moleculetype)

Atomtype corresponding to dummy atoms should be "PHA" (stands for PHAntom atomtype).
I recommend to use heavy hydrogen to increase calculation stability.


Files:

amber.*.itp       AMBER atom type (OW/HW), TIP3P and Joung-Cheatham monovalent ions (Na/Cl).
charmm.*.itp      CHARMM atom type (OT/HT), TIPs3P and SOD/CLA.

TODO(me):
Atom types for positive / negative monovalent ions in OPLS are opls_407 / opls_401
