import argparse
import sys
import math
import copy
import cmath
import numpy
import mdtraj

def parse_and_generate_top(args):
    target_molecule = args.mol
    atom_ix_begin = {}
    atom_ix_end = {}
    with \
            open(args.topology) as fh, \
            open(args.output_charging, "w") as ofhcharging, \
            open(args.output_ligand_q0, "w") as ofhlig, \
            open(args.output_complex_q0, "w") as ofhcmplx:
        bondtype_of_atomtype = {}
        atomtype_info = {}
        dummy_atomtypes = set()
        dihtype = {}
        molecule = None
        sectiontype = None
        coulombrule = None
        residues = None
        atomnames = None
        moleculetype_section_ex_atoms = []
        natom_per_mol = {}
        atom_ix_begins = [0]

        def write_cloned_ligand():
            nonlocal moleculetype_section_ex_atoms
            if molecule == target_molecule:
                # clone the moleculetype section with modified charges
                for l in moleculetype_section_ex_atoms:
                    ofhcharging.write(l)
            moleculetype_section_ex_atoms = []

        def write_line(data, out_type):
            if "Q" in out_type:
                ofhcharging.write(data)
            if "L" in out_type:
                ofhlig.write(data)
            if "C" in out_type:
                ofhcmplx.write(data)
        
        for lraw in fh:
            while lraw.endswith("\\\n"):
                lraw = "%s %s" % (lraw[:-2], next(fh)) # this is an undocumented feature in grompp, space is implicit
            ltmp = lraw.split(';', 1)
            if len(ltmp) == 1:
                l = ltmp[0]
                comment = ""
            else:
                l = ltmp[0]
                comment = ";" + ltmp[1]
            l = l.strip()
            ls = l.split()
            if l.startswith('#'):
                sys.stderr.write("The topology file is not preprocessed")
                sys.exit(1)
            if l.startswith('['):
                sectiontype = ls[1]
                moleculetype_section_ex_atoms.append(lraw)
                write_line(lraw, "QLC")
                continue

            # blank line
            if len(ls) == 0:
                write_line(lraw, "QLC")
                continue

            if sectiontype is None:
                pass
            elif sectiontype == 'defaults':
                combrule = int(ls[1])
            elif sectiontype == 'atomtypes':
                if len(ls) < 6:
                    raise RuntimeError("Atomtype contains < 6 fields")
                # here everything is mess but toppush.cpp is actually super mess 
                if len(ls[5]) == 1 and ls[5][0].isalpha():
                    have_bonded_type = True
                    have_atomic_number = True
                elif len(ls[3]) == 1 and ls[3][0].isalpha():
                    have_bonded_type = False
                    have_atomic_number = False
                else:
                    have_bonded_type = ls[1][0].isalpha()
                    have_atomic_number = not have_bonded_type

                atomtype = ls[0]
                (mass, charge, particle, sigc6, epsc12) = ls[1 + int(have_bonded_type) + int(have_atomic_number):]

                if have_bonded_type:
                    bondtype = ls[1]
                else:
                    bondtype = atomtype
                if have_atomic_number:
                    atomic_ix = 1 + int(have_bonded_type)
                    atomic_number = int(ls[atomic_ix])
                else:
                    atomic_number = 0 # ??
                    
                # store this because we use in [ atoms ] section
                bondtype_of_atomtype[atomtype] = bondtype

                mass = float(mass)
                charge = float(charge)
                sigc6 = float(sigc6)
                epsc12 = float(epsc12)
                is_dummy = epsc12 == 0.

                atomtype_info[atomtype] = (charge, mass)
                if is_dummy:
                    dummy_atomtypes.add(atomtype)

                # in WD-type topology generator module it is not necessary to print scaled one
                # and the original line is also appended after this big if-elif-else block
            elif sectiontype == 'moleculetype':
                write_cloned_ligand()
                molecule = ls[0]
                # These None are sentinels for 1-origin access
                bondtype_list = [None]
                scaled = [None]
                atomnames = []
                residues = []
                moleculetype_section_ex_atoms = []
                if molecule == target_molecule:
                    moleculetype_section_ex_atoms.append("[ moleculetype ]\n")
                    moleculetype_section_ex_atoms.append("%s_ref %s\n" % (molecule, ls[1]))
                assert molecule not in natom_per_mol
                natom_per_mol[molecule] = 0
                write_line(lraw, "QLC")
                continue
            elif sectiontype == 'atoms':
                aindex = int(ls[0])
                atomtype = ls[1]
                bondtype_list.append(bondtype_of_atomtype[atomtype])

                # charge & mass is optional parameters, oof...
                if len(ls) > 6:
                    charge = float(ls[6])
                else:
                    (charge, _) = atomtype_info[atomtype]
                if len(ls) > 7:
                    mass = float(ls[7])
                else:
                    (_, mass) = atomtype_info[atomtype]
                if len(ls) > 8:
                    atomtypeB = ls[8]
                    if atomtypeB.endswith(args.suffix):
                        sys.stderr.write("Warning: typeB in atom line should not be suffixed (molecule %s, atom index %d)\n" % (molecule, aindex))
                    (chargeB, massB) = atomtype_info[atomtypeB]
                    fep = True
                    if is_scaled and not atomtypeB.endswith(args.suffix):
                        #atomtypeB += args.suffix
                        pass
                else:
                    fep = False
                    chargeB = 0.
                    massB = 0.
                    atomtypeB = None
                if len(ls) > 9:
                    chargeB = float(ls[9])
                if len(ls) > 10:
                    massB = float(ls[10])
                if fep:
                    raise RuntimeError("FEP-ed topology is unsupported")

                _resnr = ls[2]
                resid = ls[3]
                residues.append(resid)
                atomname = ls[4]
                atomnames.append(atomname)
                
                if molecule == target_molecule:
                    write_line("%5d %4s %4s %4s %4s %5s %16.8e %16.8e %4s %16.8e %16.8e%s\n" % 
                            # charging: charge 0->alpha
                            (aindex, atomtype, ls[2], ls[3], ls[4], ls[5], 0, mass, atomtype, charge, mass, comment.rstrip()),
                            "Q")
                    write_line("%5d %4s %4s %4s %4s %5s %16.8e %16.8e%s\n" % 
                            # ligand / complex: charge 0 (fixed)
                            (aindex, atomtype, ls[2], ls[3], ls[4], ls[5], 0, mass, comment.rstrip()),
                            "LC")
                    # cloned: charge alpha->0, added to charging phase later
                    moleculetype_section_ex_atoms.append("%5d %4s %4s %4s %4s %5s %16.8e %16.8e %4s %16.8e %16.8e\n" %
                            (aindex, atomtype, ls[2], ls[3], ls[4], ls[5], charge, mass, atomtype, 0, mass))
                else:
                    # all [ moleculetype ] sections appear in all topology (only [ molecules ] will be changed)
                    write_line("%5d %4s %4s %4s %4s %5s %16.8e %16.8e%s\n" % 
                            (aindex, atomtype, ls[2], ls[3], ls[4], ls[5], charge, mass, comment.rstrip()),
                            "QLC")

                natom_per_mol[molecule] += 1
                continue
            elif sectiontype == "system":
                write_cloned_ligand()
            elif sectiontype == "molecules":
                molname = ls[0]
                nmol = int(ls[1])
                if molname == target_molecule:
                    assert nmol == 1
                    write_line(lraw, "QLC")
                else:
                    write_line(lraw, "QC")
                prevend = atom_ix_begins[-1]
                newend = prevend + natom_per_mol[molname] * nmol
                atom_ix_begins.append(newend)
                atom_ix_begin[molname] = prevend
                atom_ix_end[molname] = newend
                continue
            if molecule != None:
                moleculetype_section_ex_atoms.append(lraw)
            # With a few exceptions, we just print as-is. when write_line should not be executed, do "continue" at the end of if-clauses.
            write_line(lraw, "QLC")
        
        # at the very bottom of charging topology add reference system
        write_line("%s_ref 1\n" % target_molecule, "Q")

    print("Topology generation completed. atoms = [", atom_ix_begin[target_molecule], ",", atom_ix_end[target_molecule], ")")
    print("generating index addenda")
    with open(args.output_index, "w") as ofh:
        print("[ grp-complex ]", file=ofh)
        for i in range(atom_ix_begins[0], atom_ix_begins[-1]):
            print(i + 1, file=ofh)
        print(file=ofh)
        print("[ grp-lig ]", file=ofh)
        for i in range(atom_ix_begins[-1], atom_ix_begins[-1] + atom_ix_end[target_molecule] - atom_ix_begin[target_molecule]):
            print(i + 1, file=ofh)
        print(file=ofh)
    print("Completed")

    return (atom_ix_begin[target_molecule], atom_ix_end[target_molecule])

def estimate_radius(substructure):
    com = mdtraj.compute_center_of_mass(substructure)
    ddsq = numpy.amax(numpy.sum((substructure.xyz[0, :, :] - com[:, numpy.newaxis, :]) ** 2, axis=2))
    return math.sqrt(ddsq)

def generate_coords(structure, ix_begin, ix_end, dbuffer):
    rprotlig = estimate_radius(structure)
    lig = structure.atom_slice([i for i in range(ix_begin, ix_end)])
    rlig = estimate_radius(lig)
    protlig_com = mdtraj.compute_center_of_mass(structure)
    lig_com = mdtraj.compute_center_of_mass(lig)
    print("r(complex) =", rprotlig, "r(lig) =", rlig)

    # place in X direction, in (buffer/2, protlig, buffer, lig, buffer/2) order
    box_dist = numpy.array([rprotlig * 2 + rlig * 2 + dbuffer * 2, rprotlig * 2 + dbuffer, rprotlig * 2 + dbuffer])
    new_protlig_com = numpy.array([dbuffer / 2 + rprotlig, dbuffer / 2 + rprotlig, dbuffer / 2 + rprotlig])
    new_lig_com = numpy.array([dbuffer * 3 / 2 + rprotlig * 2 + rlig, dbuffer / 2 + rprotlig, dbuffer / 2 + rprotlig])

    newprotligcrd = structure.xyz[:, :, :] - protlig_com[:, numpy.newaxis, :] + new_protlig_com[numpy.newaxis, numpy.newaxis, :]
    print(newprotligcrd.shape)
    newligcrd = lig.xyz[:, :, :] - lig_com[:, numpy.newaxis, :] + new_lig_com[numpy.newaxis, numpy.newaxis, :]
    print(newligcrd.shape)
    new_coords = numpy.concatenate((newprotligcrd, newligcrd), axis=1)

    return (new_coords, box_dist, new_protlig_com, new_lig_com)

def generate_pdb(args, ix_begin, ix_end):
    print("Loading PDB")
    target_molecule = args.mol
    structure = mdtraj.load(args.structure)
    topology = structure.topology

    print("Calcing new coordinates")
    (newcoords, box_dist, new_protlig_com, new_lig_com) = generate_coords(structure, ix_begin, ix_end, args.distance)
    # new com information is necessary for generating restraints

    print("Generating new PDB structure")
    ligtopology = mdtraj.Topology()
    ligchain = ligtopology.add_chain()
    prevres = None
    ligres = None
    for aix in range(ix_begin, ix_end):
        atom = topology.atom(aix)
        res = atom.residue
        if prevres != res:
            ligres = ligtopology.add_residue(res.name, ligchain, res.resSeq, res.segment_id)
        prevres = res
        ligtopology.add_atom(atom.name, atom.element, ligres)

    print("Writing new PDB structure")
    newtopology = topology.join(ligtopology)
    newstructure = mdtraj.Trajectory(newcoords, newtopology, unitcell_lengths = box_dist, unitcell_angles=[90.0, 90.0, 90.0])
    newstructure.save_pdb(args.output_structure)

    print("PDB generation completed")
    return (new_protlig_com, new_lig_com)

def init_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--topology", type=str, required=True, help="Input topology file")
    parser.add_argument("--structure", type=str, required=True, help="Input PDB formatted file. molecules are assumed to be whole and clustered, and chains are not splitted.")
    parser.add_argument("--output-structure", type=str, required=True, help="output PDB file")
    parser.add_argument("--mol", type=str, required=True, help="Target molecule name")
    parser.add_argument("--distance", type=float, default=1.0, help="Buffer region distance")
    parser.add_argument("--output-charging", type=str, required=True, help="Output topology file with distant charging")
    parser.add_argument("--output-ligand-q0", type=str, required=True, help="Output topology file for the ligand with 0-charged.")
    parser.add_argument("--output-complex-q0", type=str, required=True, help="Output topology file for the complex with 0-charged.")
    parser.add_argument("--output-com-info", type=str, required=True, help="Output COM position used in simulation")
    parser.add_argument("--output-index", type=str, required=True, help="Output index file addendum for com pulling")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    ix_begin, ix_end = parse_and_generate_top(args)
    (new_protlig_com, new_lig_com) = generate_pdb(args, ix_begin, ix_end)
    print("Generating cominfo")
    with open(args.output_com_info, "w") as ofh:
        print("# complex", file=ofh)
        print(*new_protlig_com, file=ofh)
        print("# ligand", file=ofh)
        print(*new_lig_com, file=ofh)
    print("Successfully Generated")
