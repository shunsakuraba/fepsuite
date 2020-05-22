import argparse
import os.path
import math

cycle_contribution = {
        "charging": ("charging", 1.0),
        "restrain": ("restrain", 1.0),
        "annihilation-lig": ("annihilation", 1.0),
        "annihilation-complex": ("annihilation", -1.0)
        }
kcal_of_kj = 0.2390
gasconstant = 0.00831446  # kJ / mol / K 

def read_bar(basedir):
    ret = {}
    for f in cycle_contribution:
        with open(os.path.join(basedir, "%s.bar.log" % f)) as fh:
            for l in fh:
                ls = l.split()
                if len(ls) == 0:
                    continue
                if ls[0] != "total":
                    continue
                part = cycle_contribution[f][0]
                contribution = cycle_contribution[f][1]
                energy = float(ls[5])
                error = float(ls[7])
                if part not in ret:
                    ret[part] = (contribution * energy, error ** 2)
                else:
                    curene, curvar = ret[part]
                    curene += contribution * energy
                    curvar += error ** 2
                    ret[part] = (curene, curvar)
    return ret

def read_restr(args):
    with open(args.restrinfo) as fh:
        ls = fh.readlines()
        ls = [x for x in ls if not x.startswith("#")]
        anchor_atoms = [int(x) + 1 for x in ls[0].split()] # +1 for itp being 1-origin
        avgs = [float(x) for x in ls[1].split()]

    v0 = 1.66 # 1M standard state in nm^3
    RT = args.temp * gasconstant
    mdeltaf = math.log(8 * math.pi ** 2 * v0) + 0.5 * math.log(args.distance_spring) + 1.0 * math.log(args.angle_spring) + 1.5 * math.log(args.dihedral_spring) \
        - 2. * math.log(avgs[0]) - math.log(math.sin(avgs[1])) - math.log(math.sin(avgs[2])) - 3. * math.log(2 * math.pi * RT) # avogadro const in R cancels with kJ/mols in springs consts
    return - RT * mdeltaf

def parse_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--basedir", default=".", type=str, help="Base directory")
    parser.add_argument("--restrinfo", required=True, type=str, help="Restrain info")
    parser.add_argument("--distance-spring", type=float, default=418.68, help="Spring constant to apply (kJ/mol/nm^2)")
    parser.add_argument("--angle-spring", type=float, default=4.1868, help="Weight to the angle stdev (kJ/mol/rad^2)")
    parser.add_argument("--dihedral-spring", type=float, default=4.1868, help="Weight to the dihedral stdev (kJ/mol/rad^2)")
    parser.add_argument("--temp", type=float, default=300.0, help="Temperature")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    bar_res = read_bar(args.basedir)
    analytical_e = read_restr(args)
    restre, var = bar_res["restrain"] 
    bar_res["restrain"] = (restre - analytical_e, var) # ArVBA: P...L -> P+L, thus NEGATIVE contribution
    total = 0.
    totalvar = 0.
    for mode in sorted(bar_res.keys()):
        ene, var = bar_res[mode]
        print(mode, "%.3f" % (kcal_of_kj * ene), "%.3f" % (kcal_of_kj * math.sqrt(var)))
        total += ene
        totalvar += var
    print("total", "%.3f" % (kcal_of_kj * total), "%.3f" % (kcal_of_kj * math.sqrt(totalvar)))

