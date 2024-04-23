import argparse
import os.path
import math

cycle_contribution = {
        "charging-lig": ("bar", "charging", 1.0),
        "charging-complex": ("bar", "charging", -1.0),
        "restraint-analytical": ("restrcor", "restraint", 1.0),
        "restraint": ("bar", "restraint", 1.0),
        "annihilation-lig": ("bar", "annihilation", 1.0),
        "annihilation-complex": ("bar", "annihilation", -1.0),
        "lr-lig": ("lrc", "long-range-correction", -1.0),
        "lr-annihilation-lig": ("lrc", "long-range-correction", 1.0),
        "lr-complex": ("lrc", "long-range-correction", 1.0),
        "lr-annihilation-complex": ("lrc", "long-range-correction", -1.0),
        "charge-correction-complex": ("cc", "charge-correction", 1.0),
        "charge-correction-ligand": ("cc", "charge-correction", -1.0)
        }
kcal_of_kj = 0.23900574
gasconstant = 0.00831446  # kJ / mol / K 

def read_bar(basedir):
    ret = {}
    individual = {}
    for f, (ty, part, coeff) in cycle_contribution.items():
        if ty != "bar":
            continue
        with open(os.path.join(basedir, "%s.bar.log" % f)) as fh:
            for l in fh:
                ls = l.split()
                if len(ls) == 0:
                    continue
                if ls[0] != "total":
                    continue
                energy = float(ls[5])
                error = float(ls[7])
                individual[f] = (coeff * energy, error ** 2)
                if part not in ret:
                    ret[part] = (coeff * energy, error ** 2)
                else:
                    curene, curvar = ret[part]
                    curene += coeff * energy
                    curvar += error ** 2
                    ret[part] = (curene, curvar)
    return (ret, individual)

def read_update_lrc(basedir, agg, individual):
    for f, (ty, part, coeff) in cycle_contribution.items():
        if ty != "lrc":
            continue
        with open(os.path.join(basedir, "%s.lrc.txt" % f)) as fh:
            result = fh.readlines()[-1]
            ls = [float(x) for x in result.split()]
            df = ls[0] * coeff
            dfstd = ls[1]
            dfvar = dfstd ** 2
            individual[f] = (df, dfvar)
            if part not in agg:
                agg[part] = (0., 0.)
            curene, curvar = agg[part]
            curene += df
            curvar += dfvar
            agg[part] = (curene, curvar)

def read_restr(args):
    with open(args.restrinfo) as fh:
        ls = fh.readlines()
        ls = [x for x in ls if not x.startswith("#")]
        anchor_atoms = [int(x) + 1 for x in ls[0].split()] # +1 for itp being 1-origin
        avgs = [float(x) for x in ls[1].split()]

    v0 = 1.6605391 # 1M standard state in nm^3
    RT = args.temp * gasconstant # kJ/mol
    mdeltaf = math.log(8 * math.pi ** 2 * v0) + 0.5 * math.log(args.distance_spring) + 1.0 * math.log(args.angle_spring) + 1.5 * math.log(args.dihedral_spring) \
        - 2. * math.log(avgs[0]) - math.log(math.sin(avgs[1])) - math.log(math.sin(avgs[2])) - 3. * math.log(2 * math.pi * RT) # avogadro const in R cancels with kJ/mols in springs consts
    return - RT * mdeltaf # kJ/mol

def calc_charge_correction(fsummary):
    if os.path.exists(fsummary):
        dat = []
        with open(fsummary) as fh:
            for l in fh:
                if l.startswith("#"):
                    continue
                ls = l.split()
                dat.append(float(ls[-1]))
        if len(dat) == 0:
            avg = 0.
            var = 0.
        else:
            avg = sum(dat) / len(dat)
            var = 0.
            for d in dat:
                var += (d - avg) ** 2
            if len(dat) <= 1:
                var = 0.
            else:
                var /= (len(dat) - 1)
        return (avg, var)
    else:
        return (0., 0.)

def parse_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--basedir", default=".", type=str, help="Base directory")
    parser.add_argument("--restrinfo", required=True, type=str, help="Restrain info")
    parser.add_argument("--distance-spring", type=float, default=4184.0, help="Spring constant to apply (kJ/mol/nm^2)")
    parser.add_argument("--angle-spring", type=float, default=41.848, help="Weight to the angle stdev (kJ/mol/rad^2)")
    parser.add_argument("--dihedral-spring", type=float, default=41.84, help="Weight to the dihedral stdev (kJ/mol/rad^2)")
    parser.add_argument("--temp", type=float, default=300.0, help="Temperature")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    (bar_res, individual_bar_res) = read_bar(args.basedir)
    read_update_lrc(args.basedir, bar_res, individual_bar_res)
    analytical_e = read_restr(args)
    restre, var = bar_res["restraint"] 
    bar_res["restraint"] = (restre - analytical_e, var) # ArVBA in Boresch 2003 paper corresponds to  P...L -> P+L, thus NEGATIVE contribution
    individual_bar_res["restraint-analytical"] = (-analytical_e, 0)
    individual_bar_res["charge-correction-complex"] = calc_charge_correction(f"{args.basedir}/charge_correction.txt")
    ligcontr = calc_charge_correction(f"{args.basedir}/charge_correction_lig.txt")
    individual_bar_res["charge-correction-ligand"] = (-ligcontr[0], ligcontr[1])

    total = 0.
    totalvar = 0.
    for mode in sorted(individual_bar_res.keys()):
        ene, var = individual_bar_res[mode]
        print(mode, "%.3f" % (kcal_of_kj * ene), "%.3f" % (kcal_of_kj * math.sqrt(var)))
        total += ene
        totalvar += var
    print("----")
    for mode in sorted(bar_res.keys()):
        ene, var = bar_res[mode]
        print(mode, "%.3f" % (kcal_of_kj * ene), "%.3f" % (kcal_of_kj * math.sqrt(var)))
    print("----")
    print("total", "%.3f" % (kcal_of_kj * total), "%.3f" % (kcal_of_kj * math.sqrt(totalvar)))

