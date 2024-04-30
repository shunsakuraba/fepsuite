import argparse
import os.path
import math
import re
import sys

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
        "charge-correction-complex": ("cc", "charge-correction", -1.0), # same as charging-complex
        "charge-correction-ligand": ("cc", "charge-correction", 1.0) # same as charging-lig
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

def avg(ls):
    return sum(ls) / len(ls)

def calc_charge_correction_impl(values):
    UNIT = 0
    VALUE = 1

    # vacuum permittivity in e^2/(kJ/mol)/nm.
    # use kJ/mol at this moment (later converted to kcal at the end)
    eps0 = 8.8541878128e-12 / 1.602176634e-19**2 / 6.02214076e+23 / 1e9 * 1e3

    retvals = []

    for isample in range(len(values["V"])):
        latticea = values["lengths"][isample][VALUE]
        xi_LS = values["xi_LS"][isample][VALUE]
        epsS = values["epsS"][isample][VALUE]

        # Note Rocklin corrected charging FE while we need to correct discharging FE (on complex state)
        q_initial = values["QS"][isample][VALUE]
        q_final = q_initial - values["QL"][isample][VALUE]
        # latticea is used because xi_LS is scaled by this number
        dgnet_usv = - xi_LS / (8. * math.pi * eps0 * epsS) * (q_final**2 - q_initial**2) / latticea

        ip_v = values["IP/V"][isample][VALUE]
        il_v = values["IL/V"][isample][VALUE]

        dgrip = ip_v * q_final - (ip_v + il_v) * q_initial

        nsol = values["Ns"][isample][VALUE]
        gammas = values["gammaS"][isample][VALUE]
        volume = values["V"][isample][VALUE]

        dgdsc = - gammas * nsol / (6 * eps0 * volume) * (q_final - q_initial)
        
        dgtot = dgnet_usv + dgrip + dgdsc
        print(f" CC {isample} dgnet_usv {dgnet_usv:10.5f}  dgrip {dgrip:10.5f}  dgdsc {dgdsc:10.5f}  dgtot {dgtot:10.5f}",
              file=sys.stderr)

        retvals.append(dgtot)
    mean = avg(retvals)
    devs = sum([(x - mean) ** 2 for x in retvals])
    var = devs / (len(retvals) - 1)

    return (mean, var)

def calc_charge_correction(fdetail):
    if os.path.exists(fdetail):
        dat = []
        with open(fdetail) as fh:
            values = {}
            pat = re.compile(r"\s*(?P<name>[A-Za-z_/]+)(?:\[(?P<unit>.+)\])?:\s*(?P<value>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)")
            for l in fh:
                m = pat.match(l.strip())
                if m is None:
                    continue
                name = m.group("name")
                unit = m.group("unit")
                value = float(m.group("value"))
                if name not in values:
                    values[name] = []
                values[name].append((unit, value))
        #print(values)
        print("""If you use this program please read and cite:
Accurate Calculation of Relative Binding Free Energies between Ligands with Different Net Charges.
Wei Chen, Yuqing Deng, Ellery Russell, Yujie Wu, Robert Abel, and Lingle Wang
Journal of Chemical Theory and Computation, 14, 6346 (2018).""",
              file=sys.stderr)
        (charge_corr, var) = calc_charge_correction_impl(values)
        return (charge_corr, var)
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

    cfactor = cycle_contribution["charge-correction-complex"][2]
    compcontr = calc_charge_correction(f"{args.basedir}/charge_correction.txt")
    ccc = (compcontr[0] * cfactor, compcontr[1] * cfactor ** 2)
    individual_bar_res["charge-correction-complex"] = ccc
    cfactor = cycle_contribution["charge-correction-ligand"][2]
    ligcontr = calc_charge_correction(f"{args.basedir}/charge_correction_lig.txt")
    ccl = (ligcontr[0] * cfactor, ligcontr[1] * cfactor ** 2)
    individual_bar_res["charge-correction-ligand"] = ccl
    for (m, vnew) in [ccc, ccl]:
        (e, v) = bar_res["charging"]
        e += m
        v += vnew
        bar_res["charging"] = (e, v)

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

