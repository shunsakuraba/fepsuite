import re
import numpy
import os.path
import argparse
import math

modes = ["charging", "restraint", "annihilation-lig", "annihilation-complex"]

def parse_repl_ex(f):
    findpat_lam = re.compile(r"\s*fep-lambdas\s*=")
    findpat = re.compile("Repl +average probabilities:")
    ret_lambda = None
    with open(f) as fh:
        for l in fh:
            if findpat_lam.match(l):
                lambdas = l.split("=")[1]
                if lambdas.strip() == "TRUE":
                    continue
                ret_lambda = [float(l) for l in lambdas.split()]
            if findpat.match(l):
                break
        l = next(fh) # Replica #s
        l = next(fh) # actual exchange probabilities
        exs = [float(x) for x in l.split()[1:]]
    if ret_lambda is None:
        raise RuntimeError("lambda does not appear in log files")
    return (ret_lambda, exs)

def optimize_state_from_exprobs(n, prevlambda, exprobs):
    nlexval = [- math.log(max(e, 0.01)) for e in exprobs]

    cumuval = [] 
    cur = 0.0
    for i in range(n - 1):
        cumuval.append(cur)
        cur += nlexval[i]
    cumuval.append(cur)

    print("cumu = ", cumuval)

    prop = cumuval[-1] / (n - 1.0)
    result = [0.0]
    for i in range(1, n - 1):
        level = float(i) * prop
        #print("level = ", level)
        ix = 0
        for j in range(0, n):
            if cumuval[j] > level:
                ix = j - 1
                break
        #print("ix = ", ix)
        remain = level - cumuval[ix]
        newlam = remain / (cumuval[ix + 1] - cumuval[ix]) * (prevlambda[ix + 1] - prevlambda[ix]) + prevlambda[ix]
        result.append(newlam)
    result.append(1)
    print("after opt:", result)
    return result

def update_params(oldcoord, newcoords, nth):
    blend = (0.7 ** nth)
    retcoord = []
    for (i, x) in enumerate(newcoords):
        v = ((1. - blend) * oldcoord[i] +
             blend * x)
        retcoord.append(v)
    return retcoord

def update_lambda(_mode, n, update, prerunphase):
    prevlambda, exprob = parse_repl_ex(update)
    assert len(prevlambda) == n
    assert len(exprob) == (n - 1)

    newlambda = optimize_state_from_exprobs(n, prevlambda, exprob)
    return update_params(prevlambda, newlambda, prerunphase)

def lambda_schedule(mode, n, update, prerunphase):
    if update is not None:
        return update_lambda(mode, n, update, prerunphase)
    # return lambda values for each mode
    # actually for all cases we just return with equal distribution
    if mode == "charging":
        return numpy.linspace(0., 1., n, endpoint=True)
    elif mode == "restraint":
        return numpy.linspace(0., 1., n, endpoint=True)
    elif mode == "annihilation":
        # We tested SSC(2) with alpha=0.2 but it was *very* bad indeed, tested with replica exchange rate
        # SSC(2) may be only applicable with charge sc + vdw sc
        # SSC(2), Lee et al. JCTC 16 5512 (2020).
        #x = numpy.linspace(0., 1., n, endpoint=True)
        #return (6. * x**2 - 15. * x + 10.) * x ** 3
        #x = numpy.linspace(0., 1., n, endpoint=True)
        #return (-2. * x ** 3 + 3. * x ** 2)
        #x = numpy.linspace(0., 1., n, endpoint=True)
        #return 0.75 * x + 0.25 * (-2. * x ** 3 + 3. * x ** 2)
        # With alpha=0.2, this gives a good starting point.
        x = numpy.linspace(0., 1., n, endpoint=True)
        return (1. - (x - 1.) ** 2)

    raise RuntimeError("Unsupported mode")

def generate_decoupling(args):
    mode_target = args.mode
    mode = mode_target.split('-')[0]

    # read template file, small enough
    with open(os.path.join(args.template_dir, mode_target + ".mdp")) as fh:
        template = fh.read()
    if args.additional_mdp is not None:
        # add here if special handling needed
        if False:
            removelines = re.compile("")
        else:
            removelines = re.compile("dummy-never-match")
        with open(args.additional_mdp) as fh:
            for l in fh:
                if not removelines.match(l):
                    template += l
    # just feed into format function

    lambda_values = lambda_schedule(mode, args.N, args.update, args.update_nth)
    lambda_values_str = " ".join(["%.4f" % l for l in lambda_values])
    params = {
            "lambdas_formatted": lambda_values_str,
            "group_mol": args.mol_name,
            }
    for i in range(args.N):
        params["lambda_state"] = i
        with open(os.path.join(args.output_mdp, "%s-%d.mdp" % (mode_target, i)), "w") as ofh:
            ofh.write(template.format(**params))

def init_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mol-name", type=str, default="MOL", help="target molecul name")
    parser.add_argument("--output-mdp", type=str, default='.', help="mdp output directory")
    parser.add_argument("--additional-mdp", type=str, default=None, help="Additional mdp file (typically pull)")
    parser.add_argument("--mode", type=str, required=True, help="Stage (charging, restrain, annihilation-lig, or annihilation-complex)")
    parser.add_argument("--update", type=str, help="Update lambda values with this log file")
    parser.add_argument("--update-nth", type=int, help="pre-run phase no (1-origin)")
    parser.add_argument("-N", type=int, required=True, help="Number of states")
    parser.add_argument("--template-dir", type=str, required=True, help="Template directory")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    generate_decoupling(args)



