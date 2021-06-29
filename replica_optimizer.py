#!/usr/bin/python3

# Ad-hock replica exchange parameter optimizer

import sys
import math
import argparse
import re
import os
import os.path
import collections
import subprocess

replica_state_file = "replica_states"
State = collections.namedtuple("State", ["mode", "n", "mid_integer_points", "coordinates"])
RestParams = collections.namedtuple("RestParams", ["temp0", "temp", "general_lambda", "charge_lambda_A", "charge_lambda_B", "vdw_lambda", "end_restrain_A", "end_restrain_B"])

def get_script_dir():
    return os.path.dirname(sys.argv[0])

def parse_repl_ex(f):
    findpat = re.compile("Repl  average probabilities:")
    with open(f) as fh:
        for l in fh:
            if findpat.match(l):
                break
        l = next(fh) # Replica #s
        l = next(fh) # actual exchange probabilities
        exs = [float(x) for x in l.split()[1:]]
    return exs

def load_state(basedir):
    statefile = os.path.join(basedir, replica_state_file)
    with open(statefile) as fh:
        l = next(fh)
        ls = l.split()
        assert ls[0] == "mode" and ls[1] == "="
        mode = ls[2]

        l = next(fh)
        ls = l.split()
        assert ls[0] == "n" and ls[1] == "="
        n = int(ls[2])

        l = next(fh)
        ls = l.split()
        assert ls[0] == "mid-points" and ls[1] == "="
        mid_points = int(ls[2])
        
        l = next(fh)
        ls = l.split()
        assert ls[0] == "coordinates" and ls[1] == "="
        coordinates = [float(x) for x in ls[2:]]
    state = State(mode=mode, n=n,
                  mid_integer_points=mid_points,
                  coordinates=coordinates)
    return state

def save_state(state, basedir):
    statefile = os.path.join(basedir, replica_state_file)
    assert type(state.mode) == str
    assert type(state.n) == int
    assert type(state.mid_integer_points) == int
    assert type(state.coordinates) == list
    assert len(state.coordinates) == state.n
    assert abs(state.coordinates[0] - 0.0) < 1e-8
    assert abs(state.coordinates[-1] - (state.mid_integer_points + 1)) < 1e-8
    with open(statefile, 'w') as fh:
        print("mode = %s" % state.mode, file=fh)
        print("n = %d" % state.n, file=fh)
        print("mid-points = %d" % state.mid_integer_points, file=fh)
        print("coordinates = %s" % (" ".join([str(x) for x in state.coordinates])),
              file=fh)

def reset_midpoints(coords, intstates):
    # (idx, prev, new)
    anchor = [(0, 0., 0.)]
    for i in range(1, intstates):
        nearest = -1
        dist = 1e9
        for (j, x) in enumerate(coords):
            r = abs(x - i)
            if r < dist:
                dist = r
                nearest = j
        anchor.append((nearest, coords[nearest], float(i)))
        coords[nearest] = float(i)

    anchor.append((len(coords) - 1, float(intstates), float(intstates)))
    # rescale values according to anchors
    for i in range(1, intstates):
        (b, prevx, newx) = anchor[i - 1]
        (e, prevy, newy) = anchor[i]
        for j in range(b + 1, e):
            coords[j] = (coords[j] - prevx) / (prevy - prevx) * (newy - newx) + newx

    assert abs(coords[0] - 0.0) < 1e-8
    assert abs(coords[-1] - intstates) < 1e-8

    # fix to be just that value
    coords[0] = 0.0
    coords[-1] = float(intstates)
    return coords

def init_state(nstate, intermeds):
    coords = []
    for i in range(nstate):
        #x = 0.5 * (1 - math.cos(float(i) / (nstate - 1) * 180.0 * math.pi / 180.0))
        x = float(i) / (nstate - 1)
        coords.append(x * intermeds)
    coords = reset_midpoints(coords, intermeds)
    print(coords)
    return coords
            
def do_init(nstate, mode, basedir):
    if mode == "feprest":
        state_trans = 3
        coords = init_state(nstate, state_trans)
        # use
        state = State(mode = mode,
                      n = nstate,
                      mid_integer_points = state_trans - 1,
                      coordinates = coords)
    save_state(state, basedir)

def optimize_state_from_exprobs(state, probs):
    nlexval = [- math.log(max(e, 0.01)) for e in probs]

    cumuval = [] 
    cur = 0.0
    for i in range(state.n - 1):
        cumuval.append(cur)
        cur += nlexval[i]
    cumuval.append(cur)
    #print("cumu = ", cumuval)

    prop = cumuval[-1] / (state.n - 1.0)
    result = [0.0]
    for i in range(1, state.n - 1):
        level = float(i) * prop
        #print("level = ", level)
        ix = 0
        for j in range(0, state.n):
            if cumuval[j] > level:
                ix = j - 1
                break
        #print("ix = ", ix)
        remain = level - cumuval[ix]
        newlam = remain / (cumuval[ix + 1] - cumuval[ix]) * (state.coordinates[ix + 1] - state.coordinates[ix]) + state.coordinates[ix]
        #print("lam = ", newlam)
        result.append(newlam)
    result.append(float(state.mid_integer_points + 1))
    return result

def update_params(oldcoord, newcoords):
    blend = 0.5
    retcoord = []
    for (i, x) in enumerate(newcoords):
        v = ((1. - blend) * oldcoord[i] +
             blend * x)
        retcoord.append(v)
    return retcoord

def do_optimize_step(logfile, basedir):
    state = load_state(basedir)
    exprobs = parse_repl_ex(logfile)
    new_coords = optimize_state_from_exprobs(state, exprobs)
    new_coords = update_params(state.coordinates, new_coords)
    new_coords = reset_midpoints(new_coords, state.mid_integer_points + 1)

    state = state._replace(coordinates = new_coords)
    save_state(state, basedir)

def get_parameter_lists(state):
    paramlist = []
    if state.mode == "feprest":
        basetemp = 300 # real temperature does not matter
        maxtemp = 1200
        rest_nmax = (state.n - 1) // 2 # [0 .. rest_nmax]
        for i in range(state.n):
            rest_state = i
            if rest_state > rest_nmax:
                rest_state = state.n - 1 - i
            assert rest_state <= rest_nmax
            temp = basetemp * math.exp(math.log(maxtemp / basetemp) / rest_nmax
                                       * rest_state)
            c = state.coordinates[i]
            assert state.mid_integer_points == 2
            maxlambda = 3.0
            general_lambda = c / maxlambda
            if c <= 1.0:
                cx = c
                vdw_lambda = 0.0
                charge_lambda_A = cx
                charge_lambda_B = 0.0
                end_restrain_A = 0.0
                end_restrain_B = math.sin(cx * 90.0 * math.pi / 180.0)
            elif c > 1.0 and c < 2.0:
                cx = c - 1.0
                vdw_lambda = cx
                charge_lambda_A = 1.0
                charge_lambda_B = 0.0
                end_restrain_A = 0.0
                end_restrain_B = 1.0
            else:
                cx = c - 2.0
                vdw_lambda = 1.0
                charge_lambda_A = 1.0
                charge_lambda_B = cx
                end_restrain_A = 1. - math.cos(cx * 90.0 * math.pi / 180.0)
                end_restrain_B = 1.0

            paramlist.append(RestParams(temp0=basetemp, temp=temp,
                                        general_lambda=general_lambda,
                                        charge_lambda_A=charge_lambda_A,
                                        charge_lambda_B=charge_lambda_B,
                                        vdw_lambda=vdw_lambda,
                                        end_restrain_A=end_restrain_A,
                                        end_restrain_B=end_restrain_B))
            
    return paramlist

def call_rest2py(top_pp, top_generate, param):
    rest2py = os.path.join(get_script_dir(), "rest2py.py")
    python = "python3"
    
    # Python 3.5 may not be available, so keep using the old interface.
    # (Tsubame3, I'm looking at you!)
    subprocess.check_call([python, rest2py,
                           "--unify-charge",
                           "--charge-lambda", "%.7f" % param.general_lambda,
                           "--charge-lambda-A", "%.7f" % param.charge_lambda_A,
                           "--charge-lambda-B", "%.7f" % param.charge_lambda_B,
                           "--end-restrain-dihedralA", "%.7f" % param.end_restrain_A,
                           "--end-restrain-dihedralB", "%.7f" % param.end_restrain_B,
                           "--temp0", "%.7f" % param.temp0,
                           "--temp", "%.7f" % param.temp,
                           top_pp, top_generate])

def replace_mdp_template(mdp_template, mdp_generate, param):
    with open(mdp_template) as fh, open(mdp_generate, "w") as ofh:
        general_lambda = "%.7f" % param.general_lambda
        vdw_lambda = "%.7f" % param.vdw_lambda
        for l in fh:
            l = l.replace("%LAMBDA%", general_lambda)
            l = l.replace("%VDWLAMBDA%", vdw_lambda)
            ofh.write(l)

def do_update_topology(top_pp, top_generate, basedir):
    state = load_state(basedir)
    params = get_parameter_lists(state)
    if state.mode == "feprest":
        for i in range(state.n):
            call_rest2py(top_pp, top_generate % i, params[i])
    
def do_update_mdp(mdp_template, mdp_generate, basedir):
    state = load_state(basedir)
    params = get_parameter_lists(state)
    if state.mode == "feprest":
        for i in range(state.n):
            replace_mdp_template(mdp_template, mdp_generate % i, params[i])

def parse_args():
    parser = argparse.ArgumentParser(description="""Optimize the replica exchange parameters for FEP+REST2 calculation. REST2 parameters are indeed not optimized at all.

Typical usage:
    replica_optimizer.py init 32 feprest
    replica_optimizer.py optimize foo.log
    replica_optimizer.py update-topology foo_pp.top bar%d.top
    replica_optimizer.py update-mdp foo.mdp bar%d.mdp
    """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('command', metavar='COMMAND', type=str,
                        help='init/optimize/update-topology')
    parser.add_argument('arguments', action='store', type=str, nargs='+',
                        help="arguments")
    parser.add_argument('--basedir', action='store', type=str, default=os.getcwd(),
                        help="Base directory to save state (default: current dir)")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    if args.command == 'init':
        [N, mode] = args.arguments
        N = int(N)
        do_init(N, mode, basedir)
    elif args.command == 'optimize':
        [logfile] = args.arguments
        do_optimize_step(logfile, basedir)
    elif args.command == 'update-topology':
        [top_pp, top_generate] = args.arguments
        do_update_topology(top_pp, top_generate, basedir)
    elif args.command == 'update-mdp':
        [mdp_template, mdp_generate] = args.arguments
        do_update_mdp(mdp_template, mdp_generate, basedir)
    else:
        print("Error: unknown command", file=sys.stderr)
        sys.exit(1)



