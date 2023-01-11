import numpy
import mdtraj
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Output index file for analysis")

    parser.add_argument("-i", dest="input", type=str, required=True, help="Input structure file")
    parser.add_argument("-o", dest="output", type=str, required=True, help="Output index file")

    return parser.parse_args()

def gen_file(instr, outndx):
    structure = mdtraj.load(instr)

    com = mdtraj.compute_center_of_mass(structure)[0, :]

    displs = structure.xyz[0, :, :] - com[numpy.newaxis, :]
    dist2 = numpy.sum(displs[:, :] * displs[:, :], axis=1)
    closest = numpy.argmin(dist2)

    with open(outndx, "w") as ofh:
        print("[ centering ]", file=ofh)
        print("%d" % (closest + 1), file=ofh)
        print(file=ofh)

        print("[ output ]", file=ofh)
        for i in range(structure.n_atoms):
            print("%d" % (i + 1), file=ofh)
        print(file=ofh)

if __name__ == "__main__":
    args = parse_args()

    gen_file(args.input, args.output)

