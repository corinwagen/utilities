import copy, cctk, argparse, sys, re
import numpy as np

# usage: python generate_conformations.py molecule.gjf molecule

parser = argparse.ArgumentParser(prog="generate_conformations.py")
parser.add_argument("--procs", "-p", type=int, default=16, help="Number of processors to use.")
parser.add_argument('-c', '--constraints', nargs='+', type=int, help="Atom numbers to freeze.")
parser.add_argument("file", help="Molecule to perform conformation search on.")
parser.add_argument("new_filename", help="Prefix for generated pairs.")

args = vars(parser.parse_args(sys.argv[1:]))
args["new_filename"] = re.sub(".gjf$", "", args["new_filename"])

print("freezing")
print(args["constraints"])

f = cctk.GaussianFile.read_file(args["file"])
molecule = f.get_molecule()

e = molecule.csearch(nprocs=args["procs"], constraints=args["constraints"])

for i, m in enumerate(e.molecule_list()):
    f.write_file(f"{args['new_filename']}_c{i:03d}.gjf", m)

