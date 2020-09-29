import copy, cctk, argparse, sys, re
import numpy as np

# usage: python generate_conformations.py molecule.gjf molecule

parser = argparse.ArgumentParser(prog="generate_conformations.py")
parser.add_argument("--procs", "-p", type=int, default=16, help="Number of processors to use.")
parser.add_argument('-c', '--constraints', nargs='+', type=int, help="Atom numbers to freeze.")
parser.add_argument('-l', '--logfile', type=str, default="crest.out", help="Logfile to write output to.")
parser.add_argument('--noncovalent', action="store_true", help="Sets parameters to ``noncovalent`` in crest.")
parser.add_argument('-C', '--current', action="store_true", help="Store intermediate files in current directory.")
parser.add_argument("file", help="Molecule to perform conformation search on.")
parser.add_argument("new_filename", help="Prefix for generated pairs.")

args = vars(parser.parse_args(sys.argv[1:]))
args["new_filename"] = re.sub(".gjf$", "", args["new_filename"])

f = cctk.GaussianFile.read_file(args["file"])
molecule = f.get_molecule()

tmp = True
if args["current"]:
    tmp = False

e = molecule.csearch(nprocs=args["procs"], constraints=args["constraints"], logfile=args["logfile"], noncovalent=args["noncovalent"], use_tempdir=tmp)

for i, m in enumerate(e.molecule_list()):
    f.write_file(f"{args['new_filename']}_c{i:03d}.gjf", m)

