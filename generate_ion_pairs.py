import copy, cctk, argparse, sys, re
import numpy as np

# usage: python generate_ion_pairs.py -n 100 -r 10 cation.out anion.out cation_anion

parser = argparse.ArgumentParser(prog="generate_ion_pairs.py")
parser.add_argument("--num", "-n", type=int, default=25, help="Number of conformations to generate.")
parser.add_argument("--radius", "-r", type=float, default=8, help="Radius of imaginary sphere to distribute pairs on.")
parser.add_argument("--optimize", "-o", action="store_true", help="Whether or not to pre-optimize pairs with ``xtb``.")
parser.add_argument("file1", help="Molecule to place in the center (path to Gaussian .gjf or .out file).")
parser.add_argument("file2", help="Molecule to place around the outside of the imaginary sphere (path to Gaussian .gjf or .out file).")
parser.add_argument("new_filename", help="Prefix for generated pairs.")

args = vars(parser.parse_args(sys.argv[1:]))

def spherical_random(radius=1):
    """
    Generates a random point on a sphere of radius ``radius``.
    """
    v = np.random.normal(size=3)
    v = v/np.linalg.norm(v)
    return v * radius

f1 = cctk.GaussianFile.read_file(args["file1"])
m1 = f1.get_molecule().center()
m2 = cctk.GaussianFile.read_file(args["file2"]).get_molecule().center()

args["new_filename"] = re.sub(".gjf$", "", args["new_filename"])

for i in range(args["num"]):
    m = copy.deepcopy(m1)
    x = copy.deepcopy(m2)

    x.translate_molecule(spherical_random(radius=args["radius"]))
    x.rotate_molecule(np.array([1,0,0]), np.random.random()*360)
    x.rotate_molecule(np.array([0,1,0]), np.random.random()*360)
    x.rotate_molecule(np.array([0,0,1]), np.random.random()*360)

    mx = cctk.Molecule.combine_molecules(m, x)

    if args["optimize"]:
        mx.optimize()

    f1.write_file(f"{args['new_filename']}_c{i:03d}.gjf", mx)
