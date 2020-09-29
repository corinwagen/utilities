import cctk, argparse, sys

parser = argparse.ArgumentParser(prog="submit_from_crest.py")
parser.add_argument('-c', '--constraints', nargs='+', type=int, help="Atom numbers to freeze.")
parser.add_argument("--prefix", "-p", type=str, help="Prefix for generated files.")
parser.add_argument("conformer_file")

args = vars(parser.parse_args(sys.argv[1:]))

confs = cctk.XYZFile.read_ensemble(args["conformer_file"])
print(f"{len(confs)} conformers read.")

footer = ""
for x in args["constraints"]:
    for y in args["constraints"]:
        if int(x) <= int(y):
            continue
        footer += f"B {int(x)} {int(y)} F\n"

footer += """
--Link1--
%nprocshared=16
%mem=32GB
%chk=opt.chk
#p opt=(ts,calcfc,noeigentest,maxstep=10) freq=noraman b3lyp/6-31g(d) empiricaldispersion=gd3bj geom=allcheck guess=read

"""

for idx, mol in enumerate(confs.molecules):
    cctk.GaussianFile.write_molecule_to_file(
        filename=f"{args['prefix']}_c{idx:03d}.gjf",
        molecule=mol,
        link0={"nprocshared": 16, "mem": "32GB", "chk": "opt.chk"},
        route_card="#p opt=modredundant b3lyp/6-31g(d) empiricaldispersion=gd3bj",
        footer=footer,
    )
