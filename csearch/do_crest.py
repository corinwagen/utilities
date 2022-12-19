import yaml, cctk, sys, os

filename = sys.argv[1]

settings = dict()
with open(filename, "r+") as f:
    settings = yaml.safe_load(f)

mol = cctk.GaussianFile.read_file(settings["input_geom"]).get_molecule()

freeze_atoms = set()
if "constraints" in settings:
    for constraint_row in settings["constraints"].values():
        atom1, atom2, dist = constraint_row.split(" ")
        atom1 = int(atom1)
        atom2 = int(atom2)

        if dist != "auto":
            dist = float(dist)
            mol.set_distance(atom1, atom2, dist)

        freeze_atoms.add(atom1)
        freeze_atoms.add(atom2)

logfile = "crest.log"
if "logfile" in settings:
    logfile = settings["logfile"]

noncovalent = False
if "noncovalent" in settings:
    noncovalent = settings["noncovalent"]

output_route_card = "#p sp m062x/6-31g(d)"
if "output_route_card" in settings:
    output_route_card = settings["output_route_card"]

assert "directory" in settings
if not os.path.exists(settings["directory"]):
    os.mkdir(settings["directory"])
os.chdir(settings["directory"])
e = mol.csearch(constraints=list(freeze_atoms), nprocs=8, noncovalent=noncovalent, logfile=logfile, use_tempdir=False, gfn="ff", additional_flags="--squick")
os.chdir("..")

output_name = filename.replace(".yaml", "_crest")
for i, m in enumerate(e.molecule_list()):
    cctk.GaussianFile.write_molecule_to_file(f"{output_name}_c{i:05d}.gjf", m, route_card=output_route_card)
