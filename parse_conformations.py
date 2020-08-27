import sys, re, glob, cctk, argparse, copy
import numpy as np
import pandas as pd
import multiprocessing as mp
from multiprocessing import freeze_support

from tabulate import tabulate
from tqdm import tqdm

# usage: python parse_conformations.py -c eneryg_cutoff prefix_of_new_files path/to/output/*.out

def read(file):
    if re.match("slurm", file):
        return
    try:
        output_file = cctk.GaussianFile.read_file(file)
        if isinstance(output_file, list):
            if len(output_file) > 0:
                output_file = output_file[-1]
            else:
                return
        return output_file
    except Exception as e:
        print(f"Error reading {file}\n{e}")
        return

def main():
    e = cctk.ConformationalEnsemble()

    parser = argparse.ArgumentParser(prog="parse_conformations.py")
    parser.add_argument("--cutoff", "-c", type=float, default=15)
    parser.add_argument("--rmsd_cutoff", "-C", type=float, default=0.5)
    parser.add_argument("prefix", type=str)
    parser.add_argument("files", nargs='+')
    args = vars(parser.parse_args(sys.argv[1:]))

    print("\n\033[3mreading files:\033[0m")

    pool = mp.Pool(processes=16)
    for output_file in tqdm(pool.imap(read, args['files']), total=len(args['files'])):
        molecule = output_file.get_molecule()
        # there has got to be a better way to do this
        e.add_molecule(*list(output_file.ensemble.items())[-1])
        e[molecule, "iters"] = len(output_file.ensemble)
        e[molecule, "success"] = output_file.successful_terminations
        e[molecule, "imaginary"] = output_file.imaginaries()

    if len(e) == 0:
        print("no jobs to analyze!")
        exit()

    print(f"{len(e)} files read.")
    e, rmsds = e.eliminate_redundant(RMSD_cutoff=args['rmsd_cutoff'], return_RMSD=True)
    print(f"{len(e)} distinct conformations identified (RMSD cutoff of {args['rmsd_cutoff']:.2f}).")

    print("\n\033[3manalysis:\033[0m")

    property_names = ["filename", "new_filename", "rmsd", "iters", "energy", "success", "imaginary"]
    values = e[:, property_names]
    if not isinstance(values[0], list):
        values = [values]

    df = pd.DataFrame(values, columns=property_names).fillna(0)
    df["rmsd"] = rmsds
    df["new_filename"] = [f"{args['prefix']}_c{i:03}.gjf" for i in range(len(df))]
    df["rel_energy"] = (df.energy - df.energy.min()) * 627.509469
    df.sort_values("rel_energy", inplace=True)

    normal_col_names = copy.deepcopy(df.columns)

    df.columns = [f"\033[1m{c}\033[0m" for c in df.columns]
    print(tabulate(df, headers="keys", tablefmt="presto", floatfmt=".5f"))
    df.columns = normal_col_names

    template_file = read(args['files'][0])

    print("\n\033[3moutput:\033[0m")
    for m, n, e in zip(e.molecule_list(), df["new_filename"], df["rel_energy"]):
        if e <= args["cutoff"]:
            template_file.write_file(n, molecule=m)
            print(f"Wrote {n}")
        else:
            print(f"Skipping {n} and all subsequent rows: relative energy of {e:.2f} exceeds cutoff of {args['cutoff']:.2f}")
            break

if __name__ == '__main__':
    freeze_support()
    main()

