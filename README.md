# utilities
_useful scripts for computational chemistry_

Most scripts will probably require [_cctk_](https://cctk.readthedocs.io/en/latest/) and a working Python 3.7+ installation.

### Contents:

#### generate_ion_pairs.py

Given two input molecules, randomly arranges them in relation to each other. Allows for unbiased generation of different conformations of loosely-bound outer-sphere complexes (e.g. ion pairs).

#### generate_conformations.py

Given an input structure, uses [``crest``](https://xtb-docs.readthedocs.io/en/latest/crest.html) to generate new conformations and outputs Gaussian input files.

#### parse_conformations.py

Given a list of output files, automatically sorts them by energy, removes duplicate conformations, and prints the low-energy conformations to new input files.

_Corin Wagen_
