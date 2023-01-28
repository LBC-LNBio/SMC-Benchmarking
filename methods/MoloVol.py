#!/usr/env/bin python3

import os
import subprocess
from typing import List, Union


def _pymol(molecule: str, large_probe: float, basedir: str = "."):
    # Base name
    basename = os.path.basename(molecule).strip(".pdb").strip(".xyz")

    # Write script to visualization
    with open(os.path.join(basedir, f"{basename}-pymol2.py"), "w") as f:
        f.write("import pymol\n")
        f.write("from pymol import cmd, stored\n\n")
        f.write('pymol.finish_launching(["pymol", "-q"])\n\n')
        f.write(f'cmd.load("../../../hosts/{basename}.pdb", quiet=False)\n')
        f.write(
            f'cmd.load("{basename}_surface-map_grid-0.6_rad1-1.4_rad2-{large_probe:.1f}_full-structure.dx", quiet=False)\n\n'
        )
        f.write(f'cmd.isomesh("cavities", "{basename}_surface-map_grid-0.6_rad1-1.4_rad2-{large_probe:.1f}_full-structure", level=3.0) # Isolated cavity\n')
        f.write("cmd.orient()\n\n")
        f.write(f'cmd.save("{f"{basename}.pse"}")\n')

def _run_MoloVol(
    molecule: str,
    grid_spacing: float,
    small_probe: float,
    large_probe: float,
    basedir: str = ".",
):
    # Basename
    basename = os.path.basename(molecule).strip(".pdb")

    # Basedir
    basedir = os.path.join(basedir, basename)
    os.makedirs(basedir, exist_ok=True)

    # Run MoloVol
    command = f"molovol -fs {molecule} -g {grid_spacing} -r {small_probe} -r2 {large_probe} -do {basedir} -xt -xr -xc -ht"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    # pymol
    _pymol(molecule, large_probe, basedir)


def run(
    molecules: List[str],
    grid_spacings: Union[float, List[float]] = 0.2,
    small_probes: Union[float, List[float]] = 1.4,
    large_probes: Union[float, List[float]] = 10.0,
):
    # Check arguments
    if type(molecules) not in [list]:
        raise TypeError("`molecules` must be a list of PDB or XYZ files.")
    if any(
        [
            not (molecule.endswith(".pdb") or molecule.endswith(".xyz"))
            for molecule in molecules
        ]
    ):
        raise ValueError("`molecules` must contain only PDB or XYZ files.")
    if any([not os.path.isfile(molecule) for molecule in molecules]):
        raise ValueError("`molecules` must contain valid PDB or XYZ files.")
    if type(grid_spacings) not in [float, int, list]:
        raise TypeError("`grid_spacings` must be a float or a list of them.")
    elif type(grid_spacings) is float:
        grid_spacings = [grid_spacings] * len(molecules)
    elif len(grid_spacings) != len(molecules):
        raise Exception("`grid_spacings` must have the same length as `molecules`.")
    if type(small_probes) not in [float, int, list]:
        raise TypeError("`small_probes` must be a float or a list of them.")
    elif type(small_probes) is float:
        small_probes = [small_probes] * len(molecules)
    elif len(small_probes) != len(molecules):
        raise Exception("`small_probes` must have the same length as `molecules`.")
    if type(large_probes) not in [float, int, list]:
        raise TypeError("`large_probes` must be a float or a list of them.")
    elif type(large_probes) is float:
        large_probes = [large_probes] * len(molecules)
    elif len(large_probes) != len(molecules):
        raise Exception("`large_probes` must have the same length as `molecules`.")

    # Create basedir
    basedir = "./results/MoloVol"
    os.makedirs(basedir, exist_ok=True)

    # Instructions to visualize results
    with open(os.path.join(basedir, "INSTRUCTIONS.md"), "w") as f:
        f.write("# Instructions\n\n")
        f.write("To visualiaze MoloVol results for any cage, run:\n\n")
        f.write("```bash\npython {ID}-pymol2.py\n```\n")
        f.write("So, for cage A1, inside results/MoloVol/A1 directory, run:\n\n")
        f.write("```bash\npython A1-pymol2.py\n```\n")

    # Run MoloVol for all files
    print("> MoloVol (v1.0.0)")
    for molecule, gs, sp, lp in zip(
        molecules, grid_spacings, small_probes, large_probes
    ):
        print(molecule)
        _run_MoloVol(molecule, gs, sp, lp, basedir)
