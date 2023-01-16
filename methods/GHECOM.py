#!/usr/env/bin python3
import os
import subprocess
from typing import List, Union


def _pymol(molecule: str, pocket: str, gw: float, basedir: str = "."):
    # Base name
    basename = os.path.basename(molecule).strip(".pdb").strip('.pocketness')

    pocket = os.path.basename(pocket).strip(".pdb")
    molecule = os.path.basename(molecule).strip(".pdb")

    # Write script to visualization
    with open(
        os.path.join(basedir, f"{basename}-pymol2.py"), "w"
    ) as f:
        f.write("import pymol\n")
        f.write("from pymol import cmd, stored\n\n")
        f.write('pymol.finish_launching(["pymol", "-q"])\n\n')
        f.write(f'cmd.load("{molecule}.pdb", quiet=False)\n\n')
        f.write(f'cmd.load("{pocket}.pdb", quiet=False)\n')

        f.write(f'cmd.hide("sticks", "{pocket}")\n')
        f.write(f'cmd.show("spheres", "{pocket}")\n')
        f.write(f'cmd.alter("resn GRD", "vdw={gw/2}")\n')
        f.write(f'cmd.alter("resn CEN", "vdw=0.1")\n')
        f.write(f"cmd.rebuild()\n\n")
        f.write("stored.b = []\n")
        f.write(f'cmd.iterate("(obj {pocket})", "stored.b.append(b)")\n')
        f.write(
            f'cmd.spectrum("b", "blue_white_red", "{pocket}", [min(stored.b), max(stored.b)])\n'
        )
        f.write(
            f'cmd.ramp_new("Rinaccess", "{pocket}", [min(stored.b), max(stored.b)], ["blue", "white", "red"])\n\n'
        )
        f.write("stored.b = []\n")
        f.write(f'cmd.iterate("(obj {molecule})", "stored.b.append(b)")\n')
        f.write(
            f'cmd.spectrum("b", "blue_white_red", "{molecule}", [min(stored.b), max(stored.b)])\n'
        )
        f.write(
            f'cmd.ramp_new("Pocketness", "{molecule}", [min(stored.b), max(stored.b)], ["blue", "white", "red"])\n\n'
        )
        f.write("cmd.orient()\n\n")
        f.write(f'cmd.save("{f"{basename}.pse"}")\n')


def _run_GHECOM(
    molecule: str,
    gw: float = 0.8,
    rlx: float = 10.0,
    basedir: str = ".",
) -> None:
    # Basename
    basename = os.path.basename(molecule).strip(".pdb")

    # Basedir
    basedir = os.path.join(basedir, basename)
    os.makedirs(basedir, exist_ok=True)

    # Run GHECOM
    command = f"ghecom -M M -ipdb {molecule} -opocpdb {os.path.join(basedir, f'{basename}.pocket.pdb')} -opdb {os.path.join(basedir, f'{basename}.pocketness.pdb')} -ores {os.path.join(basedir, f'{basename}.res')} -atmhet B -gw {gw} -rlx {rlx}"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    # Write pymol visualization
    _pymol(
        os.path.join(basedir, f"{basename}.pocketness.pdb"),
        os.path.join(basedir, f"{basename}.pocket.pdb"),
        gw,
        basedir,
    )


def run(
    molecules: List[str],
    gws: Union[float, List[float]] = 0.8,
    rlxs: Union[float, List[float]] = 10.0,
) -> None:
    # Check arguments
    if type(molecules) not in [list]:
        raise TypeError("`molecules` must be a list of PDB files.")
    if any([not molecule.endswith(".pdb") for molecule in molecules]):
        raise ValueError("`molecules` must contain only PDB files.")
    if any([not os.path.isfile(molecule) for molecule in molecules]):
        raise ValueError("`molecules` must contain valid PDB files.")
    if type(gws) not in [float, list]:
        raise TypeError("`gws` must be a float or a list of them.")
    elif type(gws) is float:
        gws = [gws] * len(molecules)
    elif len(gws) != len(molecules):
        raise Exception("`gws` must have the same length as `molecules`.")
    if type(rlxs) not in [float, list]:
        raise TypeError("`rlxs` must be a float or a list of them.")
    elif type(rlxs) is float:
        rlxs = [rlxs] * len(molecules)
    elif len(rlxs) != len(molecules):
        raise Exception("`rlxs` must have the same length as `molecules`.")

    # Basedir
    basedir = "./results/GHECOM"
    os.makedirs(basedir, exist_ok=True)

    # Instructions to visualize results
    with open(os.path.join(basedir, "INSTRUCTIONS.md"), "w") as f:
        f.write("# Instructions\n\n")
        f.write("To visualiaze GHECOM results for any cage, run:\n\n")
        f.write("```bash\npython {ID}-pymol2.py\n```\n")
        f.write("So, for cage A1, inside results/GHECOM/A1 directory, run:\n\n")
        f.write("```bash\npython A1-pymol2.py\n```\n")

    # Run GHECOM
    print("> GHECOM (21/07/2020)")
    for molecule, gw, rlx in zip(molecules, gws, rlxs):
        print(molecule)
        _run_GHECOM(molecule, gw, rlx, basedir)
