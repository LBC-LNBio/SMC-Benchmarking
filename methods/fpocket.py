#!/usr/bin/env python3
import os
from typing import List, Union
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_select import fpocket_select


def _pymol(molecule: str, cavity: str, atm: str, basedir: "."):
    # Base name
    basename = os.path.basename(molecule).strip(".pdb")

    molecule = os.path.realpath(molecule)
    cavity = os.path.realpath(cavity)
    atm = os.path.realpath(atm)

    # Write script to visualization
    with open(os.path.join(basedir, f"{basename}-pymol2.py"), "w") as f:
        f.write("import pymol\n")
        f.write("from pymol import cmd, stored\n\n")
        f.write('pymol.finish_launching(["pymol", "-q"])\n\n')
        f.write(f'cmd.load("../../../hosts/{basename}.pdb", quiet=False)\n')
        f.write(f'cmd.load("{os.path.basename(cavity)}", quiet=False)\n')
        f.write(f'cmd.load("{os.path.basename(atm)}", quiet=False)\n')
        f.write("\n")
        f.write(
            f'cmd.hide("sticks", "{os.path.basename(cavity).replace(".pqr", "")}")\n'
        )
        f.write(
            f'cmd.show("surface", "{os.path.basename(cavity).replace(".pqr", "")}")\n'
        )
        f.write(
            f'cmd.color("white", "{os.path.basename(cavity).replace(".pqr", "")}")\n'
        )
        f.write(f"cmd.rebuild()\n\n")
        f.write("cmd.orient()\n\n")
        f.write(f'cmd.save("{f"{basename}.pse"}")\n')


def _run_fpocket(
    pdb: str,
    min_radius: float = 3.4,
    max_radius: float = 6.2,
    num_spheres: int = 15,
    selection: int = 1,
    basedir: str = ".",
) -> None:
    # Base Name
    basename = os.path.basename(pdb).replace(".pdb", "")

    # Basedir
    basedir = os.path.join(basedir, basename)
    os.makedirs(basedir, exist_ok=True)

    # Output file
    output = os.path.join(basedir, f"{basename}.zip")

    # Run fpocket
    if pdb.find("O2.pdb") < 0:
        fpocket_run(
            input_pdb_path=pdb,
            output_pockets_zip=output,
            output_summary=os.path.join(basedir, f"summary.json"),
            properties={
                "min_radius": min_radius,
                "max_radius": max_radius,
                "num_spheres": num_spheres,
            },
        )

        # Select
        fpocket_select(
            input_pockets_zip=output,
            output_pocket_pdb=os.path.join(basedir, f"pocket{selection}_atm.pdb"),
            output_pocket_pqr=os.path.join(basedir, f"pocket{selection}_vert.pqr"),
            properties={"pocket": selection},
        )

        # Write pymol visualization
        _pymol(
            pdb,
            os.path.join(basedir, f"pocket{selection}_vert.pqr"),
            os.path.join(basedir, f"pocket{selection}_atm.pdb"),
            basedir,
        )


def run(
    molecules: List[str],
    min_radius: Union[float, List[float]] = 3.2,
    max_radius: Union[float, List[float]] = 6.4,
    num_spheres: Union[int, List[int]] = 15,
    selection: Union[int, List[int]] = 1,
) -> None:
    # Check arguments
    if type(molecules) not in [list]:
        raise TypeError("`molecules` must be a list of PDB files.")
    if any([not molecule.endswith(".pdb") for molecule in molecules]):
        raise ValueError("`molecules` must contain only PDB files.")
    if any([not os.path.isfile(molecule) for molecule in molecules]):
        raise ValueError("`molecules` must contain valid PDB files.")
    if type(min_radius) not in [float, list]:
        raise TypeError("`min_radius` must be a float or a list of them.")
    elif type(min_radius) is float:
        min_radius = [min_radius] * len(molecules)
    elif len(min_radius) != len(molecules):
        raise Exception("`min_radius` must have the same length as `molecules`.")
    if type(max_radius) not in [float, list]:
        raise TypeError("`max_radius` must be a float or a list of them.")
    elif type(max_radius) is float:
        max_radius = [max_radius] * len(molecules)
    elif len(max_radius) != len(molecules):
        raise Exception("`max_radius` must have the same length as `molecules`.")
    if type(num_spheres) not in [int, list]:
        raise TypeError("`num_spheres` must be a integer or a list of them.")
    elif type(num_spheres) is int:
        num_spheres = [num_spheres] * len(molecules)
    elif len(num_spheres) != len(molecules):
        raise Exception("`num_spheres` must have the same length as `molecules`.")
    if type(selection) not in [int, list]:
        raise TypeError("`selection` must be a integer or a list of them.")
    elif type(selection) is int:
        selection = [selection] * len(molecules)
    elif len(selection) != len(molecules):
        raise Exception("`selection` must have the same length as `molecules`.")

    # Create basedir
    basedir = "./results/fpocket"
    os.makedirs(basedir, exist_ok=True)

    # Instructions to visualize results
    with open(os.path.join(basedir, "INSTRUCTIONS.md"), "w") as f:
        f.write("# Instructions\n\n")
        f.write("To visualiaze fpocket results for any cage, run:\n\n")
        f.write("```bash\npython {ID}-pymol2.py\n```\n")
        f.write("So, for cage A1, inside results/fpocket/A1 directory, run:\n\n")
        f.write("```bash\npython A1-pymol2.py\n```\n")

    # Run fpocket for all files
    print("> fpocket (v3.1.4.2)")
    for molecule, m, M, n, s in zip(
        molecules, min_radius, max_radius, num_spheres, selection
    ):
        print(molecule)
        _run_fpocket(molecule, m, M, n, s, basedir)

    # Remove log and error files
    for fn in os.listdir():
        if fn.endswith(".out") or fn.endswith(".err"):
            os.remove(fn)
