#!/usr/env/bin python3
import os
from typing import Any, Dict, List

import numpy
import pywindow
import toml


def _get_PDB_line(
    atom_number: int,
    residue_name: str,
    coordinates: numpy.ndarray,
    radius: float,
):
    field_name = "HETATM".ljust(6)
    atom_number = f"{atom_number}".rjust(5)
    atom_name = "PS1".rjust(4)
    residue_name = residue_name.ljust(3)
    chain = "A".rjust(1)
    residue_number = str(1).rjust(4)
    x = f"{coordinates[0]:8.3f}".rjust(8)
    y = f"{coordinates[1]:8.3f}".rjust(8)
    z = f"{coordinates[2]:8.3f}".rjust(8)
    occupancy = "0.0".rjust(6)
    b = f"{radius:6.2f}".ljust(6)
    element = "P".rjust(12)

    return f"{field_name}{atom_number} {atom_name} {residue_name} {chain}{residue_number}    {x}{y}{z}{occupancy}{b}{element}\n"


def _pdb(
    molecule: str,
    windows: Dict[str, Any],
    volume: float,
    center_of_mass: float,
    radius: float,
    basedir: str = ".",
):
    # Base name
    basename = os.path.basename(molecule).strip(".pdb").strip(".xyz")

    # Write pdb file
    with open(os.path.join(basedir, f"{basename}.pywindow.pdb"), "w") as f:
        f.write("HEADER\n")
        f.write(f"HEADER 0  - {'Void volume'.ljust(34)}: {volume:6.2f}\n")
        for nwindow, window in windows.items():
            f.write(
                f"HEADER {nwindow}  - {f'Window {nwindow} diameter'.ljust(34)}: {window[3]*2:.4f}\n"
            )
        f.write(_get_PDB_line(0, "COM", center_of_mass, radius))
        for nwindow, window in windows.items():
            f.write(_get_PDB_line(nwindow, f"W{nwindow}", window[0:3], window[3]))


def _pymol(molecule: str, windows: Dict[str, Any], basedir: str = "."):
    # Base name
    basename = os.path.basename(molecule).strip(".pdb").strip(".xyz")

    # Write script to visualization
    with open(os.path.join(basedir, f"{basename}-pymol2.py"), "w") as f:
        f.write("import pymol\n")
        f.write("from pymol import cmd, stored\n\n")
        f.write('pymol.finish_launching(["pymol", "-q"])\n\n')
        f.write(f'cmd.load("../../../hosts/{basename}.pdb", quiet=False)\n')
        f.write(f'cmd.load("{basename}.pywindow.pdb", quiet=False)\n\n')
        f.write("stored.b = []\n")
        residues = ["COM"] + [f"W{nwindow+1}" for nwindow in range(len(windows))]
        f.write(
            f"cmd.iterate(\"(resn {'+'.join(residues)})\", \"stored.b.append(b)\")\n"
        )
        for i, residue in enumerate(residues):
            f.write(f'cmd.alter("resn {residue}", "vdw=stored.b[{i}]")\n')

        f.write("cmd.rebuild()\n\n")
        f.write(
            f'cmd.spectrum("b", "blue_white_red", "{basename}.pywindow", [0, max(stored.b)])\n'
        )
        f.write(
            f'cmd.ramp_new("radius", "{basename}.pywindow", [0, max(stored.b)], ["blue", "white", "red"])\n'
        )
        f.write("cmd.orient()\n\n")
        f.write(f'cmd.save("{f"{basename}.pse"}")\n')


def _post_processing(
    molecule: str, results: Dict[str, Any], basedir: str = "."
) -> None:
    # Base name
    basename = os.path.basename(molecule).strip(".pdb").strip(".xyz")

    # Get void diameter, volume and center of mass
    diameter = numpy.array(results["pore_diameter_opt"]["diameter"])
    radius = diameter / 2
    volume = numpy.array(results["pore_volume_opt"])
    center_of_mass = numpy.array(results["pore_diameter_opt"]["centre_of_mass"])

    # Get windows radius and center of mass
    nwindows = len(results["windows"]["diameters"])
    windows = {}
    for window in range(nwindows):
        # Center of mass
        windows[window + 1] = numpy.array(results["windows"]["centre_of_mass"][window])
        # Radius
        windows[window + 1] = numpy.append(
            windows[window + 1], float(results["windows"]["diameters"][window]) / 2
        )

    # Write PDB file
    _pdb(molecule, windows, volume, center_of_mass, radius, basedir)

    # Write PyMOL script to visualization
    _pymol(molecule, windows, basedir)

    # Write results
    with open(os.path.join(basedir, f"{basename}.results.toml"), "w") as f:
        f.write('title = "pywindow results file"\n\n')
        toml.dump(results, f)


def _run_pywindow(molecule: str, basedir: str = ".") -> None:
    # Base name
    basename = os.path.basename(molecule).strip(".pdb").strip(".xyz")

    # Basedir
    os.makedirs(os.path.join(basedir, basename), exist_ok=True)

    # Load molecular system
    molsys = pywindow.MolecularSystem.load_file(molecule)

    # Convert molecular system to molecule
    mol = molsys.system_to_molecule()

    # Perform full_analysis
    results = mol.full_analysis()

    # Post-processing
    results = _post_processing(molecule, results, os.path.join(basedir, basename))


def run(molecules: List[str]) -> None:
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

    # Create basedir
    basedir = "./results/pywindow"
    os.makedirs(basedir, exist_ok=True)

    # Instructions to visualize results
    with open(os.path.join(basedir, "INSTRUCTIONS.md"), "w") as f:
        f.write("# Instructions\n\n")
        f.write("To visualiaze pywindow results for any cage, run:\n\n")
        f.write("```bash\npython {ID}-pymol2.py\n```\n")
        f.write("So, for cage A1, inside results/pywindow/A1 directory, run:\n\n")
        f.write("```bash\npython A1-pymol2.py\n```\n")

    # Run pywindow
    print("> pywindow (v0.0.4)")
    for molecule in molecules:
        print(molecule)
        _run_pywindow(molecule, basedir)
