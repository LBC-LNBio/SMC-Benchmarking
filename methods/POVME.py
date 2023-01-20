#!/usr/env/bin python3
import os
import subprocess
from typing import List, Union


def _pymol(
    molecule: str,
    inclusion: str,
    cavity: str,
    surface: str,
    grid_spacing: float,
    basedir: str = ".",
):
    # Base name
    basename = os.path.basename(molecule).strip(".pdb")

    # Write script to visualization
    with open(os.path.join(basedir, f"{basename}-pymol2.py"), "w") as f:
        f.write("import pymol\n")
        f.write("from pymol import cmd, stored\n\n")
        f.write('pymol.finish_launching(["pymol", "-q"])\n\n')
        f.write(f'cmd.load("../../../hosts/{basename}.pdb", quiet=False)\n')
        f.write(f'cmd.load("{inclusion}", quiet=False)\n')
        f.write(f'cmd.load("{cavity}", quiet=False)\n')
        f.write(f'cmd.load("{surface}", quiet=False)\n')
        f.write("\n")
        f.write(f'cmd.show("dots", "{os.path.basename(cavity).strip(".dx")}")\n\n')
        f.write(
            f'cmd.alter("obj {os.path.basename(inclusion).strip(".pdb")}", "vdw={grid_spacing/2}")\n'
        )
        f.write(
            f'cmd.alter("obj {os.path.basename(surface).strip(".pdb")}", "vdw={grid_spacing/2}")\n'
        )
        f.write(f"cmd.rebuild()\n\n")
        f.write("cmd.orient()\n\n")
        f.write(f'cmd.save("{f"{basename}.pse"}")\n')


def _input_file(
    molecule: str,
    inclusion_region: List[float],
    grid_spacing: float,
    contiguous_points_criteria: int,
    basedir: str,
) -> None:
    # Basename
    basename = os.path.basename(molecule).strip(".pdb")

    with open(os.path.join(basedir, f"{basename}.ini"), "w") as f:
        f.write(f"# POVME 3.0 {molecule} Input File\n\n")
        f.write(f"PDBFileName {os.path.realpath(molecule)}\n")
        f.write(f"GridSpacing {grid_spacing:.2f}\n")
        if len(inclusion_region) == 4:  # Inclusion Sphere
            f.write(
                f"InclusionSphere {inclusion_region[0]} {inclusion_region[1]} {inclusion_region[2]} {inclusion_region[3]}\n"
            )
        elif len(inclusion_region) == 6:  # Inclusion Box
            f.write(
                f"InclusionBox {inclusion_region[0]} {inclusion_region[1]} {inclusion_region[2]} {inclusion_region[3]} {inclusion_region[4]} {inclusion_region[5]}\n"
            )
        f.write("DistanceCutoff 1.09\n")
        f.write("ConvexHullExclusion max\n")
        f.write(f"ContiguousPointsCriteria {contiguous_points_criteria}\n")
        f.write("NumProcessors 12\n")
        f.write(
            f"OutputFilenamePrefix {os.path.join(os.path.realpath(basedir), f'{basename}_')}\n"
        )


def _run_POVME(
    molecule: str,
    inclusion_region: List[float],
    grid_spacing: float = 1.0,
    contiguous_points_criteria: int = 3,
    basedir: str = ".",
) -> None:
    # Basename
    basename = os.path.basename(molecule).strip(".pdb")

    # Basedir
    basedir = os.path.join(basedir, basename)
    os.makedirs(basedir, exist_ok=True)

    # Prepare input file
    _input_file(
        molecule, inclusion_region, grid_spacing, contiguous_points_criteria, basedir
    )

    # Run POVME3
    command = f"POVME3.py {os.path.join(basedir, basename)}.ini"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    # Write pymol visualization
    inclusion = os.path.join(f"{basename}_frameInfo", f"{basename}_inclusion.pdb")
    cavity = os.path.join(f"{basename}_frameInfo", f"{basename}_volumetric_density.dx")
    surface = os.path.join(f"{basename}_frameInfo", f"{basename}_frame_1_surface.pdb")
    _pymol(molecule, inclusion, cavity, surface, grid_spacing, basedir)


def run(
    molecules: List[str],
    inclusion_regions: List[List[float]],
    grid_spacings: Union[float, List[float]] = 1.0,
    contiguous_points_criterias: Union[int, List[int]] = 3,
) -> None:
    # Check arguments
    if type(molecules) not in [list]:
        raise TypeError("`molecules` must be a list of PDB files.")
    if any([not molecule.endswith(".pdb") for molecule in molecules]):
        raise ValueError("`molecules` must contain only PDB files.")
    if any([not os.path.isfile(molecule) for molecule in molecules]):
        raise ValueError("`molecules` must contain valid PDB files.")
    if type(inclusion_regions) not in [list]:
        raise TypeError("`inclusion_region` must be a list of coordinates.")
    elif len(inclusion_regions) != len(molecules):
        raise Exception("`inclusion_regions` must have the same length as `molecules`.")
    if type(grid_spacings) not in [float, int, list]:
        raise TypeError("`grid_spacings` must be a float or a list of them.")
    elif type(grid_spacings) is float:
        grid_spacings = [grid_spacings] * len(molecules)
    elif len(grid_spacings) != len(molecules):
        raise Exception("`grid_spacings` must have the same length as `molecules`.")
    if type(contiguous_points_criterias) not in [float, int, list]:
        raise TypeError(
            "`contiguous_points_criterias` must be an integer or a list of them."
        )
    elif type(contiguous_points_criterias) is (int or float):
        contiguous_points_criterias = [contiguous_points_criterias] * len(molecules)
    elif len(contiguous_points_criterias) != len(contiguous_points_criterias):
        raise Exception(
            "`contiguous_points_criterias` must have the same length as `molecules`."
        )

    # Basedir
    basedir = "./results/POVME"
    os.makedirs(basedir, exist_ok=True)

    # Instructions to visualize results
    with open(os.path.join(basedir, "INSTRUCTIONS.md"), "w") as f:
        f.write("# Instructions\n\n")
        f.write("To visualiaze POVME results for any cage, run:\n\n")
        f.write("```bash\npython {ID}-pymol2.py\n```\n")
        f.write("So, for cage A1, inside results/POVME/A1 directory, run:\n\n")
        f.write("```bash\npython A1-pymol2.py\n```\n")

    # Run POVME
    print("> POVME (v3.0.35)")
    for molecule, inclusion_region, grid_spacing, contiguous_points_criteria in zip(
        molecules, inclusion_regions, grid_spacings, contiguous_points_criterias
    ):
        print(molecule)
        _run_POVME(
            molecule,
            inclusion_region,
            grid_spacing,
            contiguous_points_criteria,
            basedir,
        )
