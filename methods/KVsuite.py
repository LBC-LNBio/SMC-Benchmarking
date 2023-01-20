#!/usr/env/bin python3
import os
import subprocess
from typing import Any, List, Optional, Union

import pyKVFinder
import toml


def _split(parameter: Union[List[float], float]) -> List[float]:
    if type(parameter) in [list]:
        p1, p2 = parameter
    else:
        p1 = p2 = parameter

    return p1, p2


def _pymol(
    molecule: str, cavity: str, surface: Optional[str], step: float, basedir: str = "."
):
    # Base name
    basename = os.path.basename(molecule).strip(".pdb")

    # Write script to visualization
    with open(os.path.join(basedir, f"{basename}-pymol2.py"), "w") as f:
        f.write("import pymol\n")
        f.write("from pymol import cmd, stored\n\n")
        f.write('pymol.finish_launching(["pymol", "-q"])\n\n')
        f.write(f'cmd.load("{molecule}", quiet=False)\n')
        f.write(f'cmd.load("{cavity}", quiet=False)\n')
        if surface is not None:
            f.write(f'cmd.load("{surface}", quiet=False)\n')
        f.write("\n")
        f.write(f'cmd.hide("sticks", "{os.path.basename(cavity).strip(".pdb")}")\n')
        f.write(f'cmd.show("spheres", "{os.path.basename(cavity).strip(".pdb")}")\n')
        f.write(
            f'cmd.alter("obj {os.path.basename(cavity).strip(".pdb")}", "vdw={step/2}")\n'
        )
        if surface is not None:
            f.write(
                f'cmd.alter("obj {os.path.basename(surface).strip(".pdb")}", "vdw={step/2}")\n'
            )
        f.write(f"cmd.rebuild()\n\n")
        f.write("stored.b = []\n")
        f.write(
            f'cmd.iterate("(obj {os.path.basename(cavity).strip(".pdb")})", "stored.b.append(b)")\n'
        )
        f.write(
            f'cmd.spectrum("b", "blue_white_red", "{os.path.basename(cavity).strip(".pdb")}", [min(stored.b), max(stored.b)])\n'
        )
        f.write(
            f'cmd.ramp_new("Depth", "{os.path.basename(cavity).strip(".pdb")}", [min(stored.b), max(stored.b)], ["blue", "white", "red"])\n\n'
        )
        f.write("cmd.orient()\n\n")
        f.write(f'cmd.save("{f"{basename}.pse"}")\n')


def _run_pyKVFinder(
    molecule: str,
    step: float = 0.6,
    probe_out: float = 8.0,
    removal_distance: float = 2.4,
    volume_cutoff: float = 5.0,
    basedir: str = ".",
) -> None:
    # Base Name
    basename = os.path.basename(molecule).strip(".pdb")

    # Basedir
    basedir = os.path.join(basedir, basename)
    os.makedirs(basedir, exist_ok=True)

    # Atomic information
    atomic = pyKVFinder.read_pdb(molecule)

    # Grid dimensions
    vertices = pyKVFinder.get_vertices(atomic, probe_out=probe_out)

    # Cavity detection
    _, cavities = pyKVFinder.detect(
        atomic,
        vertices,
        step=step,
        probe_out=probe_out,
        removal_distance=removal_distance,
        volume_cutoff=volume_cutoff,
    )

    # Spatial characterization
    surface, volume, area = pyKVFinder.spatial(cavities, step=step)

    # Depth characterization
    depths, max_depth, avg_depth = pyKVFinder.depth(cavities, step=step)

    # Export cavities
    output = os.path.join(basedir, f"{basename}.KVFinder.output.pdb")
    pyKVFinder.export(output, cavities, surface, vertices, B=depths, step=step)
    osurf = os.path.join(basedir, f"{basename}.surface.pdb")
    surface_representation = (cavities == 0).astype(int) * 2
    pyKVFinder.export(
        osurf, surface_representation * 2, None, vertices, B=depths, step=step
    )

    # Write results
    oresults = os.path.join(basedir, f"{basename}.KVFinder.results.toml")
    pyKVFinder.write_results(
        oresults,
        molecule,
        None,
        output,
        volume=volume,
        area=area,
        max_depth=max_depth,
        avg_depth=avg_depth,
        step=step,
    )

    # Write pymol visualization
    _pymol(
        f"../../../{molecule}",
        os.path.basename(output),
        os.path.basename(osurf),
        step,
        basedir,
    )


def _run_parKVFinder(
    molecule: str,
    step: float = 0.6,
    probe_out: float = 8.0,
    removal_distance: float = 2.4,
    volume_cutoff: float = 5.0,
    basedir: str = ".",
) -> None:
    # Base Name
    basename = os.path.basename(molecule).strip(".pdb")

    # Basedir
    basedir = os.path.join(basedir, basename)
    os.makedirs(basedir, exist_ok=True)

    # Check if parameters exist
    if os.path.exists("parameters.toml"):
        os.remove("parameters.toml")

    # Write parameters file
    with open("parameters.toml", "w") as f:
        f.write("# TOML configuration file for parKVFinder software.\n")
        f.write('\ntitle = "parKVFinder parameters file"\n')
        f.write("\n[FILES_PATH]\n")
        f.write("# The path of van der Waals radii dictionary for parKVFinder.\n")
        f.write(
            f"dictionary = \"{os.path.join(os.getenv('KVFinder_PATH'), 'dictionary')}\"\n"
        )
        f.write("# The path of the input PDB file.\n")
        f.write(f'pdb = "{os.path.realpath(molecule)}"\n')
        f.write("# The path of the output directory.\n")
        f.write(f'output = "{os.path.realpath(basedir)}"\n')
        f.write("# Base name for output files.\n")
        f.write(f"base_name = \"{os.path.basename(molecule).strip('.pdb')}\"\n")
        f.write("# Path for the ligand's PDB file.\n")
        f.write(f'ligand = "-"\n')

        f.write("\n[SETTINGS]\n")
        f.write("# Settings for parKVFinder software.\n")
        f.write("\n[SETTINGS.modes]\n")
        f.write("# Whole protein mode defines the search space as the whole protein.\n")
        f.write(f"whole_protein_mode = true\n")
        f.write(
            "# Box adjustment mode defines the search space as the box drawn in PyMOL.\n"
        )
        f.write(f"box_mode = false\n")
        f.write(
            "# Resolution mode implicitly sets the step size (grid spacing) of the 3D grid.\n"
        )
        f.write(
            "# If set to High, sets a voxel volume of 0.2. If set to Medium, sets a voxel volume of 0.1. If set to Low, it sets a voxel volume of 0.01. If set to Off, the step size must be set explicitly.\n"
        )
        f.write(f'resolution_mode = "Off"\n')
        f.write(
            "# Surface mode defines the type of surface representation to be applied, van der Waals molecular surface (true) or solvent accessible surface (false).\n"
        )
        f.write(f"surface_mode = true\n")
        f.write(
            "# Cavity representation defines whether cavities are exported to the output PDB file as filled cavities (true) or filtered cavities (false).\n"
        )
        f.write(f"kvp_mode = true\n")
        f.write(
            "# Ligand adjustment mode defines the search space around the ligand.\n"
        )
        f.write(f"ligand_mode = false\n")

        f.write("\n[SETTINGS.step_size]\n")
        f.write(
            "# Sets the 3D grid spacing. It directly affects accuracy and runtime.\n"
        )
        f.write(f"step_size = {step:.2f}\n")

        f.write("\n[SETTINGS.probes]\n")
        f.write(
            "# parKVFinder works with a dual probe system. A smaller probe, called Probe In, and a bigger one, called Probe Out, rolls around the protein.\n"
        )
        f.write(
            "# Points reached by the Probe In, but not the Probe Out are considered cavity points.\n"
        )
        f.write("# Sets the Probe In diameter. Default: 1.4 angstroms.\n")
        f.write(f"probe_in = 1.40\n")
        f.write("# Sets the Probe Out diameter. Default: 4.0 angstroms.\n")
        f.write(f"probe_out = {probe_out:.2f}\n")

        f.write("\n[SETTINGS.cutoffs]\n")
        f.write(
            "# Sets a volume cutoff for the detected cavities. Default: 5.0 angstroms.\n"
        )
        f.write(f"volume_cutoff = {volume_cutoff:.2f}\n")
        f.write(
            "# Sets a distance cutoff for a search space around the ligand in ligand adjustment mode. Default: 5.0 angstroms.\n"
        )
        f.write(f"ligand_cutoff = 5.0\n")
        f.write(
            "# Sets a removal distance for the cavity frontier, which is defined by comparing Probe In and Probe Out surfaces. Default: 2.4 angstroms.\n"
        )
        f.write(f"removal_distance = {removal_distance:.2f}\n")

        f.write("\n[SETTINGS.visiblebox]\n")
        f.write(
            "# Coordinates of the vertices that define the visible 3D grid. Only four points are required to define the search space.\n\n"
        )
        visiblebox = {
            "p1": {"x": 0.0, "y": 0.0, "z": 0.0},
            "p2": {"x": 0.0, "y": 0.0, "z": 0.0},
            "p3": {"x": 0.0, "y": 0.0, "z": 0.0},
            "p4": {"x": 0.0, "y": 0.0, "z": 0.0},
        }
        internalbox = {
            "p1": {"x": -probe_out, "y": -probe_out, "z": -probe_out},
            "p2": {"x": probe_out, "y": -probe_out, "z": -probe_out},
            "p3": {"x": -probe_out, "y": probe_out, "z": -probe_out},
            "p4": {"x": -probe_out, "y": -probe_out, "z": probe_out},
        }
        d = {"SETTINGS": {"visiblebox": visiblebox, "internalbox": internalbox}}
        toml.dump(o=d, f=f)

    # Run parKVFinder
    command = "parKVFinder -p parameters.toml"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    # Write pymol visualization
    _pymol(
        f"../../../../../{molecule}",
        f"{basename}.KVFinder.output.pdb",
        None,
        step,
        os.path.join(basedir, "KV_Files", basename),
    )

    # Remove parameters.toml
    os.remove("parameters.toml")


def run(
    molecules: List[str],
    step: Union[float, List[Any]] = 0.6,
    probe_out: Union[float, List[Any]] = 8.0,
    removal_distance: Union[float, List[Any]] = 2.4,
    volume_cutoff: Union[float, List[Any]] = 5.0,
) -> None:
    # Check arguments
    if type(molecules) not in [list]:
        raise TypeError("`molecules` must be a list of PDB files.")
    if any([not molecule.endswith(".pdb") for molecule in molecules]):
        raise ValueError("`molecules` must contain only PDB files.")
    if any([not os.path.isfile(molecule) for molecule in molecules]):
        raise ValueError("`pdbs` must contain valid PDB files.")
    if type(probe_out) not in [float, int, list]:
        raise TypeError("`probe_out` must be a float or a list of them.")
    elif type(probe_out) is float:
        probe_out = [probe_out] * len(molecules)
    elif len(probe_out) != len(molecules):
        raise Exception("`probe_out` must have the same length as `molecules`.")
    if type(step) not in [float, int, list]:
        raise TypeError("`step` must be a float or a list of them.")
    elif type(step) is float:
        step = [step] * len(molecules)
    elif len(step) != len(molecules):
        raise Exception("`step` must have the same length as `molecules`.")
    if type(removal_distance) not in [float, int, list]:
        raise TypeError("`removal_distance` must be a float or a list of them.")
    elif type(removal_distance) is float:
        removal_distance = [removal_distance] * len(molecules)
    elif len(removal_distance) != len(molecules):
        raise Exception("`removal_distance` must have the same length as `molecules`.")
    if type(volume_cutoff) not in [float, int, list]:
        raise TypeError("`volume_cutoff` must be a float or a list of them.")
    elif type(volume_cutoff) is float:
        volume_cutoff = [volume_cutoff] * len(molecules)
    elif len(volume_cutoff) != len(molecules):
        raise Exception("`volume_cutoff` must have the same length as `molecules`.")

    # Create basedir
    basedir = "./results/KVFinder-suite"
    os.makedirs(basedir, exist_ok=True)

    # Instructions to visualize results
    with open(os.path.join(basedir, "INSTRUCTIONS.md"), "w") as f:
        f.write("# Instructions\n\n")
        f.write("## pyKVFinder\n\n")
        f.write("To visualiaze pyKVFinder results for any cage, run:\n\n")
        f.write("```bash\npython {ID}-pymol2.py\n```\n")
        f.write("So, for cage A1, inside results/KVFinder-suite/A1 directory, run:\n\n")
        f.write("```bash\npython A1-pymol2.py\n```\n\n")
        f.write("## parKVFinder\n\n")
        f.write("To visualiaze parKVFinder results for any cage, run:\n\n")
        f.write("```bash\npython {ID}-pymol2.py\n```\n")
        f.write(
            "So, for cage A1, inside results/KVFinder-suite/A1/KV_Files directory, run:\n\n"
        )
        f.write("```bash\npython A1-pymol2.py\n```\n")

    # Run KVFinder-suite for all files
    print("> pyKVFinder (v0.4.5)")
    for molecule, s, po, rd, vc in zip(
        molecules, step, probe_out, removal_distance, volume_cutoff
    ):
        print(molecule)

        # Prepare different parameters for pyKVFinder and parKVFinder
        s, _ = _split(s)
        po, _ = _split(po)
        rd, _ = _split(rd)
        vc, _ = _split(vc)

        # Run pyKVFinder
        _run_pyKVFinder(molecule, s, po, rd, vc, basedir)

    print("> parKVFinder (v1.1.4)")
    for molecule, s, po, rd, vc in zip(
        molecules, step, probe_out, removal_distance, volume_cutoff
    ):
        print(molecule)

        # Prepare different parameters for pyKVFinder and parKVFinder
        _, s = _split(s)
        _, po = _split(po)
        _, rd = _split(rd)
        _, vc = _split(vc)

        # Run parKVFinder
        _run_parKVFinder(molecule, s, po, rd, vc, basedir)
