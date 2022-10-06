import os
import subprocess
from typing import List, Union
import toml
import pyKVFinder


def _run_pyKVFinder(
    pdb: str,
    step: float = 0.6,
    probe_out: float = 8.0,
    removal_distance: float = 2.4,
    volume_cutoff: float = 5.0,
    basedir: str = ".",
):
    # Base Name
    basename = os.path.basename(pdb).replace(".pdb", "")

    # Atomic information
    atomic = pyKVFinder.read_pdb(pdb)
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
    opdb = os.path.join(basedir, f"{basename}.KVFinder.output.pdb")
    pyKVFinder.export(opdb, cavities, surface, vertices, B=depths, step=step)
    osurf = os.path.join(basedir, f"{basename}.surface.pdb")
    surface_representation = (cavities == 0).astype(int) * 2
    pyKVFinder.export(
        osurf, surface_representation * 2, None, vertices, B=depths, step=step
    )

    # Write results
    oresults = os.path.join(basedir, f"{basename}.KVFinder.results.toml")
    pyKVFinder.write_results(
        oresults,
        pdb,
        None,
        opdb,
        volume=volume,
        area=area,
        max_depth=max_depth,
        avg_depth=avg_depth,
        step=step,
    )


def _run_parKVFinder(
    pdb: str,
    step: float = 0.6,
    probe_out: float = 8.0,
    removal_distance: float = 2.4,
    volume_cutoff: float = 5.0,
    basedir: str = ".",
):
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
        f.write(f'pdb = "{os.path.realpath(pdb)}"\n')
        f.write("# The path of the output directory.\n")
        f.write(f'output = "{os.path.realpath(basedir)}"\n')
        f.write("# Base name for output files.\n")
        f.write(f"base_name = \"{os.path.basename(pdb).replace('.pdb', '')}\"\n")
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

    # Remove parameters.toml
    os.remove("parameters.toml")


def run(
    pdbs: List[str],
    step: Union[float, List[float]] = 0.6,
    probe_out: Union[float, List[float]] = 8.0,
    removal_distance: Union[float, List[float]] = 2.4,
    volume_cutoff: Union[float, List[float]] = 5.0,
) -> None:
    # Check arguments
    if type(pdbs) not in [list]:
        raise TypeError("`pdbs` must be a list of PDB files.")
    if any([not pdb.endswith(".pdb") for pdb in pdbs]):
        raise ValueError("`pdbs` must contain only PDB files.")
    if any([not os.path.isfile(pdb) for pdb in pdbs]):
        raise ValueError("`pdbs` must contain valid PDB files.")
    if type(probe_out) not in [float, list]:
        raise TypeError("`probe_out` must be a float or a list of them.")
    elif type(probe_out) is float:
        probe_out = [probe_out] * len(pdbs)
    elif len(probe_out) != len(pdbs):
        raise Exception("`probe_out` must have the same length as `pdbs`.")
    if type(step) not in [float, list]:
        raise TypeError("`step` must be a float or a list of them.")
    elif type(step) is float:
        step = [step] * len(pdbs)
    elif len(step) != len(pdbs):
        raise Exception("`step` must have the same length as `pdbs`.")

    # Create basedir
    basedir = "./results/KVFinder-suite"
    os.makedirs(basedir, exist_ok=True)

    # Run KVFinder-suite for all files
    print("> pyKVFinder (v0.4.4)")
    for pdb, s, po, rd, vc in zip(
        pdbs, step, probe_out, removal_distance, volume_cutoff
    ):
        _run_pyKVFinder(pdb, s, po, rd, vc, basedir)

    print("> parKVFinder (v1.1.4)")
    for pdb, s, po, rd, vc in zip(
        pdbs, step, probe_out, removal_distance, volume_cutoff
    ):

        _run_parKVFinder(pdb, s, po, rd, vc, basedir)
