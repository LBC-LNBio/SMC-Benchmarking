import os
from typing import List, Union
from openbabel import pybel
import pyKVFinder

def xyz2pdb(xyzs: List[str]) -> List[str]:
    """Convert XYZ to PDB files.

    Parameters
    ----------
    xyzs : List[Union[str, pathlib.Path]]
        List of XYZ files.

    Returns
    -------
    pdbs : List[str]
        List of PDB files.
    """
    # Check arguments
    if type(xyzs) not in [list]:
        raise TypeError("`xyzs` must be a list of XYZ files.")
    if any([not xyz.endswith(".xyz") for xyz in xyzs]):
        raise ValueError("`xyzs` must contain only XYZ files.")
    if any([not os.path.isfile(xyz) for xyz in xyzs]):
        raise ValueError("`xyzs` must contain valid XYZ files.")

    # Prepare filenames
    pdbs = [xyz.replace("xyz", "pdb") for xyz in xyzs]

    # Convert XYZ to PDB
    for xyz, pdb in zip(xyzs, pdbs):
        mol = next(pybel.readfile("xyz", xyz))
        mol.write("pdb", pdb, overwrite=True)

    return pdbs


def run_KVFinder_suite(
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
    for pdb, s, po, rd, vc in zip(pdbs, step, probe_out, removal_distance, volume_cutoff):
        # Base Name
        basename = os.path.basename(pdb).replace(".pdb", "")

        # Atomic information
        atomic = pyKVFinder.read_pdb(pdb)
        # Grid dimensions
        vertices = pyKVFinder.get_vertices(atomic, probe_out=po)

        # Cavity detection
        _, cavities = pyKVFinder.detect(
            atomic, vertices, step=s, probe_out=po, removal_distance=rd, volume_cutoff=vc
        )

        # Spatial characterization
        surface, volume, area = pyKVFinder.spatial(cavities, step=s)
        # Depth characterization
        depths, max_depth, avg_depth = pyKVFinder.depth(cavities, step=s)

        # Export cavities
        opdb = os.path.join(basedir
            , f"{basename}.KVFinder.output.pdb"
        )
        pyKVFinder.export(opdb, cavities, surface, vertices, B=depths, step=s)
        osurf = os.path.join(basedir, f"{basename}.surface.pdb")
        surface_representation = (cavities == 0).astype(int) * 2
        pyKVFinder.export(osurf, surface_representation * 2, None, vertices, B=depths, step=s)

        # Write results
        oresults = os.path.join(
            basedir, f"{basename}.KVFinder.results.toml"
        )
        pyKVFinder.write_results(
            oresults,
            pdb,
            None,
            opdb,
            volume=volume,
            area=area,
            max_depth=max_depth,
            avg_depth=avg_depth,
            step=s,
        )




if __name__ == "__main__":
    print("[==> Converting XYZ to PDB")

    # Get XYZ files in data
    xyzs = [
        os.path.join("./data", f)
        for f in sorted(os.listdir("./data"))
        if f.endswith(".xyz")
    ]

    # Convert XYZ to PDB files
    pdbs = xyz2pdb(xyzs)

    print("[==> Benchmarking methods: ")
    print("> pyKVFinder (v0.4.4)")
    run_KVFinder_suite(pdbs,
    step=0.25,
    probe_out=[8.0, 20.0, 20.0, 40.0, 20.0, 20.0, 20.0, 20.0],
    removal_distance=[0.0, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4],
    volume_cutoff=[800.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    )
