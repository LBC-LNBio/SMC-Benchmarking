import os
from typing import List
from openbabel import pybel
from methods import KVsuite, fpocket


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
        # Convert residue names from UNL to UNK
        for atom in mol.atoms:
            atom.residue.OBResidue.SetName("UNK")
        mol.write("pdb", pdb, overwrite=True)

    return pdbs


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
    # KVFinder suite
    KVsuite.run(
        pdbs,
        step=0.25,
        probe_out=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0],
        removal_distance=[0.75, 2.0, 1.75, 2.0, 2.0, 1.0, 1.25, 3.5, 2.0, 1.5, 1.25, 0.5, 2.0, 1.5],
        volume_cutoff=[80.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 80.0, 20.0, 5.0, 5.0],
    )

    # Fpocket
    pdbs = ["./data/A1.pdb"]
    fpocket.run(pdbs, min_radius=[3.2], max_radius=[6.4], num_spheres=[15])

    # McVol

    # pywindow

    # POVME

    # GHECOM

    # CAVER

    # MoloVol

    # PyVol
