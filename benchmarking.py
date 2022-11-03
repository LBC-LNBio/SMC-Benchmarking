#!/usr/bin/env python3
import os
from typing import List
from openbabel import pybel
from methods import KVsuite, fpocket, pywindow, GHECOM, POVME


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
     # Get XYZ files in guests
    xyzs = [
        os.path.join("./guests", f)
        for f in sorted(os.listdir("./guests"))
        if f.endswith(".xyz")
    ]
 
    # Convert guests XYZ to PDB files
    xyz2pdb(xyzs)
 
    # Get XYZ files in hosts
    xyzs = [
        os.path.join("./hosts", f)
        for f in sorted(os.listdir("./hosts"))
        if f.endswith(".xyz")
    ]

    # Convert hosts XYZ to PDB files
    pdbs = xyz2pdb(xyzs)

    print("[==> Benchmarking methods: ")
    # KVFinder suite
    # parKVFinder documentation: https://github.com/LBC-LNBio/parKVFinder/wiki
    # pyKVFinder documentation: https://lbc-lnbio.github.io/pyKVFinder/
    KVsuite.run(
        pdbs,
        step=[0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.6, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25],
        probe_out=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0],
        removal_distance=[0.75, 2.0, 1.75, 2.0, 2.0, 1.0, 1.25, 3.5, 2.0, 1.5, 1.25, 0.5, 2.0, 1.75],
        volume_cutoff=[80.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 80.0, 20.0, 5.0, 20.0],
    )

    # Fpocket
    # User manual: https://fpocket.sourceforge.net/manual_fpocket2.pdf
    fpocket.run(
        pdbs,
        min_radius=[3.4, 3.4, 3.4, 3.4, 3.4, 2.0, 2.0, 3.4, 3.4, 3.7, 3.4, 3.4, 3.4, 4.0],
        max_radius=[8.0, 8.0, 8.0, 8.0, 8.0, 6.2, 6.2, 40., 6.2, 8.0, 6.2, 6.2, 6.2, 7.0],
        num_spheres=[15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15],
        selection=[1, 1, 1, 1, 1, 1, 1, 37, 1, 1, 1, 1, 1, 7],
    )

    # GHECOM
    # Documentation: https://pdbj.org/ghecom/README_ghecom.html
    GHECOM.run(
        pdbs,
        gws=[0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80],
        rlxs=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.00, 10.0, 10.0, 10.0, 10.00, 10.0, 10.0],
    )

    # pywindow
    # Documentation: https://marcinmiklitz.github.io/pywindow/
    pywindow.run(pdbs)

    # POVME
    # Documentation: https://github.com/POVME/POVME3
    POVME.run(
        pdbs,
        inclusion_regions=[
            [12.171, 12.171, 12.171, 10],
            [0.000, 0.000, 20.583, 6],
            [-4.729, 18.156, 7.611, 6.5],
            [20.538, 0.000, 63.698, 6.5],
            [18.745, 18.747, 18.753, 6.5],
            [6.927, 5.681, 6.603, 10],
            [13.148, 19.722, 3.994, 6.5],
            [37.141, 17.836, 18.198, 25],
            [18.038, 7.504, 11.476, 5],
            [11.196, 6.464, 81.912, 5],
            [4, 6, 12, 16, 15, 15],
            [4.201, 12.171, 12.171, 3],
            [-4.846, 17.740, 7.438, 7],
            [6.819, 9.706, 0.000, 8],
        ],
        grid_spacings=[
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            2.0,
            1.0,
            1.0,
            1.0,
            0.5,
            1.0,
            1.0,
        ],
        contiguous_points_criterias=3,
    )

    # McVol

    # CAVER

    # MoloVol
