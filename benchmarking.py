import os
from typing import List
from openbabel import pybel
from methods import KVsuite, fpocket, pywindow, GHECOM


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
    # parKVFinder documentation: https://github.com/LBC-LNBio/parKVFinder/wiki
    # pyKVFinder documentation: https://lbc-lnbio.github.io/pyKVFinder/
    KVsuite.run(
        pdbs,
        step=0.25,
        probe_out=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0],
        removal_distance=[0.75, 2.0, 1.75, 2.0, 2.0, 1.0, 1.25, 3.5, 2.0, 1.5, 1.25, 0.5, 2.0, 1.5],
        volume_cutoff=[80.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 80.0, 20.0, 5.0, 5.0],
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

    # McVol

    # CAVER

    # MoloVol
