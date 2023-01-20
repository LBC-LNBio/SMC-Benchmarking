#!/usr/env/bin python3
import os
import re
from typing import List


from methods import GHECOM, POVME, KVsuite, MoloVol, fpocket, pywindow, utils


if __name__ == "__main__":
    # print("[==> Estimating guests' vdW volume")
    # # Get guests' files
    # guests = [os.path.join("./guests", f) for f in os.listdir("./guests")]
    # guests = utils.sorting(guests)

    # # vdW radii dictionary (ChimeraX)
    # # REFERENCE: UCSF ChimeraX: Structure visualization for researchers, educators, and developers. Pettersen EF, Goddard TD, Huang CC, Meng EC, Couch GS, Croll TI, Morris JH, Ferrin TE. Protein Sci. 2021 Jan;30(1):70-82.
    # # NOTE: https://www.cgl.ucsf.edu/chimerax/docs/user/radii.html
    # radii = [
    #     {"GEN": {"C": 1.7, "H": 1, "N": 1.625}},  # Et4N
    #     {"GEN": {"C": 1.7, "H": 1, "N": 1.625}},  # BnNMe3
    #     {"GEN": {"C": 1.7, "H": 1, "CO": 0.56}},  # CoCp2
    #     {"GEN": {"C": 1.7, "H": 1, "CO": 0.56}},  # CoCp*2
    #     {"GEN": {"B": 1.66, "F": 1.28}},  # BF4
    #     {"GEN": {"CL": 1.98, "O": 1.46}},  # ClO4
    #     {"GEN": {"C": 1.61}},  # C60
    #     {"GEN": {"C": 1.61}},  # C60
    #     {"GEN": {"C": 1.7, "H": 1, "O": 1.48}},  # Ad
    #     {"GEN": {"C": 1.7, "H": 1, "N": 1.625}},  # Et4N
    #     {"GEN": {"C": 1.7, "H": 1, "N": 1.625, 'O': 1.46}},  # AQ
    #     {"GEN": {"C": 1.7, "H": 1, "N": 1.625}},  # 2 guests
    #     {"GEN": {"C": 1.61}},  # C60
    #     {"GEN": {"C": 1.61}},  # C70
    # ]

    # # vdW radii dictionary (PyMOL v2.5.0)
    # # REFERENCE: The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.
    # # NOTE: https://pymolwiki.org/index.php/Category:Uncategorized
    # # radii = [{"C": 1.7, "H": 1.2, "N": 1.55, "O": 1.52, "B":1.85, "F":1.47, "CL": 1.75, "CO":1.8}] * 14

    # # Estimate guests' volumes
    # volumes = utils.guests.volume(guests, step=0.1, radii=radii)

    # # Add void volume to ball-shaped cages, ie C60 (B7, B8, 13) and C70 (B14)
    # voids = [
    #     "./guests/B7-C60.pdb",
    #     "./guests/B8-C60.pdb",
    #     "./guests/B13-C60.pdb",
    #     "./guests/B14-C70.pdb",
    # ]
    # radii = [radii[index] for index in [6, 7, 12, 13]]
    # additions = utils.guests.void(voids, removal_distance=0.0, volume_cutoff=10.0, radii=radii)

    # for index in additions.index:
    #     volumes.loc[index, :] = volumes.loc[index, :] + additions.loc[index, "void"]

    # # Save guests' volume
    # volumes.to_csv("results/guest-volume.csv")

    print("[==> Converting hosts from XYZ to PDB")
    # Get hosts' XYZ files
    xyzs = [
        os.path.join("./hosts", f)
        for f in sorted(os.listdir("./hosts"))
        if f.endswith(".xyz")
    ]

    # Convert hosts XYZ to PDB files
    pdbs = utils.xyz2pdb(xyzs)
    pdbs = utils.sorting(pdbs)

    print("[==> Benchmarking methods: ")
    # KVFinder suite
    # parKVFinder documentation: https://github.com/LBC-LNBio/parKVFinder/wiki
    # pyKVFinder documentation: https://lbc-lnbio.github.io/pyKVFinder/
    KVsuite.run(
        pdbs,
        step=[0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.6, 0.25, 0.25, 0.25, 0.25, 0.25],
        probe_out=[10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 6, 6, 10, 10, 20, 10, 10, 10, 20, 10],
        removal_distance=[0.75, 2, 1.75, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1.5, 1.5, 1.5, 1, 1.25, 3.5, 2, 1.5, 1.25, [0.5, 1.25], 1.75],
        volume_cutoff=[80, 5, 5, 20, 5, 5, 5, 110, 25, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 80, 5, 20],
    )

    # Fpocket
    # User manual: https://fpocket.sourceforge.net/manual_fpocket2.pdf
    fpocket.run(
        pdbs,
        min_radius=[3.4, 3.4, 3.4, 3.4, 3.4, 3, 3.4, 3.4, 3.4, 3.4, 3.4, 3.7, 3.4, 3.4, 3.4, 3.4, 2.0, 3.4, 3.4, 3.7, 3.4, 3.4, 4.0],
        max_radius=[8, 8, 8, 8, 8, 6.2, 8, 8, 8, 8, 8, 6.2, 8, 8, 8, 8, 6.2, 40, 6.2, 8, 6.2, 6.2, 7],
        num_spheres=[15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15],
        selection=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 37, 1, 1, 1, 1, 7],
    )

    # GHECOM
    # Documentation: https://pdbj.org/ghecom/README_ghecom.html
    GHECOM.run(
        pdbs,
        gws=[0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8],
        rlxs=[10, 10, 10, 10, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 10, 10, 10, 10, 10],
    )

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
            [13.179, 9.119, 7.562, 6.5],
            [1.589, 4.259, 4.046, 6.5],
            [22.614, 23.073, 15.430, 8],
            [-9.087,12.579, -10.771, 6.5],
            [18.603, 4.743, 21.449, 6.5],
            [-20.465, -20.518, -20.505, 6.5],
            [2.482, 8.476, 20.911, 6.5],
            [-1.705, 17.095, 17.441, 6.5],
            [7.977, 28.209, 9.022, 6.5],
            [15.111, 16.105, 14.164, 6.5],
            [-6.470, 29.785, 13.872, 10],
            [13.148, 19.722, 3.994, 6.5],
            [37.141, 17.836, 18.198, 25],
            [18.038, 7.504, 11.476, 5],
            [11.196, 6.464, 81.912, 5],
            [4, 6, 12, 16, 15, 15],
            [4.201, 12.171, 12.171, 3],
            [6.819, 9.706, 0.000, 8],
        ],
        grid_spacings=[1, 1, 1, 1, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 0.5, 1],
        contiguous_points_criterias=3,
    )

    # pywindow
    # Documentation: https://marcinmiklitz.github.io/pywindow/
    # WARNING: numpy<=1.23.5
    # NOTE: numpy>=1.24.0 returns:
    # "ValueError: setting an array element with a sequence. The requested array has an inhomogeneous shape after 1 dimensions. The detected shape was (3,) + inhomogeneous part"
    pywindow.run(pdbs)

    # MoloVol
    # Documentation: https://molovol.com/
    MoloVol.run(
        pdbs,
        grid_spacings=0.6,
        small_probes=1.4,
        large_probes=[8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 5, 20, 5, 5, 8, 15, 5],
    )

    # CAVER 3.0
    # Documentation: https://caver.cz
    print("> CAVER (v3.0.3)")
    print(
        "Creating output directories. Now, run CAVER v3.0.3 in pymol to compute tunnels."
    )
    os.makedirs("./results/CAVER", exist_ok=True)
    for d in pdbs:
        d = os.path.basename(d).strip(".pdb")
        os.makedirs(os.path.join("./results/CAVER", d), exist_ok=True)

    print(
        "Here, we defined the starting point as the center of mass of the target supramolecular cage\n"
    )
    print(
        "{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format(
            "Supramolecular cage",
            "Minimum Probe Radius",
            "Shell Depth",
            "Shell Radius",
            "Clustering threshold",
        )
    )
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("A1", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B1", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B2", 0.7, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B3", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B4", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B5", 0.6, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B6", 0.6, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B7", 0.7, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B8", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B9", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B10", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B11", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B12", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B13", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("B14", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("C1", 0.7, 4.0, 10.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("F1", 0.9, 4.0, 5.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("F2", 0.9, 4.0, 3.0, 5.0))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("H1", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("N1", 0.9, 4.0, 7.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("O1", 0.9, 4.0, 3.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("O2", 0.9, 4.0, 2.0, 3.5))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".format("W1", 0.9, 4.0, 3.0, 3.5))

    # CAVER Analyst 2.0
    # Documentation: https://caver.cz
    print("> CAVER Analyst 2.0 (CAVER v3.0.2)")
    print("Now, run CAVER Analyst 2.0 to compute cavity volumes, by typing:")
    print("$ bash etc/CAVER/caver_analyst2/bin/caver_analyst\n")
    print(
        "{:<20}\t{:<20}\t{:<20}".format(
            "Supramolecular cage",
            "Probe",
            "Large Probe",
        )
    )
    print("{:<20}\t{:<20}\t{:<20}".format("A1", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B1", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B2", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B3", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B4", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B5", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B6", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B7", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B8", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B9", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B10", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B11", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B12", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B13", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("B14", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("C1", 1.4, 4.0))
    print("{:<20}\t{:<20}\t{:<20}".format("F1", 1.4, 4.0))
    print("{:<20}\t{:<20}\t{:<20}".format("F2", 1.4, 5.3))
    print("{:<20}\t{:<20}\t{:<20}".format("H1", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("N1", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("O1", 1.4, 3.0))
    print("{:<20}\t{:<20}\t{:<20}".format("O2", "NA", "NA"))
    print("{:<20}\t{:<20}\t{:<20}".format("W1", 1.4, 3.0))
