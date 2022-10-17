import os
import zipfile
from typing import List, Union
from biobb_vs.fpocket.fpocket_run import fpocket_run


# fpocket_all_pockets = "fpocket_all_pockets.zip"
# fpocket_summary = "fpocket_summary.json"
# prop = {"min_radius": 3, "max_radius": 6, "num_spheres": 35}

# fpocket_run(
#     input_pdb_path="./data/A1.pdb",
#     output_pockets_zip=fpocket_all_pockets,
#     output_summary=fpocket_summary,
#     properties=prop,
# )


def _run_fpocket(
    pdb: str,
    min_radius: float = 3.2,
    max_radius: float = 6.4,
    num_spheres: int = 15,
    basedir: str = ".",
) -> None:
    # Base Name
    basename = os.path.basename(pdb).replace(".pdb", "")

    # Basedir
    os.makedirs(os.path.join(basedir, basename), exist_ok=True)

    # Output directly
    output = os.path.join(basedir, f"{basename}.zip")

    # Run fpocket
    fpocket_run(
        input_pdb_path=pdb,
        output_pockets_zip=output,
        output_summary=os.path.join(basedir, basename, f"summary.json"),
        properties={
            "min_radius": min_radius,
            "max_radius": max_radius,
            "num_spheres": num_spheres,
            "restart": True
        },
    )

    # Unzip results
    with zipfile.ZipFile(output, 'r') as zip:
        zip.extractall(os.path.join(basedir, basename))


def run(
    pdbs: List[str],
    min_radius: Union[float, List[float]] = 3.2,
    max_radius: Union[float, List[float]] = 6.4,
    num_spheres: Union[int, List[int]] = 15,
) -> None:
    # Check arguments
    if type(pdbs) not in [list]:
        raise TypeError("`pdbs` must be a list of PDB files.")
    if any([not pdb.endswith(".pdb") for pdb in pdbs]):
        raise ValueError("`pdbs` must contain only PDB files.")
    if any([not os.path.isfile(pdb) for pdb in pdbs]):
        raise ValueError("`pdbs` must contain valid PDB files.")
    if type(min_radius) not in [float, list]:
        raise TypeError("`min_radius` must be a float or a list of them.")
    elif type(min_radius) is float:
        min_radius = [min_radius] * len(pdbs)
    elif len(min_radius) != len(pdbs):
        raise Exception("`min_radius` must have the same length as `pdbs`.")
    if type(max_radius) not in [float, list]:
        raise TypeError("`max_radius` must be a float or a list of them.")
    elif type(max_radius) is float:
        max_radius = [max_radius] * len(pdbs)
    elif len(max_radius) != len(pdbs):
        raise Exception("`max_radius` must have the same length as `pdbs`.")
    if type(num_spheres) not in [float, list]:
        raise TypeError("`num_spheres` must be a float or a list of them.")
    elif type(num_spheres) is float:
        num_spheres = [num_spheres] * len(pdbs)
    elif len(num_spheres) != len(pdbs):
        raise Exception("`num_spheres` must have the same length as `pdbs`.")

    # Create basedir
    basedir = "./results/fpocket"
    os.makedirs(basedir, exist_ok=True)

    # Run fpocket for all files
    print("> fpocket (v3.1.4.2)")
    for pdb, m, M, n in zip(pdbs, min_radius, max_radius, num_spheres):
        _run_fpocket(pdb, m, M, n, basedir)

    # Remove zip, log and error files
    for fn in os.listdir(basedir):
        if fn.endswith('.zip'):
            os.remove(os.path.join(basedir, fn))
    for fn in os.listdir():
        if fn.endswith('.out') or fn.endswith('.err'):
            os.remove(fn)
