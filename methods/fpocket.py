import os
import zipfile
from typing import List, Union
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_select import fpocket_select


def _run_fpocket(
    pdb: str,
    min_radius: float = 3.4,
    max_radius: float = 6.2,
    num_spheres: int = 15,
    selection: int = 1,
    basedir: str = ".",
) -> None:
    # Base Name
    basename = os.path.basename(pdb).replace(".pdb", "")

    # Basedir
    os.makedirs(os.path.join(basedir, basename), exist_ok=True)

    # Output directly
    output = os.path.join(basedir, f"{basename}.zip")

    # Run fpocket
    if pdb.find("O2.pdb") < 0:
        fpocket_run(
            input_pdb_path=pdb,
            output_pockets_zip=output,
            output_summary=os.path.join(basedir, basename, f"summary.json"),
            properties={
                "min_radius": min_radius,
                "max_radius": max_radius,
                "num_spheres": num_spheres,
            },
        )

        # Select
        fpocket_select(
            input_pockets_zip=output,
            output_pocket_pdb=os.path.join(basedir, basename, f"pocket{selection}_atm.pdb"),
            output_pocket_pqr=os.path.join(basedir, basename, f"pocket{selection}_vert.pqr"),
            properties={"pocket": selection},
        )


def run(
    pdbs: List[str],
    min_radius: Union[float, List[float]] = 3.2,
    max_radius: Union[float, List[float]] = 6.4,
    num_spheres: Union[int, List[int]] = 15,
    selection: Union[int, List[int]] = 1,
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
    if type(num_spheres) not in [int, list]:
        raise TypeError("`num_spheres` must be a integer or a list of them.")
    elif type(num_spheres) is int:
        num_spheres = [num_spheres] * len(pdbs)
    elif len(num_spheres) != len(pdbs):
        raise Exception("`num_spheres` must have the same length as `pdbs`.")
    if type(selection) not in [int, list]:
        raise TypeError("`selection` must be a integer or a list of them.")
    elif type(selection) is int:
        selection = [selection] * len(pdbs)
    elif len(selection) != len(pdbs):
        raise Exception("`selection` must have the same length as `pdbs`.")

    # Create basedir
    basedir = "./results/fpocket"
    os.makedirs(basedir, exist_ok=True)

    # Run fpocket for all files
    print("> fpocket (v3.1.4.2)")
    for pdb, m, M, n, s in zip(pdbs, min_radius, max_radius, num_spheres, selection):
        _run_fpocket(pdb, m, M, n, s, basedir)

    # Remove log and error files
    for fn in os.listdir():
        if fn.endswith(".out") or fn.endswith(".err"):
            os.remove(fn)
