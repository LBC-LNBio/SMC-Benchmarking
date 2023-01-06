#!/usr/bin/env python3
import math
import os
import pathlib
import re
import sys
from typing import Any, List, Optional, Union

import _pyKVFinder
import numpy
import pandas
import toml
from pyKVFinder import export, get_vertices, read_pdb, read_vdw, read_xyz
from pyKVFinder.grid import _get_dimensions, _get_sincos
from pyKVFinder.utils import VDW

sys.path.append("../")
from benchmarking import xyz2pdb
from methods.KVsuite import _run_pyKVFinder, _split


def atoi(text: str) -> int:
    return int(text) if text.isdigit() else text


def natural_keys(text: str) -> List[str]:
    return [atoi(c) for c in re.split(r"(\d+)", text)]


class Molecule(object):
    def __init__(
        self,
        molecule: Union[str, pathlib.Path],
        radii: Union[str, pathlib.Path] = None,
        step: Union[str, pathlib.Path] = 0.6,
        model: Optional[int] = None,
        nthreads: Optional[int] = None,
        verbose: bool = False,
    ):
        # Attributes
        self._grid = None
        self._padding = None
        self._probe = None
        self._representation = None
        self._vertices = None
        self._dim = None
        self._rotation = None
        self.verbose = verbose

        # Molecule
        if type(molecule) not in [str, pathlib.Path]:
            raise TypeError("`fn` must be a string or a pathlib.Path.")
        self._molecule = os.path.realpath(molecule)

        # van der Waals radii
        if self.verbose:
            print("> Loading van der Waals radii")
        if radii is not None:
            self._radii = read_vdw(radii)
        else:
            self._radii = read_vdw(VDW)

        # Step
        if type(step) not in [int, float]:
            raise TypeError("`step` must be a postive real number.")
        elif step <= 0.0:
            raise ValueError("`step` must be a positive real number.")
        else:
            self._step = step

        # Atomic information
        if self.verbose:
            print("> Reading molecule coordinates")
        if molecule.endswith(".pdb"):
            self._atomic = read_pdb(molecule, self.radii, model)
        elif molecule.endswith(".xyz"):
            self._atomic = read_xyz(molecule, self.radii)
        else:
            raise TypeError("`molecule` must have .pdb or .xyz extension.")

        # Number of threads
        if nthreads is not None:
            if type(nthreads) not in [int]:
                raise TypeError("`nthreads` must be a positive integer.")
            elif nthreads <= 0:
                raise ValueError("`nthreads` must be a positive integer.")
            else:
                self.nthreads = nthreads
        else:
            self.nthreads = os.cpu_count() - 1

    @property
    def molecule(self):
        return self._molecule

    @property
    def radii(self):
        return self._radii

    @property
    def step(self):
        return self._step

    @property
    def atomic(self):
        return self._atomic

    @property
    def xyzr(self):
        return self._atomic[:, 4:].astype(numpy.float64)

    @property
    def vertices(self):
        return self._vertices

    @property
    def p1(self):
        if self._vertices is not None:
            return self._vertices[0]

    @property
    def p2(self):
        if self._vertices is not None:
            return self._vertices[1]

    @property
    def p3(self):
        if self._vertices is not None:
            return self._vertices[2]

    @property
    def p4(self):
        if self._vertices is not None:
            return self._vertices[3]

    @property
    def nx(self):
        if self._dim is not None:
            return self._dim[0]

    @property
    def ny(self):
        if self._dim is not None:
            return self._dim[1]

    @property
    def nz(self):
        if self._dim is not None:
            return self._dim[2]

    @property
    def dim(self):
        return self._dim

    @property
    def rotation(self):
        return self._rotation

    @property
    def padding(self):
        return self._padding

    @property
    def probe(self):
        return self._probe

    @property
    def representation(self):
        return self._representation

    @property
    def grid(self):
        return self._grid

    def _set_grid(self, padding: Optional[float]):
        # Padding
        if padding is not None:
            if type(padding) not in [int, float, numpy.float64]:
                raise TypeError("`step` must be a non-negative real number.")
            elif padding < 0.0:
                raise ValueError("`step` must be a non-negative real number.")
            else:
                self._padding = padding
        else:
            self._padding = self._get_padding()

        # 3D grid
        if self.verbose:
            print("> Calculating 3D grid")
        self._vertices = get_vertices(self.atomic, self.padding, self.step)
        self._dim = _get_dimensions(self.vertices, self.step)
        self._rotation = _get_sincos(self.vertices)
        if self.verbose:
            print(f"p1: {self.vertices[0]}")
            print(f"p2: {self.vertices[1]}")
            print(f"p3: {self.vertices[2]}")
            print(f"p4: {self.vertices[3]}")
            print("nx: {}, ny: {}, nz: {}".format(*self.dim))
            print("sina: {}, sinb: {}, cosa: {}, cosb: {}".format(*self.rotation))

    def _get_padding(self):
        return (self._atomic[:, 4:7].astype(float).ptp(axis=0).max() / 10).round(
            decimals=1
        )

    def vdw(self, padding: Optional[float] = None):
        from _pyKVFinder import _fill_receptor as fill

        # Attributes
        self._representation = "vdW"
        self._probe = None

        # Define 3D grid
        self._set_grid(padding)

        # van der Waals atoms (hard sphere model) to grid
        if self.verbose:
            print("> Inserting atoms with van der Waals radii into 3D grid")
        self._grid = _pyKVFinder._fill_receptor(
            math.prod(self.dim),
            self.nx,
            self.ny,
            self.nz,
            self.xyzr,
            self.p1,
            self.rotation,
            self.step,
            0.0,
            False,
            self.nthreads,
            self.verbose,
        ).reshape(self.nx, self.ny, self.nz)

    def surface(
        self, probe: float = 1.4, surface: str = "SES", padding: Optional[float] = None
    ):
        # Probe
        if type(probe) not in [int, float, numpy.float64]:
            raise TypeError("`probe_out` must be a non-negative real number.")
        elif probe < 0.0:
            raise ValueError("`probe_out` must be a non-negative real number.")
        self._probe = probe

        # Surface
        if surface == "SES":
            if self.verbose:
                print("> Surface representation: Solvent Excluded Surface (SES).")
            self._representation = surface
            surface = True
        elif surface == "SAS":
            if self.verbose:
                print("> Surface representation: Solvent Accessible Surface (SAS).")
            self._representation = surface
            surface = False
        else:
            raise ValueError(f"`surface` must be SAS or SES, not {surface}.")

        # Define 3D grid
        self._set_grid(padding)

        # Molecular surface (SES or SAS) to grid
        self._grid = _pyKVFinder._fill_receptor(
            math.prod(self.dim),
            self.nx,
            self.ny,
            self.nz,
            self.xyzr,
            self.p1,
            self.rotation,
            self.step,
            self.probe,
            surface,
            self.nthreads,
            self.verbose,
        ).reshape(self.nx, self.ny, self.nz)

    def preview(self, **kwargs):
        if self.grid is not None:
            from plotly.express import scatter_3d

            x, y, z = numpy.nonzero(self.grid == 0)
            fig = scatter_3d(x=x, y=y, z=z, **kwargs)
            fig.update_layout(
                scene_xaxis_showticklabels=False,
                scene_yaxis_showticklabels=False,
                scene_zaxis_showticklabels=False,
            )
            fig.show()

    def save(
        self,
        fn: Union[str, pathlib.Path] = "molecule.pdb",
    ):
        # Filename (fn)
        if type(fn) not in [str, pathlib.Path]:
            raise TypeError("`fn` must be a string or a pathlib.Path.")
        os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

        # Save grid to PDB file
        export(fn, (self.grid == 0).astype(numpy.int32), None, self.vertices, self.step)


def estimate_cavity_volume(
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
    if type(probe_out) not in [float, list]:
        raise TypeError("`probe_out` must be a float or a list of them.")
    elif type(probe_out) is float:
        probe_out = [probe_out] * len(molecules)
    elif len(probe_out) != len(molecules):
        raise Exception("`probe_out` must have the same length as `molecules`.")
    if type(step) not in [float, list]:
        raise TypeError("`step` must be a float or a list of them.")
    elif type(step) is float:
        step = [step] * len(molecules)
    elif len(step) != len(molecules):
        raise Exception("`step` must have the same length as `molecules`.")
    if type(removal_distance) not in [float, list]:
        raise TypeError("`removal_distance` must be a float or a list of them.")
    elif type(removal_distance) is float:
        removal_distance = [removal_distance] * len(molecules)
    elif len(removal_distance) != len(molecules):
        raise Exception("`removal_distance` must have the same length as `molecules`.")
    if type(volume_cutoff) not in [float, list]:
        raise TypeError("`volume_cutoff` must be a float or a list of them.")
    elif type(volume_cutoff) is float:
        volume_cutoff = [volume_cutoff] * len(molecules)
    elif len(volume_cutoff) != len(molecules):
        raise Exception("`volume_cutoff` must have the same length as `molecules`.")

    # Run pyKVFinder for all files
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


def estimate_guests_volume(molecules: List[str], step: float = 0.1) -> pandas.DataFrame:
    # Calculate volumes
    volumes = {}
    for molecule in molecules:
        print(f"> {os.path.basename(molecule).rstrip('.pdb')}")

        # Create empty values
        volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"] = numpy.zeros(3)

        # Load molecule
        mol = Molecule(molecule, step=step)

        # van der Waals volume
        mol.vdw(padding=3)
        volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"][
            0
        ] = _pyKVFinder._volume((mol.grid == 0).astype(numpy.int32) * 2, step, 1, 16)
        mol.save(
            f"results/{os.path.basename(molecule).rstrip('.pdb').split('-')[0]}/{os.path.basename(molecule).rstrip('.pdb')}.vdw.pdb"
        )

        # Surface Excluded Surface (SES) volume
        mol.surface(probe=1.4, surface="SES", padding=3)
        volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"][
            1
        ] = _pyKVFinder._volume((mol.grid == 0).astype(numpy.int32) * 2, step, 1, 16)
        mol.save(
            f"results/{os.path.basename(molecule).rstrip('.pdb').split('-')[0]}/{os.path.basename(molecule).rstrip('.pdb')}.ses.pdb"
        )

        # Solvent Accessible Surface (SAS) volume
        mol.surface(probe=1.4, surface="SAS", padding=3)
        volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"][
            2
        ] = _pyKVFinder._volume((mol.grid == 0).astype(numpy.int32) * 2, step, 1, 16)
        mol.save(
            f"results/{os.path.basename(molecule).rstrip('.pdb').split('-')[0]}/{os.path.basename(molecule).rstrip('.pdb')}.sas.pdb"
        )

    # Convert to pandas
    volumes = pandas.DataFrame(volumes, index=["vdW", "SES", "SAS"]).T

    return volumes


def _process_cavity_volume(molecules: List[str], basedir: str) -> pandas.DataFrame:

    volumes = {}
    for molecule in molecules:
        basename = os.path.basename(molecule).rstrip(".pdb")
        results = os.path.join(basedir, basename, f"{basename}.KVFinder.results.toml")

        with open(results, "r") as f:
            volume = toml.load(f)
            volumes[basename] = volume["RESULTS"]["VOLUME"]["KAA"]

    # Convert to pandas
    volumes = pandas.DataFrame(volumes, index=["Cavity"]).T

    return volumes


if __name__ == "__main__":
    # Create basedir
    basedir = "./results"
    os.makedirs(basedir, exist_ok=True)

    print("[==> Estimate hosts' cavity volume")
    cages = [
        os.path.join("./hosts", f)
        for f in sorted(os.listdir("./hosts"))
        if f.endswith(".pdb")
    ]
    cages = sorted(cages, key=lambda x: int(x.split("B")[-1].replace(".pdb", "")))

    # Estimate hosts' cavity volumes
    estimate_cavity_volume(
        cages,
        step=0.25,
        probe_out=[10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 6, 6],
        removal_distance=[2, 1.75, 2, 2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
        volume_cutoff=[5, 5, 5, 5, 5, 5, 110, 25, 5, 5, 5, 5, 5, 5],
    )
    volumes = _process_cavity_volume(cages, basedir)
    volumes.to_csv("results/cavity-volume.csv")

    print("[==> Estimate guests' vdW volume")
    # Get XYZ files in guests
    guests = [
        os.path.join("./guests", f)
        for f in sorted(os.listdir("./guests"))
        if f.endswith(".pdb")
    ]
    guests = sorted(guests, key=lambda x: int(x.split("-")[0].split("B")[-1]))

    # Estimate guests' volumes
    volumes = estimate_guests_volume(guests).round(2)
    volumes.to_csv("results/guest-volume.csv")
