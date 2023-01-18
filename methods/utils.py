#!/usr/env/bin python3
import math
import os
import pathlib
import re
from typing import Any, List, Dict, Optional, Union

import _pyKVFinder
import numpy
import pandas
from openbabel import pybel
from pyKVFinder import (
    read_pdb,
    read_vdw,
    read_xyz,
    get_vertices,
    detect,
    spatial,
    depth,
    export,
)
from pyKVFinder.grid import _get_dimensions, _get_sincos
from pyKVFinder.utils import VDW


def atoi(self, text: str) -> int:
    return int(text) if text.isdigit() else text


def natural_keys(text: str) -> List[str]:
    return [atoi(c) for c in re.split(r"(\d+)", text)]


def sorting(molecules: List[str]) -> List[str]:
    convert = lambda text: int(text) if text.isdigit() else text
    alphanumeric = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(molecules, key=alphanumeric)


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


class Molecule(object):
    def __init__(
        self,
        molecule: Union[str, pathlib.Path],
        radii: Union[str, pathlib.Path, Dict[str, Any]] = None,
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
        if radii is None:
            # default
            self._radii = read_vdw(VDW)
        elif type(radii) in [str, pathlib.Path]:
            # vdw file
            self._radii = read_vdw(radii)
        elif type(radii) in [dict]:
            # Processed dictionary
            self._radii = radii

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
        export(
            fn, (self.grid == 0).astype(numpy.int32) * 2, None, self.vertices, self.step
        )


class guests(object):
    @staticmethod
    def volume(
        molecules: List[str], step: float = 0.1, radii: Optional[List[Any]] = None
    ) -> pandas.DataFrame:
        # Check arguments
        if type(molecules) not in [list]:
            raise TypeError("`molecules` must be a list of PDB files.")
        if any([not molecule.endswith(".pdb") for molecule in molecules]):
            raise ValueError("`molecules` must contain only PDB files.")
        if any([not os.path.isfile(molecule) for molecule in molecules]):
            raise ValueError("`pdbs` must contain valid PDB files.")
        if type(step) not in [float]:
            raise TypeError("`step` must be a positive float.")
        if radii is not None:
            if type(radii) not in [list]:
                raise TypeError(
                    "`radii` must be a list of vdW radii files or dictionaries."
                )

        # Calculate volumes
        volumes = {}
        for molecule, r in zip(molecules, radii):
            print(f"> Guest volume: {os.path.basename(molecule).rstrip('.pdb')}")

            # Create empty values
            volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"] = numpy.zeros(3)

            # Load molecule
            mol = Molecule(molecule, step=step, radii=r)

            # van der Waals volume
            mol.vdw(padding=3)
            volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"][
                0
            ] = _pyKVFinder._volume(
                (mol.grid == 0).astype(numpy.int32) * 2, step, 1, 16
            )
            mol.save(
                f"results/guests/{os.path.basename(molecule).rstrip('.pdb').split('-')[0]}/{os.path.basename(molecule).rstrip('.pdb')}.vdw.pdb"
            )

            # Surface Excluded Surface (SES) volume
            mol.surface(probe=1.4, surface="SES", padding=3)
            volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"][
                1
            ] = _pyKVFinder._volume(
                (mol.grid == 0).astype(numpy.int32) * 2, step, 1, 16
            )
            mol.save(
                f"results/guests/{os.path.basename(molecule).rstrip('.pdb').split('-')[0]}/{os.path.basename(molecule).rstrip('.pdb')}.ses.pdb"
            )

            # Solvent Accessible Surface (SAS) volume
            mol.surface(probe=1.4, surface="SAS", padding=3)
            volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"][
                2
            ] = _pyKVFinder._volume(
                (mol.grid == 0).astype(numpy.int32) * 2, step, 1, 16
            )
            mol.save(
                f"results/guests/{os.path.basename(molecule).rstrip('.pdb').split('-')[0]}/{os.path.basename(molecule).rstrip('.pdb')}.sas.pdb"
            )

        # Convert to pandas
        volumes = pandas.DataFrame(volumes, index=["vdW", "SES", "SAS"]).T

        return volumes.round(2)

    @staticmethod
    def void(
        molecules: List[str],
        probe_out: float = 4.0,
        step: float = 0.1,
        radii: Optional[List[Any]] = None,
        **kwargs,
    ) -> pandas.DataFrame:
        # Check arguments
        if type(molecules) not in [list]:
            raise TypeError("`molecules` must be a list of PDB files.")
        if any([not molecule.endswith(".pdb") for molecule in molecules]):
            raise ValueError("`molecules` must contain only PDB files.")
        if any([not os.path.isfile(molecule) for molecule in molecules]):
            raise ValueError("`pdbs` must contain valid PDB files.")
        if type(probe_out) not in [float]:
            raise TypeError("`probe_out` must be a float.")
        if type(step) not in [float]:
            raise TypeError("`step` must be a positive float.")
        if radii is not None:
            if type(radii) not in [list]:
                raise TypeError(
                    "`radii` must be a list of vdW radii files or dictionaries."
                )

        # Calculate volumes
        volumes = {}
        for molecule, r in zip(molecules, radii):
            print(f"> Void Volume: {os.path.basename(molecule).rstrip('.pdb')}")

            # Atomic information
            atomic = read_pdb(molecule, r)

            # Grid dimensions
            vertices = get_vertices(atomic, probe_out=probe_out)

            # Cavity detection
            _, cavities = detect(
                atomic,
                vertices,
                step=step,
                probe_out=probe_out,
                **kwargs,
            )

            # Spatial characterization
            surface, volume, _ = spatial(cavities, step=step)

            # Depth characterization
            depths, _, _ = depth(cavities, step=step)

            # Export cavities
            export(
                f"results/guests/{os.path.basename(molecule).rstrip('.pdb').split('-')[0]}/{os.path.basename(molecule).rstrip('.pdb')}.void.pdb",
                cavities,
                surface,
                vertices,
                step,
                B=depths,
            )

            # Save volume
            volumes[f"{os.path.basename(molecule).rstrip('.pdb')}"] = max(
                volume.values()
            )

        # Convert to pandas
        volumes = pandas.DataFrame(volumes, index=["void"]).T

        return volumes.round(2)
