# <ins>S</ins>upra<ins>M</ins>olecular <ins>C</ins>ages <ins>Benchmarking</ins>

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/LBC-LNBio/moc-benchmarking) 
![GitHub](https://img.shields.io/github/license/LBC-LNBio/moc-benchmarking)

A set of well-known supramolecular cages were selected from the literature to evaluate and validate benchmarking methods in the supramolecular cage context.

## Benchmarking methods

- KVFinder suite ([parKVFinder](https://doi.org/10.1016/j.softx.2020.100606) and [pyKVFinder](https://doi.org/10.1186/s12859-021-04519-4));
- [Fpocket](https://doi.org/10.1186/1471-2105-10-168);
- [McVol](https://doi.org/10.1007/s00894-009-0541-y);
- [pywindow](https://doi.org/10.1021/acs.jcim.8b00490 );
- [POVME](https://doi.org/10.1021/acs.jctc.7b00500);
- [GHECOM](https://doi.org/10.1002/prot.22639);
- [CAVER](https://doi.org/10.1093/bioinformatics/bty386);
- [MoloVol](https://doi.org/10.1107/S1600576722004988);
- [PyVol](https://doi.org/10.1101/816702).

## Requirements

- [openbabel](https://pypi.org/project/openbabel/) (v3.1.1.1): a chemical toolbox designed to speak the many languages of chemical data. A wrapper that is automatically using the SWIG package and provide access to almost all of the Open Babel interfaces via Python.

- [pybel](https://pypi.org/project/pybel/) (v0.15.5): a pure Python package for parsing and handling biological encoded in Biological Expression Language. A lightweight wrapper around openbabel module, that provides a more convenient and Pythonic ways to access the Open Babel toolkit.

To install Open Babel and PyBEL, run:

```bash
sudo apt install python3-openbabel
```

To use Open Babel and PyBEL in Python ecosystem, run:

```python
from openbabel import openbabel
from openbabel import pybel
```

### Python packages

- [toml](https://pypi.org/project/toml) (v0.10.2): a Python library for parsing and creating TOML.

- [pyKVFinder](https://pypi.org/project/toml) (v0.4.4): a Python package for detecting and characterizing biomolecular cavities.

- [biobb_vs](https://pypi.org/project/biobb_vs) (v3.8.1): a collection to perform virtual screening studies. Biobb (BioExcel building blocks) packages are Python building blocks that create new layer of compatibility and interoperability over popular bioinformatics tools.

To install toml, pyKVFinder, biobb_vs, run:

```bash
pip install -r requirements.txt
```

## License

The software is licensed under the terms of the MIT License and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
