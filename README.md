# <ins>S</ins>upra<ins>M</ins>olecular <ins>C</ins>ages <ins>Benchmarking</ins>

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/LBC-LNBio/SMC-Benchmarking) 
![GitHub](https://img.shields.io/github/license/LBC-LNBio/SMC-Benchmarking)

A set of well-known supramolecular cages were selected from the literature to evaluate and validate benchmarking methods in the supramolecular cage context.

## Benchmarking methods

- KVFinder suite ([parKVFinder](https://doi.org/10.1016/j.softx.2020.100606) and [pyKVFinder](https://doi.org/10.1186/s12859-021-04519-4));
- [Fpocket](https://doi.org/10.1186/1471-2105-10-168);
- [pywindow](https://doi.org/10.1021/acs.jcim.8b00490);
- [POVME](https://doi.org/10.1021/acs.jctc.7b00500);
- [GHECOM](https://doi.org/10.1002/prot.22639);
- [CAVER](https://doi.org/10.1093/bioinformatics/bty386);
- [MoloVol](https://doi.org/10.1107/S1600576722004988).

## System requirements

- [Python 3](https://www.python.org) (v3.10.6): Python is a high-level, interpreted, interactive and object-oriented scripting language. 

- [Python 2](https://www.python.org) (v2.7.18): Python is a high-level, interpreted, interactive and object-oriented scripting language. 

- [openbabel](https://pypi.org/project/openbabel/) (v3.1.1.1): a chemical toolbox designed to speak the many languages of chemical data. A wrapper that is automatically using the SWIG package and provide access to almost all of the Open Babel interfaces via Python.

- [pybel](https://pypi.org/project/pybel/) (v0.15.5): a pure Python package for parsing and handling biological encoded in Biological Expression Language. A lightweight wrapper around openbabel module, that provides a more convenient and Pythonic ways to access the Open Babel toolkit.

- [GHECOM](https://pdbj.org/ghecom/) (21/07/2020): a software for finding multi-scale pockets on protein surfaces using mathematical morphology.

- [MoloVol](https://github.com/molovol/MoloVol) (v1.0.0): a software for volume and surface computations of single molecules and crystallographic unit cells.

- [CAVER](https://www.caver.cz) (v3.0.2): is a software tool for analysis and visualization of tunnels and channels in protein structures.

To install Python 3, Python 2, Pip 3, Pip 2, Open Babel and PyBEL, run:

```bash
sudo apt install python3 python3-pip python3-openbabel python2 python-pip
```

To use Open Babel and PyBEL in Python 3 ecosystem, run:

```python3
from openbabel import openbabel
from openbabel import pybel
```

To install GHECOM, run:

```bash
mkdir -p etc/
mkdir -p etc/ghecom
wget -O etc/ghecom.tar.gz https://pdbj.org/ghecom/cgi-bin/dwnld_src_file.cgi?filename=ghecom-src-20211201.tar.gz
tar -xf etc/ghecom.tar.gz -C etc/ghecom
cd etc/ghecom/src
make
cd ../../../
sudo ln -s `pwd`/etc/ghecom/ghecom /usr/local/bin/ghecom
```

To install MoloVol, run:

```bash
mkdir -p etc/MoloVol
wget -O etc/MoloVol/MoloVol.deb https://github.com/molovol/MoloVol/releases/download/v1.0.0/MoloVol_debian_x84-64_v1.0.deb
wget -O etc/MoloVol/Guide.pdf https://github.com/molovol/MoloVol/releases/download/v1.0.0/MoloVol_Quick_Start_Guide_v1.0.pdf
wget -O etc/MoloVol/Manual.pdf https://github.com/molovol/MoloVol/releases/download/v1.0.0/MoloVol_User_Manual_v1.0.pdf
dpkg -i etc/MoloVol/MoloVol.deb
```

To install CAVER, follow:

1. Run:

```bash
mkdir -p etc/CAVER
wget -O etc/CAVER/caver-pymol-plugin-3.0.3.zip https://github.com/loschmidt/caver-pymol-plugin/archive/v3.0.3.zip
unzip etc/CAVER/caver-pymol-plugin-3.0.3.zip -d etc/CAVER/caver-pymol-plugin-3.0.3
```

2. Start PyMOL
3. Open Plugin manager (Plugin -> Plugin manager)
4. Install plugin (Install new plugin -> Choose file... and select ./etc/CAVER/caver-pymol-plugin-3.0.3/caver-pymol-plugin-3.0.3/Caver3/__init__.py)
5. Restart PyMOL

To install CAVER Analyst 2.0, run:

```bash
mkdir -p etc/CAVER
wget -O etc/CAVER/CAVER-Analyst-2.zip https://www.caver.cz/fil/download/caver_analyst/caver_analyst2_beta2.zip
unzip etc/CAVER/CAVER-Analyst-2.zip -d etc/CAVER
```


## Python 3 requirements

- [toml](https://pypi.org/project/toml) (v0.10.2): a Python library for parsing and creating TOML.

- [pyKVFinder](https://pypi.org/project/pyKVFinder) (v0.4.4): a Python package for detecting and characterizing biomolecular cavities.

- [biobb_vs](https://pypi.org/project/biobb_vs) (v3.8.1): a collection to perform virtual screening studies. Biobb (BioExcel building blocks) packages are Python building blocks that create new layer of compatibility and interoperability over popular bioinformatics tools. Fpocket software (v3.1.4.2) is available inside this Python package.

- [pywindow](https://github.com/marcinmiklitz/pywindow) (v0.0.4): a Python package for the analysis of structural properties of molecular pores, porous organic cages, MOFs and metalorganic cages.

To install toml, pyKVFinder, biobb_vs, run:

```bash
pip3 install -r requirements3.txt
```

To install pywindow, run:

```bash
mkdir -p etc
cd etc
git clone https://github.com/marcinmiklitz/pywindow
cd pywindow/
pip install .
cd ../../
```
## Python 2 requirements

- [povme](https://pypi.org/project/povme/) (v3.0.35): is a Python package characterizing the shape and flexibility of a cavity in detail.

To install povme, run:

```bash
pip2 install -r requirements2.txt
```

## Benchmarking

To execute the benchmarking, run:

```bash
python3 benchmarking.py
```

## License

The software is licensed under the terms of the MIT License and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
