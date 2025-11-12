# Vasp File Inspector (vfi)

**Vasp File Inspector (vfi)** is a lightweight Python CLI utility for analyzing **VASP structure files** (`POSCAR` and `CONTCAR`)  
It can print and save detailed information about **atomic positions**, **bonding networks**, **unit cell geometry**, **unit cell symmetry**, and **neighbor relationships**.

---

## Installation

Clone and install directly from GitHub:

```bash
git clone https://github.com/joegonzalezMSL/VaspFileInspector.git
cd VaspFileInspector
pip install .
```

Requires Python ≥ 3.9 and the following dependencies, all of which are handled by the pyproject.toml:

`numpy >= 1.24`

`spglib >= 2.0`

## Usage
`vfi [-h] [-a] [-b] [-c] [-n] [-p] [-r RCUT] [-s] [-t SYMPREC] [-v] [--debug] [--version] FILE`

```
**Argument**                            **Description**                                                                
| ------------------------------------| ----------------------------------------------------------------------------  |
| `FILE`                              | Input file containing structural data to process                              |
|                                     |                                                                               |
**Option**                              **Description**                                                                   |
| ----------------------------------- | ----------------------------------------------------------------------------- |
| `-h`, `--help`                      | Show help message and exit                                                    |
| `-a`, `--atoms`                     | Print atomic-scale info (ntotal, types, etc.)                                 |
| `-b`, `--bonds`                     | Print bonding info (total, species connectivity)                              |
| `-c`, `--cell`                      | Print unit cell parameters (a, b, c, volume, etc.)                            |
| `-n`, `--neighbors`                 | Print bonding info and nearest-neighbor data (recursive search)               |
| `-p`, `--primitive`                 | Reduce conventional cell to primitive if possible                             |
| `-r RCUT`, `--radius=RCUT`          | Search radius for considering atoms as bonded (default = `0.0 Å`)             |
| `-s`, `--save`                      | Save computed data to files (`.bonds`, `.atoms`, `.cell`) instead of printing |
| `-t SYMPREC`, `--tolerance=SYMPREC` | Precision in determining symmetry (default = `0.05 Å`)                        |
| `-v`                                | Increase verbosity                                                            |
| `--debug`                           | Print extensive diagnostic information (`-vvvv` equivalent)                   |
| `--version`                         | Show version number and exit                                                  |
| ----------------------------------- | ----------------------------------------------------------------------------- |
```

## Examples
```
# Print all structure information
vfi POSCAR > structure-data.nfo

# Reduce to primitive cell
vfi -p POSCAR

# Print bonding and cell data
vfi -bc POSCAR > bonding-data.nfo

# Save all computed data for sns2.vasp
vfi --debug --radius=2 --neighbors --save sns2.vasp --primitive

# Detailed run with tight symmetry tolerance
vfi -vvvvvnacrs 3 --tolerance=1e-3 mos2.contcar

# Moderate verbosity and custom tolerance
vfi -vvnacrs 3 mos2.contcar -t 0.01
vfi -gnacrs3 mos2.contcar -t0.1
```

## Output
If the --save option is used, files are generated automatically using the stoichiometry of the structure:
```
P1S2H2.bonds
P1S2H2.atoms
P1S2H2.cell
```

If the --primitive option is used, a single file is generated for the primitive unitcell, tagged with the original number of atoms:
```
Si215-primitive.vasp
``` 

## Directory structure
```
VaspFileInspector/
├── pyproject.toml
├── README.md
├── src/
│   └── vaspfileinspector/
│       ├── __init__.py
│       ├── common.py
│       ├── atoms.py
│       ├── lattice.py
│       ├── neighbors.py
│       └── reader.py
│       └── cli.py
```
