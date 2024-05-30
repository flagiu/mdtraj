[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# MDtraj

This program computes statistical quantities over a Molecular Dynamics TRAJectory.

- Current limited to:
	- constant number of particles
	- uniform timestep.
	- S(q) for cubic boxes only.

- Bugs to be corrected:
	- some empty files are created during 'altbc' and 'bond_order' (maybe?)
	- Check again the MSD and NGP.

## Requirements

- Main code: Basic C++ compiler and libraries.
- Python utilities: python3 with few libraries: numpy, scipy, matplotlib.

## Installation (Linux)

Clone or Download this repository, install the path and compile with make:
```bash
gh repo clone flagiu/mdtraj
cd mdtraj/version/
bash ./install_path.sh
make
```
where 'version' can be:
- 'notype' (mono-species);
- 'mixture' (mono- or multi-species);

The best maintained version is currently 'mixture'.

## Usage and examples

Run helper message for instructions:
```bash
path-to-this-repo/version/bin/mdtraj -h
```

The subfolders python/ and shell/ contain some utility scripts (to be used before/after the main program) for plotting or extra calculations. Default units are: Angstrom, picoseconds.

The subfolder examples/ contains some example of application to real or toy systems:
- (notype) LJ particles in different 3d cell shapes: cubic, orthorombic, triclinic.
- (both) Toy particles displaced along a 3d cubic lattice with small random gaussian noise.
- (mixture) Binary LJ mixture.
- (mixture) Short sample of a Phase-Change Heterostructure.

## Future development

- More input formats for multi-species.
- More analysis tools for multi-species.
- S(q) with non-cubic boxes?
- Add example with antimony?
- Add logarithmic timestep?
- Add more input formats? GROMACS, ...
- Add structural quasi-entropy? [Oganov,Valle,2008]
- Add crystalline clusters analysis?
