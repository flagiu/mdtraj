[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# MDtraj

This program computes statistical quantities over a Molecular Dynamics TRAJectory.

## Disclaimer

This code is not bug-free, nor extensively and accurately tested by the author.

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

Most of the input comes from the **command line**, except some features requiring **parameters files** (e.g. see the -rcut command).

The subfolders python/ and shell/ contain some **utility scripts** (to be used before/after the main program) for plotting or extra calculations. Default units are: Angstrom, picoseconds.

The subfolder examples/ contains some **example** of application to real or toy systems:
- (notype) LJ particles in different 3d cell shapes: cubic, orthorombic, triclinic.
- (both) Toy particles displaced along a 3d cubic lattice with small random gaussian noise.
- (mixture) Binary LJ mixture.
- (mixture) Short sample of a Phase-Change Heterostructure.
- (mixture) Water molecules.
- (mixture) Logarithmic timestep.

**Some limitations** of the current version are:
- different frames should have the same number of particles
- S(q) and S(q,t) make sense for cubic boxes only; and they deal with particles as monospecies, ignoring their type.

## Future development ideas
- S(q),S(q,t) with non-cubic boxes?
- S(q),S(q,t) with mixtures? How? Weighted by mass?
- Add more input formats? GROMACS, ...
- Add structural quasi-entropy? [Oganov,Valle,2008]
- Add crystalline clusters analysis like pyscal?

## Acknowledgements
The development of this code was supported by the ICSC - *Centro Nazionale di Ricerca in High Performance Computing, Big Data, and Quantum Computing* funded by the European Union - NextGenerationEU.
