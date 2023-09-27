# MDtraj

This program computes statistical quantities over a Molecular Dynamics TRAJectory.

- Current limited to:
	- constant number of particles
	- uniform timestep.
	- S(q) for cubic boxes only (?).

- Under development:
	- multi-species system for .xyz, CP2K-.xyz, lammpstrj input files.

- Bugs to be corrected:
	- q_l order parameters seems to be offset by ~sqrt(2).
	- g(r) is <1 for large r in triclinic boxes.
	- Check again the NGP.
	- S(q,0) is wrong by some period-dependent factor.

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
- 'notype', 'notype_dev' (mono-species);
- 'mixture', 'mixture_dev' (multi-species);
the '_dev' versions are under developement.

## Usage and examples

Run helper message for instructions:
```bash
path-to-this-repo/version/bin/mdtraj -h
```

The subfolders python/ and shell/ contain some utility scripts (to be used before/after the main program) for plotting or extra calculations. Default units are: Angstrom, picoseconds.

The subfolder examples/ contains some example of application to real or toy systems:
- LJ particles in different 3d cell shapes: cubic, orthorombic, triclinic.
- Toy particles displaced along a 3d cubic lattice with small random gaussian noise.

## Future development

- Multi-species for more input formats.
- S(q) with non-cubic boxes?
- Add example with antimony?
- Add computation of ISF.
- Add logarithmic timestep?
- Add more input formats? GROMACS, ...
