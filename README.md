# MDtraj

This program computes statistical quantities over a Molecular Dynamics TRAJectory.

- Current limited to: constant number of particles, uniform timestep.

- Under development: multi-species system for .xyz input files.

- Bugs to be corrected: q_l order parameters seems to be offset by ~sqrt(2).

## Requirements

Basic C++ compiler and libraries.

## Installation (Linux)

Clone or Download this repository and compile with make:
```bash
gh repo clone flagiu/this-repo
cd this-repo/version
make
```
where version can be 'notype' (monospecies), 'mixture' (multi-species), with or without the suffix '_dev' (stable/developement version).

## Usage

Run helper message for instructions:
```bash
path-to/this-repo/version/bin/mdtraj -h
```

The subfolders python/ and shell/ contains some utility scripts, to be used externally.

## Future development

- Add computation of structure factor, ISF.
- Add logarithmic timestep.
- Add more input formats: GROMACS, ...?
