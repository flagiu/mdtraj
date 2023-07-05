# MDtraj

This program computes statistical quantities over a Molecular Dynamics TRAJectory.

- Current limited to: monospecies system, constant number of particles, uniform timestep.

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
where version can be 'main' (stable version), 'dev' (developement version).

## Usage

Run helper message for instructions:
```bash
path-to/this-repo/version/bin/mdtraj -h
```

Folders python/ and shell/ contains some helper scripts to be used externally.

## Future development

- Add computation of ALTBC, structure factor, ISF.
- Add multi-species.
- Add logarithmic timestep.
- Add more input formats: GROMACS, ...?
