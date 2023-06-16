# MDtraj

This program computes statistical quantities over a Molecular Dynamics TRAJectory.

Current limited to: monospecies system, constant number of particles, uniform timestep.

Bugs to be corrected: exact value of q_l order parameters.

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
