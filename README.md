# Hybrid Lubrication Granular Dynamics (HLGD)

HLGD is a simple program to simulate dense particle suspensions in shear flow (Stokesian). My intention here is to promote open science and possibly help other researchers implement something similar. The code uses OpenMP for parallelization, but the implementation is not optimzed. The neighbor search part is not clean in hindsight. Overall, the program runs fine on one node in supercomputers as of 2020-21.

To run the code, simply do the following:

1. Compile the code by `make`
2. Run the code by `./run`

The simulation parameters are mainly in `param.f90`.

### References:

> 1. Ge, Zhouyang, and Luca Brandt. "Implementation note on a minimal hybrid lubrication/granular dynamics model for dense suspensions." arXiv preprint arXiv:2005.12755 (2020).
> 2. Ge, Zhouyang, Raffaella Martone, Luca Brandt, and Mario Minale. "Irreversibility and rate dependence in sheared adhesive suspensions." Physical Review Fluids 6, no. 10 (2021): L101301.
