# Polarize

MATLAB tools for building periodic or finite-cluster polarizable systems and computing polarization energies in an AMOEBA/Thole-style model.

This branch is a refactored, minimal working version focused on a robust periodic workflow:

- read a VASP structure
- identify molecules in the unit cell
- build a supercell
- identify actual molecules in the supercell with a periodic bond graph
- determine which molecules are complete in the displayed supercell
- select reference/neighbor pairs by geometric relation
- assign charges to one or two molecules
- disable polarizability on charged molecules
- pass the resulting `polsys` into a polarization-energy calculator

## Current status

The current refactor supports:

- VASP structure import
- unit-cell molecular template construction
- graph-based site typing
- blind supercell construction
- periodic supercell molecule identification
- complete-molecule filtering in the displayed box
- same-stack / side-stack / skew neighbor selection
- pair-centered selection
- uniform charge assignment to one or two molecules
- turning off polarizability on charged molecules

This is intended to be a cleaner foundation before reintroducing more of the original codebase.

## Core idea

There are two distinct layers of molecule information in this workflow:

### 1. Unit-cell chemical template
Used for:
- local atom labels
- atom classes
- external charge-template mapping
- molecule-frame construction

### 2. Working supercell molecules
Used for:
- reference selection
- neighbor selection
- completeness checks
- active/charged molecule assignment

## Nonperiodic cell-list and MEX acceleration

Finite nonperiodic calculations now use a Cartesian **cell list** to accelerate cutoff-based neighbor search. Sites are binned into a 3D grid, nearby bins are searched through a precomputed stencil, and exact Cartesian distances are used to retain only pairs within `rcut`.

This spatial index is then used to build reusable geometric caches, including the **nonperiodic pair cache** and the **active row cache** used by the matrix-free iterative solver.

The heavy cache-building path is accelerated with a **MEX backend**, so the code can quickly:

- enumerate neighbor pairs within cutoff,
- assemble cached pair geometry once,
- and reuse that data during iterative solves.

In practice, this moves work out of the solver iterations and into a one-time preprocessing step, substantially reducing the cost of nonperiodic real-space calculations.
