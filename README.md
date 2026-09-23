# Project 2: 2D Poisson Equation Using Finite Differences

This project solves a two-dimensional Poisson equation on the unit square using the standard five-point finite difference method.

## Model problem

\[
-\Delta u = f
\]

on

\[
(0,1)\times(0,1),
\]

with zero Dirichlet boundary conditions.

The exact solution is chosen as

\[
u(x,y)=\sin(\pi x)\sin(\pi y),
\]

which gives

\[
f(x,y)=2\pi^2\sin(\pi x)\sin(\pi y).
\]

## Main numerical ideas

- 2D Cartesian grid
- five-point finite difference stencil
- sparse matrix assembly
- solution of a sparse linear system
- comparison with an exact solution
- grid-refinement and convergence study

## Files

- `main_poisson2d.m` — main simulation and visualization
- `poisson2d_fd.m` — sparse finite-difference solver
- `convergence_poisson2d.m` — convergence and runtime study
- `DAY4_NOTES.md` — Day 4 derivation and checklist

## Day 4

Run:

```matlab
main_poisson2d
```

Start with `N = 41`, inspect the numerical solution and error, then repeat with `N = 81`.

Do not run the full convergence study until the main solver and five-point stencil are understood.
