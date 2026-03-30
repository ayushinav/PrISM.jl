# PrISM.jl

# Porosity.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ayushinav.github.io/PrISM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ayushinav.github.io/PrISM.jl/dev/)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/ayushinav/PrISM.jl/Tests.yml)
[![codecov](https://codecov.io/gh/ayushinav/PrISM.jl/graph/badge.svg?token=VQM6W3DUI4)](https://codecov.io/gh/ayushinav/PrISM.jl)

`PrISM.jl` is supposed to be a high performance code for doing forward and inverse modeling in geophysics using julia. We hope to write the code structure such that any other geophysical survey can also be used and we can tend towards a joint forward and inverse modeling library.

## Forward modeling

While forward modeling typically requires solving a PDE obtained using the quasi-static approximation, in 1D, we are fortunate to have the solutions in an analytical form. Currently, this is what is supported.

Supported methods:

  - 1D Magnetotellurics (MT)
  - 1D Rayleigh waves
  - 1D Love waves
  - 1D Direct Current (DC) Resistivity

## Inverse modeling

No surprises here that we are almost always trying to solve for an under-determined system.

Deterministic schemes supported:

  - Occam
  - Nonlinear schemes using `NonlinearSolve.jl`
  - Nonlinear schemes using `Optimization.jl`

Probabilistic schemes supported:

  - MCMC with fixed grids
  - MCMC with flexible grids
  - RTO-TKO
