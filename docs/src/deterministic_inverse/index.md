# Deterministic Inversion

## Brief introduction

Inverse problems in geophysics are notoriously ill-posed with non-unique solutions. This page goes over some essentials you might require for performing deterministic inversion of the geophysical datasets.

To be on the same page, the goal is to solve the following optimization problem : 

```math
\argmin_{m} \; [f(m) - d]^T C_d [f(m) - d] \; +  \; \mu \; (m - m_0)^T C_m (m - m_0)
```

where $m$ denotes the model parameters, $d$ is the data to be inverted for, $C_d$ is the data covariance matrix, composed of errors in data and $C_m$ is model covariance matrix, constructed using prior assumptions about the model, $\mu$ is the regularization coefficient, and $m_0$ is some model you want to bias you results towards.

The first term in the above equation measures the misfit between the data and the synthetic response that would have been observed because of $m$. The second term is used to stabilize and/or constrain the solutions better. Very frequently, $C_m$ is constructed using the derivative matrix but it can be anything, e.g. a weight matrix to bias some parts of the models towards $m_0$.

Obviously, there are multiple ways to solve this non-linear equation. We provide the following few natively:

 * Occam1D : [Occam's inversion: a practical algorithm for generating smooth models from electromagnetic sounding data](https://marineemlab.ucsd.edu/steve/bio/Occam1D.pdf)

* [Optimisers](https://docs.sciml.ai/NonlinearSolve/stable/native/solvers/) from `NonlinearSolve.jl`. Popular ones in geophysical community include:
     * Gauss Newton
     * Levenberg Marquardt
     * Newton Raphson

* [Optimisers](https://docs.sciml.ai/Optimization/stable/optimization_packages/optim/#Methods) from `OptimizationOptim.jl` via `Optimization.jl`. Popular ones in geophysical community include:
    * Conjugate Gradient
    * Gradient Descent
    * LBFGS
    * Simulated Annealing
    * Particle Swarm

All the inversion capabilities are accessed using `inverse!` function, which we cover in the next pages.

!!! note
    Do note that the inversion only takes care of the parameter denoted by `m`, even though there might be other parameters, e.g. for Rayleigh wave models, we invert of shear wave velocity even though the model also requires p-wave velocities and densities.