# ======================== Occam ========================

"""
`occam_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct occam_cache{T}
    μgrid::Vector{T}
end

# ======================== NonlinearSolve ========================

"""
`nl_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct nl_cache{T1, T2}
    alg::T1
    μ::T2
end

"""
    NonlinearAlg(; alg = LevenbergMarquardt, μ = 1.0)

returns `nl_cache` that specifies which non linear solver to use for the inverse problem

## Keyword Arguments

  - `alg`: `NonlinearSolve`[@ref] algorithm to be used, defaults to LevenbergMarquardt
  - `μ` : regularization weight
"""
function NonlinearAlg(; alg=LevenbergMarquardt, μ=1.0)
    return nl_cache(alg, μ)
end

# ======================== Optimization ========================

"""
`opt_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct opt_cache{T1, T2}
    alg::T1
    μ::T2
end
"""
    OptAlg(; alg = LBFGS, μ = 1.0, kwargs...)

returns `opt_cache` that specifies which `Optimization.jl` solver to use for the inverse problem

## Keyword Arguments

  - `alg`: `Optimization`[@ref] algorithm to be used, defaults to LBFGS
  - `μ` : regularization weight
"""
function OptAlg(; alg=LBFGS, μ=1.0)
    return opt_cache(alg, μ)
end
