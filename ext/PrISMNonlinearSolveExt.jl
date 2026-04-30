module PrISMNonlinearSolveExt

using PrISM, NonlinearSolve
import PrISM: inverse!
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

# ===== NonlinearSolve.jl =========

function inverse!(mₖ::model1,
        robs::response,
        vars,
        alg_cache::nl_cache;
        params=default_params(model1),
        W=nothing,
        L=nothing,
        max_iters=30,
        χ2=1.0,
        response_fields=propertynames(robs),
        model_trans_utils=sigmoid_tf,
        response_trans_utils=nothing,
        mᵣ=nothing,
        reg_term=nothing,
        verbose::Union{Bool, Int}=true) where {
        model1 <: AbstractGeophyModel, response <: AbstractGeophyResponse}
    prec = eltype(mₖ.m)

    # n_vars = length(getfield(robs, first(response_fields)))
    n_resp = sum([length(getfield(robs, k)) for k in response_fields])

    if isnothing(response_trans_utils)
        response_trans_utils = NamedTuple{response_fields}(ntuple(i -> no_tf, length(response_fields)))
    end

    ks = Tuple([k for k in propertynames(mₖ) if k != :m])
    ps = Tuple([getfield(mₖ, k) for k in propertynames(mₖ) if k != :m])

    const_m = NamedTuple{ks}(ps)

    (W === nothing) && (W = prec.(I(n_resp)))
    (L === nothing) && (L = prec.(∂(length(mₖ.m))))

    model_type = typeof(mₖ).name.wrapper
    if isnothing(mᵣ)
        mᵣ = copy(mₖ)
        mᵣ.m .= zero(mₖ.m)
    end

    p = (model_type=model_type, m_const=const_m, model_trans_utils=model_trans_utils,
        response_trans_utils=response_trans_utils, vars=vars, params=params,
        response_fields=response_fields, W=W, μ=alg_cache.μ, r_obs=robs, L=L, mᵣ=mᵣ)

    prob = SciMLBase.NonlinearLeastSquaresProblem(
        construct_cost_function_for_nl_inv, model_trans_utils.itf.(mₖ.m), p)
    nlcache = init(prob, alg_cache.alg)
    iters = 1
    chi2 = 1e6

    resp_ = zero(robs)

    while iters <= max_iters
        m0 = merge((; m=model_trans_utils.tf.(nlcache.u)), const_m)
        model = from_nt(model_type, m0)
        forward!(resp_, model, vars, params)

        for k in response_fields
            broadcast!(getfield(response_trans_utils, k).tf, getfield(resp_, k), getfield(resp_, k))
        end

        chi2 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
            reduce(vcat, [getfield(robs, k) for k in response_fields]); W=W)
        do_verbose(iters, verbose) && println("iteration = $iters : data misfit => $chi2")

        if chi2 <= χ2 # check misfit condition
            break
        end

        step!(nlcache)
        mₖ.m .= model_trans_utils.tf.(nlcache.u)
        iters += 1
    end

    return return_code(chi2 <= χ2, (μ=alg_cache.μ,), mₖ, χ2, chi2)
end

# Not performant at the moment
function construct_cost_function_for_nl_inv(m, p)
    @unpack model_type, m_const, model_trans_utils, response_trans_utils,
    vars, params, response_fields, W, μ, r_obs, L, mᵣ = p

    m0 = merge((; m=model_trans_utils.tf.(m)), m_const)
    model = from_nt(model_type, m0)
    resp_ = forward(model, vars, params)

    for k in response_fields
        broadcast!(getfield(response_trans_utils, k).tf, getfield(resp_, k), getfield(resp_, k))
    end

    L1 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
        reduce(vcat, [getfield(r_obs, k) for k in response_fields]); W=W) * sqrt(size(W, 1))
    L2 = μ * norm(L * (model.m .- mᵣ.m))

    return [L1^2 + L2]
end

end
