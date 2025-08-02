"""
    stochastic_inverse(
        r_obs::response,
        err_resp::response,
        vars,
        alg_cache::mcmc_cache;
        model_trans_utils::NamedTuple = (m = lin_tf, h = lin_tf)
        )

function to perform sampling

### Returns

    `AbstractMCMC.jl Chain` containing the vectors used for performing samping. Note that all the variables will be named `var_inf`, for variable inferred. 
    The vector in each sample will be all the fields of the model being inferred on concatenated together.

### Variables

  - `r_obs`: `response` that needs to inverted for
  - `err_resp`: `response` variable containing the errors associated with observed response
  - `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
  - `alg_cache`: to tell the compiler what type of stochastic inversion method is to be used
  - `model_trans_utils`: A named tuple containing `transform_utils` for the fields of model that need to be scaled/modified. If not provided for any `model` field, the field won't be modified.
"""
function stochastic_inverse2(r_obs::resp1, err_resp::resp2, vars, alg_cache::mcmc_cache;
        model_trans_utils::NamedTuple=(m=lin_tf, h=lin_tf), # need to take care of this
        response_trans_utils::NamedTuple=(; ρₐ=lin_tf, ϕ=lin_tf),
        kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}
    model_fields = []
    modelD = []
    const_data = []

    # segregate the constants and the Distribution parts of the alg_cache

    for k in fieldnames(typeof(alg_cache.apriori)) # make it properynames and make alg_cache.apriori a NamedTuple
        if typeof(getfield(alg_cache.apriori, k)) <: Distribution # getfield will be replaced by getproperty
            push!(model_fields, k)
            push!(modelD, getfield(alg_cache.apriori, k))
            push!(const_data, rand(getfield(alg_cache.apriori, k)))
        else
            push!(const_data, getfield(alg_cache.apriori, k))
        end
    end

    response_fields = Symbol.([])
    for k in fieldnames(typeof(alg_cache.likelihood)) # similarly, here it will be propertynames for likelihood being a NamedTuple
        if typeof(getfield(alg_cache.likelihood, k)) <: Function
            push!(response_fields, k)
        end
    end

    # putting trans_utils together for all the fields

    trans_utils_arr = []
    for k in fieldnames(typeof(alg_cache.apriori))
        if k in keys(model_trans_utils)
            push!(trans_utils_arr, model_trans_utils[k])
        else
            push!(trans_utils_arr, lin_tf)
        end
    end

    transf_utils = (; zip([fieldnames(typeof(alg_cache.apriori))...], trans_utils_arr)...) # NamedTuple for trans_utils and defaults

    m_sample = MT.inverse(r_obs; abstract=true)(const_data...)

    # resp_s1 = forward(m_sample, vars, response_trans_utils)
    resp_type = typeof(r_obs).name.wrapper

    # resp_sample = resp_type(
    #     [Array{Real, ndims(getfield(r_obs, k))}copy((getproperty(r_obs, k)))
    #     for k in propertynames(r_obs)]...
    # )

    

    resp_sample = resp_type(
        [DiffCache(copy((getproperty(r_obs, k))), length(getproperty(r_obs, k)))
        for k in propertynames(r_obs)]...
    )

    robs = (;
        zip([fieldnames(typeof(r_obs))...],
            [getfield(r_obs, k) for k in fieldnames(typeof(r_obs))])...)

    mcmc_model = mcmc_turing2(m_sample, resp_sample, const_data, vars, robs, # ::NamedTuple
        err_resp, # ::response
        alg_cache.apriori, # ::NamedTuple
        alg_cache.likelihood; # ::responseDistribution
        response_fields=Symbol.(response_fields), model_fields=Symbol.(model_fields),
        model_trans_utils=transf_utils, response_trans_utils=response_trans_utils)

    if typeof(alg_cache.sampler) <: Turing.AdvancedVI.VariationalInference
        return vi(mcmc_model, alg_cache.sampler)
    else
        # ttfx
        # Turing.sample(mcmc_model, alg_cache.sampler, 1; verbose=false, kwargs...)
        
        # return Turing.sample(
        #     mcmc_model, alg_cache.sampler, alg_cache.n_samples; verbose=false, kwargs...)
        
        return Turing.sample(
            mcmc_model, alg_cache.sampler, alg_cache.n_samples; verbose=false, kwargs...)
    end
end

"""
    @model function mcmc_turing2(
        m_sample::model,
        vars,
        r_obs::NamedTuple,
        err_resp::MTResponse,
        mDist::mdist,
        rDist::rdist;
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(rDist))],
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mDist))],
        trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
        ) where {model <: AbstractModel, mdist <: AbstractModelDistribution, rdist <: AbstractResponseDistribution}

makes a `Turing.jl` model to perform MCMC sampling

### Variables:

  - `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
  - `r_obs`: named tuple containing the observed data, with the same keys as the fields in the corresponding `response`
  - `err_resp`: `response` variable that contains the errors
  - `mDist`: any subtype of [`AbstractModelDistribution`](@ref) contains the apriori information
  - `rDist`: any subtype of [`AbstractResponseDistribution`](@ref) contains the likelihood information

### Keyword/optional arguments

  - `response_fields`:  which fields in `response` to invert for
  - `model_fields`: fields in `model` to draw inference on
  - `trans_utils`: to transform the model field variables to and from computational (inference) domain
"""
@model function mcmc_turing2(m_sample::model,
        r_sample :: resp1,
        const_data,
        vars,
        r_obs::NamedTuple,
        err_resp::response,
        mDist::mdist,
        rDist::rdist;
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(rDist))],
        model_fields::Vector{Symbol}=[k for k in fieldnames(typeof(mDist))],
        model_trans_utils::NamedTuple=(m=lin_tf, h=lin_tf),
        response_trans_utils::NamedTuple=(ρₐ=lin_tf, ϕ=lin_tf)) where {
        model <: AbstractModel, response <: AbstractResponse, resp1 <: AbstractResponse,
        mdist <: AbstractModelDistribution, rdist <: AbstractResponseDistribution}
    m0 = (; zip([propertynames(mDist)...], const_data)...)

    for k in propertynames(mDist)
        if k in model_fields
            m0[k] ~ getproperty(mDist, k)
        end
    end

    m_sample = typeof(m_sample)([broadcast(getproperty(model_trans_utils[k], :tf), m0[k])
                                 for k in propertynames(mDist)]...)

    for k in propertynames(mDist)
        setfield!(m_sample, k, broadcast(getproperty(model_trans_utils[k], :tf), m0[k]))
    end

    r_sample_ = get_tmp2(r_sample, m0[model_fields[1]])

    # r_sample = forward(m_sample, vars, response_trans_utils)
    forward!(r_sample_, m_sample, vars, response_trans_utils)

    for k in response_fields
        r_obs[k] ~ getfield(rDist, k)(getfield(r_sample_, k), getfield(err_resp, k) .^ 2)
    end
end


function get_tmp2(resp::R, x::T) where {T, R <:AbstractGeophyResponse}
    resp_type = typeof(resp).name.wrapper

    resp_type(
        [get_tmp((getproperty(resp, k)), x)
        for k in propertynames(resp)]...
    )
end