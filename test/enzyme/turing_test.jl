using PrISM, Distributions, Turing, Enzyme, LinearAlgebra, Test
using DynamicPPL.TestUtils.AD: run_ad, ADResult
using DifferentiationInterface: AutoEnzyme

model_types = [MTModel, RWModel, LWModel]
modelD_types = [MTModelDistribution, RWModelDistribution, LWModelDistribution]
respD = [MTResponseDistribution(normal_dist, normal_dist),
    SurfaceWaveResponseDistribution(normal_dist),
    SurfaceWaveResponseDistribution(normal_dist)]

true_models = ((; m=randn(50) .* 1.0 .+ 2, h=fill(100.0, 49)),
    (; m=rand(50) .* 20e-1 .+ 3.0, h=fill(100.0, 49), vp=fill(7.5, 50), ρ=fill(2.5, 50)),
    (; m=rand(50) .* 20e-1 .+ 3, h=fill(100.0, 49), ρ=fill(2.5, 50)),
    (; m=randn(50) .* 0.01 .+ 2, h=fill(100.0, 49)))

vars = [10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(0:0.1:3), 10.0 .^ collect(-1:0.1:1)]

adtype = AutoEnzyme(; mode=Enzyme.set_runtime_activity(Reverse))

@testset "$(model_types[ik])" for ik in eachindex(model_types)
    model = from_nt(model_types[ik], true_models[ik])
    resp = forward(model, vars[ik])

    err_resp = copy(resp)
    for k in propertynames(resp)
        setproperty!(err_resp, k, getproperty(err_resp, k) .* 0.05)
    end

    mdist = from_nt(modelD_types[ik],
        (; true_models[ik]...,
            m=product_distribution([Uniform(mi * 0.8, mi * 1.2) for mi in model.m])))

    m_cache = mcmc_cache(mdist, respD[ik])
    dppl_object, _ = get_stochastic_inverse_model(resp, err_resp, vars[ik], m_cache)

    result = run_ad(dppl_object, adtype; atol=1e-6, params=rand(mdist.m))
    @test result isa ADResult
end
