@testitem "ForwardDiff AD compatibility" tags=[:forwarddiff] begin
    using PrISM, Distributions, Turing, ForwardDiff, LinearAlgebra, Test
    using DynamicPPL.TestUtils.AD: run_ad, ADResult
    using DifferentiationInterface: AutoForwardDiff
    model_types = [MTModel, RWModel, LWModel]
    modelD_types = [MTModelDistribution, RWModelDistribution, LWModelDistribution]
    respD = [MTResponseDistribution(normal_dist, normal_dist),
        SurfaceWaveResponseDistribution(normal_dist),
        SurfaceWaveResponseDistribution(normal_dist)]

    true_models = [(; m=[2.0, 1.0, 2.0], h=[1000.0, 1000.0]),
        (; m=[4000.0, 3500.0, 4000.0] ./ 1e3, h=[1000.0, 1000.0],
            vp=[7000.0, 7000.0, 7000.0] ./ 1e3, ρ=[2500.0, 2500.0, 2500.0] ./ 1e3),
        (; m=[3500.0, 3600.0, 3800.0] ./ 1e3, h=[1000.0, 1000.0],
            ρ=[2500.0, 2500.0, 2500.0] ./ 1e3)]

    vars = [10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(0:0.1:3), 10.0 .^ collect(-1:0.1:1)]

    adtype = AutoForwardDiff()

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
end
