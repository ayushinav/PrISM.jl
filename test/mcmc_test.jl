@testitem "fixed discretization" tags = [:mcmc] begin
    using Distributions, Turing, LinearAlgebra, Pigeons

    model_types = [MTModel, RWModel]
    modelD_types = [MTModelDistribution, RWModelDistribution, LWModelDistribution]
    respD = [MTResponseDistribution(normal_dist, normal_dist),
        SurfaceWaveResponseDistribution(normal_dist),
        SurfaceWaveResponseDistribution(normal_dist)]

    true_models = [(; m=[2.0, 1.0, 2.0], h=[1000.0, 1000.0]),
        (; m=[3500.0, 3000.0, 3500.0] ./ 1e3, h=[100.0, 100.0],
            vp=[7000.0, 7000.0, 7000.0] ./ 1e3, ρ=[2500.0, 2500.0] ./ 1e3),
        (; m=[3500.0, 3000.0, 3500.0] ./ 1e3, h=[100.0, 100.0], ρ=[2500.0, 2500.0] ./ 1e3)]

    vars = [10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(0:0.1:3), 10.0 .^ collect(0:0.1:3)]

    samplers = [NUTS, SliceSampler]
    @testset "$(model_types[ik]) : $(sampler_)" for ik in eachindex(model_types),
        sampler_ in samplers

        model = from_nt(model_types[ik], true_models[ik])
        resp = forward(model, vars[ik])

        err_resp = copy(resp)
        for k in propertynames(resp)
            setproperty!(err_resp, k, getproperty(err_resp, k) .* 0.05)
        end

        mdist = from_nt(modelD_types[ik],
            (; true_models[ik]...,
                m=product_distribution([Uniform(mi * 0.8, mi * 1.2) for mi in model.m])))

        n_samples = 50
        m_cache = mcmc_cache(mdist, respD[ik], n_samples, sampler_())

        mcmc_chain = stochastic_inverse(resp, err_resp, vars[ik], m_cache)

        model_list = get_model_list(mcmc_chain, mdist)
        wvec = vcat([getfield(err_resp, k) for k in propertynames(resp)]...)
        W = diagm(inv.(wvec)) .^ 2

        err = zeros(n_samples)
        for idx in 1:n_samples
            m_model = model_list[idx]
            resp_model = forward(m_model, vars[ik])

            err[idx] = χ²(
                reduce(vcat, [getfield(resp_model, k) for k in propertynames(resp)]),
                reduce(vcat, [getfield(resp, k) for k in propertynames(resp)]); W=W)
        end
        @show err

        @test sum(err[41:50]) <= (sum(err[1:10]) + 5)
    end
end

@testitem "variable discretization" tags = [:mcmc] begin
    using Distributions, Turing, LinearAlgebra, Pigeons

    model_types = [MTModel, RWModel, LWModel]
    modelD_types = [MTModelDistribution, RWModelDistribution, LWModelDistribution]
    respD = [MTResponseDistribution(normal_dist, normal_dist),
        SurfaceWaveResponseDistribution(normal_dist),
        SurfaceWaveResponseDistribution(normal_dist)]

    true_models = [(; m=[2.0, 1.0, 2.0], h=[1000.0, 1000.0]),
        (; m=[3500.0, 3000.0, 3500.0] ./ 1e3, h=[100.0, 100.0],
            vp=[7000.0, 7000.0, 7000.0] ./ 1e3, ρ=[2500.0, 2500.0] ./ 1e3),
        (; m=[3500.0, 3000.0, 3500.0] ./ 1e3, h=[100.0, 100.0], ρ=[2500.0, 2500.0] ./ 1e3)]

    vars = [10.0 .^ collect(-3:0.1:1), 10.0 .^ collect(0:0.1:3), 10.0 .^ collect(0:0.1:3)]

    samplers = [NUTS, SliceSampler]
    @testset "$(model_types[ik]) : $(sampler_)" for ik in eachindex(model_types),
        sampler_ in samplers

        model = from_nt(model_types[ik], true_models[ik])
        resp = forward(model, vars[ik])

        err_resp = copy(resp)
        for k in propertynames(resp)
            setproperty!(err_resp, k, getproperty(err_resp, k) .* 0.05)
        end

        mdist = from_nt(modelD_types[ik],
            (; true_models[ik]...,
                m=product_distribution([Uniform(mi * 0.8, mi * 1.2) for mi in model.m]),
                h=product_distribution([Uniform(hi * 0.8, hi * 1.2) for hi in model.h])))

        n_samples = 50
        m_cache = mcmc_cache(mdist, respD[ik], n_samples, sampler_())

        mcmc_chain = stochastic_inverse(resp, err_resp, vars[ik], m_cache)

        model_list = get_model_list(mcmc_chain, mdist)
        wvec = vcat([getfield(err_resp, k) for k in propertynames(resp)]...)
        W = diagm(inv.(wvec)) .^ 2

        err = zeros(n_samples)
        for idx in 1:n_samples
            m_model = model_list[idx]
            resp_model = forward(m_model, vars[ik])

            err[idx] = χ²(
                reduce(vcat, [getfield(resp_model, k) for k in propertynames(resp)]),
                reduce(vcat, [getfield(resp, k) for k in propertynames(resp)]); W=W)
        end
        @show err

        @test sum(err[41:50]) <= (sum(err[1:10]) + 5)
    end
end
