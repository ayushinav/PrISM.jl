using Pkg
Pkg.activate(".")

using MT
using Test
using Distributions
using Turing


# two-phase (multi-rp) without MAL

mix = HS1962_plus()
m_mix = two_phase_modelType(SEO3, Sifre2014, mix)
m = multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)
m = multi_rp_modelType(m_mix, anharmonic, Nothing, Nothing)

# @inferred MT.multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)


T = [1200f0] #collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = [1.0f0]
Ch2o_m = [1000.0f0]
Cco2_m = [1000.0f0]
Ch2o_opx = [100.0f0]
Cco2_m = [10.0f0]
# ϕ = [0.3f0] #[0.1f0, 0.2f0]
ϕ = [0.1f0, 0.2f0]
P = [4f0]
ρ = [3300f0]
m_MAL = 0.6f0
m_GAL = [0.4, 0.6]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx, P, ρ, m_GAL, m_MAL)

m1 = m(ps_nt)

# todo : from_nt

# plus  : Float32[-3.2222052, -2.105757, -1.5998105, -1.103057, -0.5914401]
# minus :  Float32[2.5088632, 2.6440232, 2.7609413, 2.8631563, 2.9532573]

# m = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, MT.HS_plus_multi_phase())
# m = two_phase_modelType(SEO3, Sifre2014, mix)

model = m(ps_nt)
resp = forward(model, [])
err_resp = multi_rp_response(RockphyCond(0.01 .* abs.(resp.cond.σ)), 
            RockphyElastic(
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.Vp),
                0.01 .* abs.(resp.elastic.Vs),
            ), Nothing, Nothing)

m_mixdist = two_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, mix)

# m_mixdist = two_phase_modelDistributionType(
#     SEO3Distribution, Sifre2014Distribution, HS1962_plus())
m = multi_rp_modelDistributionType(
    typeof(m_mixdist), anharmonicDistribution, Nothing, Nothing)
ps_nt_dist = (;
    T=product_distribution([Uniform(1000., 1400.)]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], # ϕ=product_distribution([Uniform(0.01, 0.2)]),
    Ch2o_opx=product_distribution([Uniform(50.0, 150.0)]),
    # ϕ = product_distribution([Uniform(0.0f0, 0.1f0)]),
    ϕ = product_distribution([Uniform(0.0f0, 0.1f0), Uniform(0.0f0, 0.3f0)]),
    P = [4f0], ρ = [3300f0],
    m_MAL = [0.3f0], #product_distribution([Uniform(0.3f0, 0.7f0)])
    m_GAL = [0.4, 0.6]
    # m_GAL = product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])
)

mdist = m(ps_nt_dist);
@show typeof(mdist)
rdist = MT.multi_rp_responseDistribution(
    RockphyCondDistribution(normal_dist),
    RockphyElasticDistribution(normal_dist, normal_dist, normal_dist, normal_dist),
    Nothing, Nothing
)

m_cache = mcmc_cache(mdist, rdist, 50, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 20, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);


# ================================================================================================
# ================================================================================================

# ================================================================================================
# ================================================================================================

# two-phase (multi-rp) with MAL

mix = MAL(0.5)
m_mix = two_phase_modelType(SEO3, Sifre2014, mix)
m = multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)

# @inferred MT.multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)


T = [1200f0] #collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = [1.0f0]
Ch2o_m = [1000.0f0]
Cco2_m = [1000.0f0]
Ch2o_opx = [100.0f0]
Cco2_m = [10.0f0]
# ϕ = [0.3f0] #[0.1f0, 0.2f0]
ϕ = [0.1f0, 0.2f0]
P = [4f0]
ρ = [3300f0]
m_MAL = 0.6f0
m_GAL = [0.4, 0.6]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx, P, ρ, m_GAL, m_MAL)

m1 = m(ps_nt)

# todo : from_nt

# plus  : Float32[-3.2222052, -2.105757, -1.5998105, -1.103057, -0.5914401]
# minus :  Float32[2.5088632, 2.6440232, 2.7609413, 2.8631563, 2.9532573]

# m = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, MT.HS_plus_multi_phase())
# m = two_phase_modelType(SEO3, Sifre2014, mix)

model = m(ps_nt)
resp = forward(model, [])
err_resp = multi_rp_response(RockphyCond(0.01 .* abs.(resp.cond.σ)), 
            RockphyElastic(
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.Vp),
                0.01 .* abs.(resp.elastic.Vs),
            ), Nothing, Nothing)

m_mixdist = two_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, MAL([0.5]))

# m_mixdist = two_phase_modelDistributionType(
#     SEO3Distribution, Sifre2014Distribution, HS1962_plus())
m = multi_rp_modelDistributionType(
    typeof(m_mixdist), anharmonicDistribution, Nothing, Nothing)
ps_nt_dist = (;
    T=product_distribution([Uniform(1000., 1400.)]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], # ϕ=product_distribution([Uniform(0.01, 0.2)]),
    Ch2o_opx=product_distribution([Uniform(50.0, 150.0)]),
    # ϕ = product_distribution([Uniform(0.0f0, 0.1f0)]),
    ϕ = product_distribution([Uniform(0.0f0, 0.1f0), Uniform(0.0f0, 0.3f0)]),
    P = [4f0], ρ = [3300f0],
    m_MAL = [0.3f0], #product_distribution([Uniform(0.3f0, 0.7f0)])
    m_GAL = [0.4, 0.6]
    # m_GAL = product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])
)

mdist = m(ps_nt_dist);
@show typeof(mdist)
rdist = MT.multi_rp_responseDistribution(
    RockphyCondDistribution(normal_dist),
    RockphyElasticDistribution(normal_dist, normal_dist, normal_dist, normal_dist),
    Nothing, Nothing
)

m_cache = mcmc_cache(mdist, rdist, 50, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 20, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);


# ================================================================================================
# ================================================================================================

# ================================================================================================
# ================================================================================================


# multi-phase (multi-rp) without GAL


mix = HS_plus_multi_phase()
# mix = GAL([0.4, 0.6])
m_mix = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, mix)
# m_mix = two_phase_modelType(SEO3, Sifre2014, mix)

# m_mix = two_phase_modelType(SEO3, Sifre2014, mix)
m = multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)

# @inferred MT.multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)


T = [1200f0] #collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = [1.0f0]
Ch2o_m = [1000.0f0]
Cco2_m = [1000.0f0]
Ch2o_opx = [100.0f0]
Cco2_m = [10.0f0]
# ϕ = [0.3f0] #[0.1f0, 0.2f0]
ϕ = [0.1f0, 0.2f0]
P = [4f0]
ρ = [3300f0]
m_MAL = 0.6f0
m_GAL = [0.4, 0.6]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx, P, ρ, m_GAL, m_MAL)

# todo : from_nt

# plus  : Float32[-3.2222052, -2.105757, -1.5998105, -1.103057, -0.5914401]
# minus :  Float32[2.5088632, 2.6440232, 2.7609413, 2.8631563, 2.9532573]

model = m(ps_nt)
resp = forward(model, [])
err_resp = multi_rp_response(RockphyCond(0.01 .* abs.(resp.cond.σ)), 
            RockphyElastic(
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.Vp),
                0.01 .* abs.(resp.elastic.Vs),
            ), Nothing, Nothing)

m_mixdist = MT.multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, Zhang2012Distribution, mix)
# m_mixdist = MT.multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, Zhang2012Distribution, GAL(product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])))

# m_mixdist = two_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, MAL([0.3f0])) #[Uniform(0.3f0, 0.7f0)])))

# m_mixdist = two_phase_modelDistributionType(
#     SEO3Distribution, Sifre2014Distribution, HS1962_plus())
m = multi_rp_modelDistributionType(
    typeof(m_mixdist), anharmonicDistribution, Nothing, Nothing)
ps_nt_dist = (;
    T=product_distribution([Uniform(1000., 1400.)]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], # ϕ=product_distribution([Uniform(0.01, 0.2)]),
    Ch2o_opx=product_distribution([Uniform(50.0, 150.0)]),
    # ϕ = product_distribution([Uniform(0.0f0, 0.1f0)]),
    ϕ = product_distribution([Uniform(0.0f0, 0.1f0), Uniform(0.0f0, 0.3f0)]),
    P = [4f0], ρ = [3300f0],
    m_MAL = [0.3f0], #product_distribution([Uniform(0.3f0, 0.7f0)])
    m_GAL = [0.4, 0.6]
    # m_GAL = product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])
)

mdist = m(ps_nt_dist);
@show typeof(mdist)
rdist = MT.multi_rp_responseDistribution(
    RockphyCondDistribution(normal_dist),
    RockphyElasticDistribution(normal_dist, normal_dist, normal_dist, normal_dist),
    Nothing, Nothing
)

m_cache = mcmc_cache(mdist, rdist, 50, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 20, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);


# ================================================================================================
# ================================================================================================

# ================================================================================================
# ================================================================================================


# multi-phase (multi-rp) with GAL


# mix = HS_plus_multi_phase()
mix = GAL([0.4, 0.6])
m_mix = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, mix)
# m_mix = two_phase_modelType(SEO3, Sifre2014, mix)

# m_mix = two_phase_modelType(SEO3, Sifre2014, mix)
m = multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)

# @inferred MT.multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)


T = [1200f0] #collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = [1.0f0]
Ch2o_m = [1000.0f0]
Cco2_m = [1000.0f0]
Ch2o_opx = [100.0f0]
Cco2_m = [10.0f0]
# ϕ = [0.3f0] #[0.1f0, 0.2f0]
ϕ = [0.1f0, 0.2f0]
P = [4f0]
ρ = [3300f0]
m_MAL = 0.6f0
m_GAL = [0.4, 0.6]

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m, Ch2o_opx, P, ρ, m_GAL, m_MAL)

# todo : from_nt

# plus  : Float32[-3.2222052, -2.105757, -1.5998105, -1.103057, -0.5914401]
# minus :  Float32[2.5088632, 2.6440232, 2.7609413, 2.8631563, 2.9532573]

model = m(ps_nt)
resp = forward(model, [])
err_resp = multi_rp_response(RockphyCond(0.01 .* abs.(resp.cond.σ)), 
            RockphyElastic(
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.G),
                0.01 .* abs.(resp.elastic.Vp),
                0.01 .* abs.(resp.elastic.Vs),
            ), Nothing, Nothing)

m_mixdist = MT.multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, Zhang2012Distribution, mix)
# m_mixdist = MT.multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, Zhang2012Distribution, GAL(product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])))

# m_mixdist = two_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, MAL([0.3f0])) #[Uniform(0.3f0, 0.7f0)])))

# m_mixdist = two_phase_modelDistributionType(
#     SEO3Distribution, Sifre2014Distribution, HS1962_plus())
m = multi_rp_modelDistributionType(
    typeof(m_mixdist), anharmonicDistribution, Nothing, Nothing)
ps_nt_dist = (;
    T=product_distribution([Uniform(1000., 1400.)]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], # ϕ=product_distribution([Uniform(0.01, 0.2)]),
    Ch2o_opx=product_distribution([Uniform(50.0, 150.0)]),
    # ϕ = product_distribution([Uniform(0.0f0, 0.1f0)]),
    ϕ = product_distribution([Uniform(0.0f0, 0.1f0), Uniform(0.0f0, 0.3f0)]),
    P = [4f0], ρ = [3300f0],
    m_MAL = [0.3f0], #product_distribution([Uniform(0.3f0, 0.7f0)])
    m_GAL = [0.4, 0.6]
    # m_GAL = product_distribution([Uniform(0.3f0, 0.7f0), Uniform(0.3f0, 0.7f0)])
)

mdist = m(ps_nt_dist);
@show typeof(mdist)
rdist = MT.multi_rp_responseDistribution(
    RockphyCondDistribution(normal_dist),
    RockphyElasticDistribution(normal_dist, normal_dist, normal_dist, normal_dist),
    Nothing, Nothing
)

m_cache = mcmc_cache(mdist, rdist, 50, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 20, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);


