using Pkg
Pkg.activate("pkgs/MT.jl/wip/")

using MT
using Distributions
using LinearAlgebra
using Turing
import Base: show



m1 = SEO3(1200.0 + 273)
m2 = Ni2011(1200.0 + 273, 4000)

forward(m1, [])
forward(m2, [])

mix1 = construct_mixing_models(
    [1300.0 + 273, 4000.0],
    [:T, :Ch2o_m],
    [0.2],
    [SEO3, Ni2011],
    HS1962_plus,
)

r_obs = forward(mix1, [])

# 1200 - 1500
#  0 - 500 ppm water

m_dist = RockphyModelDistribution(
    # MvNormal([1000., 0.3], diagm([200., 0.0001])),
    Product([Uniform(273 + 1000.0, 1800.0 + 273), Uniform(3991, 4001)]), #, Uniform(1991, 2001)]),
    [:T, :Ch2o_m],
    # [0.8, 0.2],
    Product([Uniform(0.1, 0.3)]),
    [SEO3, Ni2011],
    HS1962_plus,
)

r_dist = RockphyResponseDistribution((x, y) -> MultivariateNormal(x, y))

m_cache = mcmc_cache(m_dist, r_dist, 5_000, NUTS())

# r_obs = MT.RockphyCond([1e7])
err_resp = MT.RockphyCond([0.001])

# dist_check = MvNormal([0.4], [0.00001])
# samps = rand(dist_check, 100_000);

# mean(samps)
# std(samps)

samples = stochastic_inverse(r_obs, err_resp, [], m_cache) #, trans_utils=(m=log_tf,))


# samples.value

using Plots
Plots.plot(Ts .- 273, conds, yscale = :log10)
# r_mt = MTResponse(rand(100), rand(100))
# MT.inverse(r_mt)

# type_ = MT.inverse(r_obs, abstract = true)

# type_(
#     [1000., 0.1],
#     [:T, :Ch2o_m],
#     [0.8, 0.2],
#     [MT.SEO3, MT.Ni2011],
#     MT.HS_plus
# )
using StatsPlots
plot(samples)


# ======= occam 1d :nomelt ======

using JLD2

st = "1_yx"
kk = jldopen("nomelt/filter_data/" * st * ".jld2")
r_obs = kk["r_obs4"];
err_resp = kk["resp_err4"];
ω = kk["ω"];

idxs = sortperm(ω);
r_obs.ρₐ .= r_obs.ρₐ[idxs];
r_obs.ϕ .= r_obs.ϕ[idxs];
err_resp.ρₐ .= err_resp.ρₐ[idxs];
err_resp.ϕ .= err_resp.ϕ[idxs];
ω .= ω[idxs];


new_robs = MTResponse(r_obs.ρₐ[5:end], r_obs.ϕ[5:end])
new_err_resp = MTResponse(err_resp.ρₐ[5:end], err_resp.ϕ[5:end])
new_err_resp.ϕ .= asin(0.05) * 180 / π
new_err_resp.ρₐ .= 0.1 .* new_robs.ρₐ
new_ω = ω[5:end]

W = diagm(inv.([new_err_resp.ρₐ..., new_err_resp.ϕ...]))

z = 10.0f0 .^ collect(range(-1, 6, length = 200))
m_test = MTModel(fill(100.0, length(z)), diff(z))

inverse!(m_test, new_robs, new_ω, MT.Occam(μgrid = [1e-2, 1e6]), W = W, max_iters = 50)

plot_model(m_test, label = false, yscale = :identity, ylim = (1, 1e6))

plt_resp = prepare_plot(new_robs, new_ω, new_err_resp)

plt_resp = prepare_plot(r_obs, ω, new_err_resp)

resp_occam = forward(m_test, new_ω)
prepare_plot!(plt_resp, resp_occam, new_ω, plt_style = :plot)
plot_response(plt_resp)
