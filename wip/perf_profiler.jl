using Pkg
Pkg.activate("wip/.")

using MT
using Profile
using LinearAlgebra
using Turing
using Distributions
using Plots

Turing.setrdcache(true)

m_test = MTModel([100.0, 10.0, 1000.0], [1e3, 1e3])
f = 10 .^ range(-4, stop = 1, length = 25)
ω = vec(2π .* f)

r_obs = forward(m_test, ω)

err_phi = asin(0.01) * 180 / π .* ones(length(ω))
err_appres = 0.02 * r_obs.ρₐ
err_resp = MTResponse(err_appres, err_phi)

r_obs.ρₐ .= r_obs.ρₐ .+ err_appres
r_obs.ϕ .= r_obs.ϕ .+ err_phi

respD = MTResponseDistribution(normal_dist, normal_dist)

z = 10 .^ collect(range(1, stop = 4, length = 150))
h = diff(z)
h_bounds = [[ih / 3, ih * 3] for ih in h]

modelD = MTModelDistribution(
    product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
    # product_distribution(
    #     [Uniform(ih_bounds...) for ih_bounds in h_bounds]
    # )
    vec(h),
)

n_samples = 50
mcache = mcmc_cache(modelD, respD, n_samples, NUTS())

@time mcmc_chain_nuts = stochastic_inverse(
    r_obs,
    err_resp,
    ω,
    mcache,
    trans_utils = (m = log_tf,);
    save_state = true,
)

n_samples = 30_00

mcache = mcmc_cache(modelD, respD, n_samples, MH())

log_tf2 = transform_utils([], log10, (x) -> 10^x, (x) -> inv(x * log(10)))

# @time mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache; trans_utils=(m=log_tf,), chain_type = Any);  #, resumer_from = mcmc_chain_nuts)
# mcmc_chain_nuts.value.data[end,1:150,1]
@time mcmc_chain = stochastic_inverse(
    r_obs,
    err_resp,
    ω,
    mcache;
    trans_utils = (m = log_tf,),
    init_params = mcmc_chain_nuts.value.data[end, 1:150, 1],
);  #, resumer_from = mcmc_chain_nuts)

mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache, trans_utils = (m = log_tf,));

# 887.623747 seconds (827.63 M allocations: 51.410 GiB, 12.39% gc time)
# @time stochastic_inverse(r_obs, err_resp, ω, mcache, trans_utils=(m=log_tf,))

# @profview for i in 1:100
#     m_occam.m .= 100.0;
#     inverse!(m_occam, r_obs, ω, Occam(μgrid=[1e-2, 1e6]), W=err_resp, χ2=1.0, max_iters= 2, verbose= false)
# end


model_list = get_model_list(mcmc_chain, modelD); #, trans_utils = (m =MT.lin_tf, h = MT.lin_tf));
# model_list = [MTModel(abs.(mcmc_chain.value.data[i, 1:150, 1]), vec(h)) for i in 1:10000];

plt_mods = Plots.plot()
for i = 1:(length(model_list) > 1000 ? 1000 : length(model_list))
    plot_model!(plt_mods, model_list[i], label = false)
    # prepare_plot!(plt_resps, resp_models, ω, alpha = 0.4, label = false, color = :gray); 
end
plt_mods

plt_resps = prepare_plot(r_obs, ω, alpha = 0.0, label = false);
resp_models = forward(model_list[1], ω);

for i = 1:(length(model_list) > 500 ? 500 : length(model_list))
    forward!(resp_models, model_list[i], ω)
    prepare_plot!(plt_resps, resp_models, ω, alpha = 0.4, label = false, color = :gray)
end

plt_resp_ = plot_response(plt_resps, size = (700, 700))

prepare_plot!(
    plt_resps,
    r_obs,
    ω,
    err_resp,
    markersize = 3,
    color = :orange,
    label = "true",
);
plt_resp_ = plot_response(plt_resps, size = (700, 700))

kk = MT.pre_image(modelD, mcmc_chain);#, trans_utils = (;m=MT.lin_tf, h = MT.lin_tf));
kde_plt = get_kde_image(
    kk...;
    yscale = :log10,
    xscale = :log10,
    yflip = true,
    title = "",
    clims = (0, 0.2),
)
# kde3 = get_kde_image(kk...; yscale = :log10, xscale = :log10, yflip = true, ylim = (1, 1e6))
# kde3 = get_kde_image(kk...; yscale = :identity, xscale = :log10, yflip = true) #, ylim = (1, 5e5))


mean_std_plt_lin = get_mean_std_image(
    kk...,
    yscale = :identity,
    ylim = (1, 5e5),
    xticks = 10.0 .^ collect(-1:5),
)
mean_std_plt_log = get_mean_std_image(
    kk...,
    yscale = :log10,
    ylim = (1, 5e5),
    xticks = 10.0 .^ collect(-1:5),
)
