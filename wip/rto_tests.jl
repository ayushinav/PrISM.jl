using Pkg
# Pkg.activate("pkgs/MT.jl/wip/")
Pkg.activate(".")

using MT
using Distributions
using LinearAlgebra
using Turing

m_test = MTModel([100.0, 10.0, 100.0], [1000.0, 1000.0])

T = 10 .^ collect(range(-2, 4, length = 41));
ω = 2π .* T;

robs = forward(m_test, ω)

err_robs = MTResponse(robs.ρₐ .* 0.1, fill(asin(0.05) * 180 / π, length(T)))

plt_model = plot_model(m_test)

z_occam = 10.0 .^ collect(range(-1, 4, 100));
h_occam = diff(z_occam);
m_occam = MTModel(fill(100.0, length(z_occam)), h_occam);

W = diagm(inv.([err_robs.ρₐ..., err_robs.ϕ...])) .^ 2;

inverse!(
    m_occam,
    robs,
    ω,
    Occam(μgrid = [1e-2, 1e2]),
    W = W,
    χ2 = 1.0,
    max_iters = 50,
    verbose = true,
)
plot_model!(plt_model, m_occam)

m_rto = MTModel(fill(100.0, length(z_occam)), h_occam);

rto_c = MT.rto_cache(m_rto, [1e-2, 1e2], 50, 50, 1.0, [:ρₐ, :ϕ], false)

rto_chains = stochastic_inverse(robs, err_robs, ω, rto_c)

# rto_chains2 = Turing.Chains(rto_chains.value.data[1:90,:,:]
#     , [Symbol("m[$i]") for i in 1:(101)])

m_dist = MTModelDistribution(Product([Uniform(-1, 5) for i = 1:length(z_occam)]), h_occam)

mt_chain = Turing.Chains(
    log10.(rto_chains.value.data[:, 1:100, :]),
    [Symbol("ρ[$i]") for i = 1:100],
);

mt_model_list = get_model_list(mt_chain, m_dist);
n_samples = rto_c.n_samples


plt_mods = Plots.plot()
for i = 1:(length(mt_model_list) > 500 ? 500 : length(mt_model_list))
    plot_model!(plt_mods, mt_model_list[i], label = false)
    # prepare_plot!(plt_resps, resp_models, ω, alpha = 0.4, label = false, color = :gray); 
end

plt_mods

plot!(plt_mods, xlims = (1, 1e5))

plt_resps = prepare_plot(robs, ω, alpha = 0.0, label = false);
resp_models = forward(mt_model_list[1], ω);

for i = 1:(length(mt_model_list) > 2000 ? 2000 : length(mt_model_list))
    forward!(resp_models, mt_model_list[i], ω)
    prepare_plot!(plt_resps, resp_models, ω, alpha = 0.4, label = false, color = :gray)
end

prepare_plot!(
    plt_resps,
    robs,
    ω,
    err_robs,
    markersize = 3,
    msc = :orange,
    color = :orange,
    label = "d_obs",
);
plt_resp_ = plot_response(plt_resps, size = (1000, 700), margin = 0Plots.mm)


pre_img = pre_image(m_dist, mt_chain)
kde_img = get_kde_image(
    pre_img...,
    # true,
    xscale = :log10,
    yscale = :log10,
    yflip = true,
    clims = (0.0, 0.02),
)

trans_utils = (m = log_tf, h = MT.lin_tf);
pre_img = pre_image(m_dist, mt_chain, trans_utils = trans_utils);

pre_img = pre_image(m_dist, mt_chain) #, grid = grid_);
kde_vals, kde_img = get_kde_image(
    pre_img...,
    true,
    xscale = :log10,
    yscale = :log10,
    yflip = true,
    clims = (0.0, 0.3),
)
# kde_img = get_kde_image(pre_img..., xscale = :log10, yscale = :log10, yflip = true, clims = (0., 0.2))

plt_dist = heatmap(
    pre_img[2][:m],
    pre_img[1][:h],
    log10.(kde_vals .+ 1e-9),
    yflip = true,
    clims = (-5, 0),
    yscale = :log10,
    yticks = 10 .^ (0:6),
)


mean_std_plt_lin = get_mean_std_image(
    pre_img...,
    trans_utils = trans_utils,
    yscale = :identity,
    ylim = (1, 1e4),
    xticks = 10.0 .^ collect(-1:5),
    ylabel = "depth (km)",
)

lay_out = Plots.@layout [Plots.grid(2, 1) b{0.6w}];

plot(plt_resp_, plt_mods, plt_dist, layout = lay_out, size = (1600, 1000))


L = MT.∂(length(m_rto.m));

inv(L'L)

L'L

# μ = 
ξ_ = rand(MvNormal(fill(0.0, length(z_occam)), fill(1.0, length(z_occam))))
μ_ = 10.0^rand(Uniform(-2, 6))
kk_ = inv(sqrt(μ_)) .* inv(L) * ξ_

broadcast!((x) -> (10.0^MT.sigmoid_tf.tf(x)), kk_, kk_)
