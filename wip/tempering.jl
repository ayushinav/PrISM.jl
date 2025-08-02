using Pkg

Pkg.activate("wip/.")
using MT

using Profile
using LinearAlgebra
using BenchmarkTools
using Distributions
using Turing

@model function tempering_check(m, h, ω, n, T, r_obs, W)

    m ~ Product([Uniform(-2,5) for _ in 1:n])
    T ~ Uniform(1, 10)

    resp = forward(MTModel(m,h), ω)

    r_obs ~ MultivariateNormal([resp.ρₐ..., resp.ϕ...] , T .* W) 

end

@model function tempering_check2(h, ω, n, r_obs, W)

    m ~ Product([Uniform(-2,5) for _ in 1:n])
    T ~ Uniform(1, 10)

    resp = forward(MTModel(m,h), ω)

    r_obs ~ MultivariateNormal([resp.ρₐ..., resp.ϕ...] , T .* W) 

end

m_test = MTModel(log10.([100, 1000, 100]), [1000, 1000])
ω = 2π ./ collect(10 .^ (-2:0.1:2))
robs= forward(m_test, ω);

err_robs_vec = [0.1 .* robs.ρₐ..., (180/π * asin(0.1) .* ω./ω)...]
W_ = diagm((inv.(err_robs_vec)) .^ 2)

robs_vec = [robs.ρₐ..., robs.ϕ...]

n = 50;
h = 100. .*ones(n - 1);

model_tempering_ = tempering_check2(h, ω, n, robs_vec, W_)
sample_ = Turing.sample(model_tempering_, MH(), 50_000)

mt_chain = Turing.Chains(
        (sample_.value.data[:,1:n,:]), [Symbol("ρ[$i]") for i in 1:n]);

mDist =  MTModelDistribution(Product([Uniform(-1, 5) for i = 1:n]), h)
mt_model_list = get_model_list(mt_chain, mDist);
# mt_model_list[111] = copy(mt_model_list[1]);
n_samples = length(mt_model_list)


using Plots
plt_mods1 = Plots.plot()
for i in 1:(length(mt_model_list) > 500 ? 500 : length(mt_model_list))
    plot_model!(plt_mods1, mt_model_list[i], label = false, color = :green, alpha = 0.4);
    # prepare_plot!(plt_resps, resp_models, ω, alpha = 0.4, label = false, color = :gray); 
end

plt_mods1

plot(plt_mods1, yscale = :identity, ylim = (0, 3000), 
        label = false, xlim = (0.01, 1e6), xticks = 10 .^ (-2:1.:6), linewidth = 3, color = :black, alpha = 1.)

r2 = copy(robs);
rms_vals = zeros(length(mt_model_list))
for i in 1:length(mt_model_list)
    forward!(r2, mt_model_list[i], ω);
    rms_vals[i] = MT.χ²([r2.ρₐ..., r2.ϕ...], [robs.ρₐ..., robs.ϕ...], W=W_)
end

# histogram(rms_vals, ylabel = "pdf", xlabel = "rms", label = false, title = "rms histogram", normalize = :pdf)
# # savefig("rms_vals.png")

# histogram(sample_.value.data[:,end,1], normalize =:pdf, label = false, xlabel = "μ", ylabel = "pdf")
# savefig("μ_vals.png")

plt_resps = prepare_plot(robs, ω, alpha = 0., label = false)
resp_models = forward(mt_model_list[1], ω);

for i in 1:(length(mt_model_list) > 500 ? 500 : length(mt_model_list))
    forward!(resp_models, mt_model_list[i], ω);
    prepare_plot!(plt_resps, resp_models, ω, alpha = 0.1, label = false, color = :gray); 
end

prepare_plot!(plt_resps, robs, ω, markersize = 4, 
    msc = :orange, color =:red, label = false, falpha = 0.);
plt_resp_ = plot_response(plt_resps, size = (1000, 700), margin = 0Plots.mm)

