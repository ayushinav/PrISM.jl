using Pkg
Pkg.activate(".")

using MT
using LinearAlgebra
using BenchmarkTools
using NonlinearSolve, Optimization, OptimizationOptimJL

freqs = [55.9910   22.0020   15.4990   10.0000    5.9999    3.8751    2.5000    1.5000    0.9688    0.6250][:];
ω = 2π .* freqs
r_obs = MTResponse(
    10. .^[0.7838    0.7813    0.7378    0.7213    0.7135    0.5780    0.4802    0.3223    0.2344    0.1300][:],
    [37.3922   44.1861   48.4856   54.6212   59.8403   60.6105   61.0912   60.3664   59.8032   54.7307][:]
)

err_resp = MTResponse(
    r_obs.ρₐ .* log(10) .* [0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217][:],
    [1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324][:]
)

r_obs2 = MTResponse(
    [0.7838    0.7813    0.7378    0.7213    0.7135    0.5780    0.4802    0.3223    0.2344    0.1300][:],
    [37.3922   44.1861   48.4856   54.6212   59.8403   60.6105   61.0912   60.3664   59.8032   54.7307][:]
)

err_resp2 = MTResponse(
    [0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217][:],
    [1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324][:]
)

# err_resp.ρₐ = log(10) .* err_resp.ρₐ .* r_obs.ρₐ;

W = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2;

m_occam = MTModel(2. *ones(61), 20. *ones(60))
inverse!(m_occam, r_obs, ω, Occam(μgrid = [1e-2, 1e2]), W = W, max_iters = 50, χ2 = 1.) #, response_trans_utils = (ρₐ = log_tf,ϕ = MT.lin_tf))

m_gn = MTModel(2. *ones(61) .+ (0.001 .* randn(61)), 20. *ones(60))
inverse!(m_gn, r_obs, ω, NonlinearAlg(alg = LevenbergMarquardt, μ = 1e-2), W = W, max_iters = 100, χ2 = 1.) #, response_trans_utils = (ρₐ = log_tf,ϕ = MT.lin_tf))

m_bfgs = MTModel(2. *ones(61) .+ (0.01 .* randn(61)), 20. *ones(60))
inverse!(m_bfgs, r_obs, ω, MT.OptAlg(alg = LBFGS, μ = 1e-1), W = W, max_iters = 100, χ2 = 1.) #, response_trans_utils = (ρₐ = log_tf,ϕ = MT.lin_tf))

plt_mod = plot_model(m_occam, yscale = :identity, ylim = (0, 1200), yticks = collect(0:100:1200), 
        label = "occam", xlim = (0.1, 100), xticks = 10 .^ (-1:0.5:2))

plot_model!(plt_mod, m_gn, yscale = :identity, ylim = (0, 1200), yticks = collect(0:100:1200), 
label = "gn", xlim = (0.001, 10000), xticks = 10 .^ (-1:0.5:2))

plot_model!(plt_mod, m_bfgs, yscale = :identity, ylim = (0, 1200), yticks = collect(0:100:1200), 
label = "bfgs", xlim = (0.001, 10000), xticks = 10 .^ (-1:0.5:2))

plt_resp = prepare_plot(r_obs, ω, err_resp, plot_type = :scatter, color = "green", label = "true")
prepare_plot!(plt_resp, forward(m_occam, ω), ω, plot_type = :plot, color =  "blue", label = "occam")
prepare_plot!(plt_resp, forward(m_gn, ω), ω, plot_type = :plot, color = "red", label = "gn")
prepare_plot!(plt_resp, forward(m_bfgs, ω), ω, plot_type = :plot, color = "black", label = "bfgs")
plt_resp_ = plot_response(plt_resp)

norm(diff((m_occam.m)))
norm(diff((m_gn.m)))

m_occam = MTModel(100. *ones(61), 20. *ones(60))

@time inverse!(m_occam, r_obs, ω, Occam(μgrid = [0.1, 2.]), max_iters = 30, verbose = true)

m_occam = MTModel(100. *ones(60), 20. *ones(59))

@benchmark inverse!(m_occam, r_obs, ω, Occam(μgrid = [1e-2, 1e2]), max_iters = 30, verbose = true)

function fn1(m_occam, r_obs, ω)
    fill!(m_occam.m, 100.)
    inverse!(m_occam, r_obs, ω, Occam(μgrid = [1e-2, 1e2]), max_iters = 30, verbose = true)
end

function fn3(m_occam)
    fill!(m_occam.m, 100.)
end

fn3(m_occam)
@btime fn3(m_occam)

# 37.848 ns (0 allocations: 0 bytes)

fn1(m_occam, r_obs, ω)
@btime fn1(m_occam, r_obs, ω)

# 155.191 ms


nlc = NonlinearAlg(alg = LevenbergMarquardt, μ = 1e-2)
typeof(nlc)

nlc = NonlinearAlg(alg = LevenbergMarquardt, μ = 1e-2)
typeof(nlc)
