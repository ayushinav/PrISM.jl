#=
 TODO
reg term in occam and inverse
NonlinearSolve.jl
=#

using Pkg
Pkg.activate(".")

using MT
using LinearAlgebra
using BenchmarkTools
using Distributions
using Turing
using NonlinearSolve, Optimization, OptimizationOptimJL
using JLD2

sts = ["2_yx", "3_yx", "6_yx", "4_yx"]
# for st in sts 
path_to_kivu = "../../kivu/"
st = "2_yx"
kk = jldopen(path_to_kivu * "filter_data_err_floor/" * st * ".jld2")
r_obs = kk["r_obs4"];
err_resp = kk["resp_err4"];
ω = kk["ω"];


idxs = sortperm(ω);
r_obs.ρₐ .= r_obs.ρₐ[idxs];
r_obs.ϕ .= r_obs.ϕ[idxs];
err_resp.ρₐ .= err_resp.ρₐ[idxs];
err_resp.ϕ .= err_resp.ϕ[idxs];
ω .= ω[idxs];

W = Diagonal(vcat([inv.(getfield(err_resp, k)) for k in propertynames(err_resp)]...)) .^ 2;

z = 10 .^ collect(range(1, stop = 6, length = 100));
h = diff(z);

model_occam = MTModel(fill(100.0, length(z)), vec(h));
r2 = copy(r_obs)
function cf1(m, vars)
    # vars[:model].m .= m
    # @show propertynames(vars)
    resp = forward(MTModel(m, vars[:h]), vars[:ω])
    resp = forward(MTModel(m, vars[:h]), vars[:ω])
    L1 =
        ([resp.ρₐ..., resp.ϕ...] .- [vars[:r_obs].ρₐ..., vars[:r_obs].ϕ...])' *
        vars[:W] *
        ([resp.ρₐ..., resp.ϕ...] .- [vars[:r_obs].ρₐ..., vars[:r_obs].ϕ...])
    L2 = norm(MT.∂(length(m)) * m)
    # @show L1, L2
    return [L1..., L2]
end


prob = NonlinearLeastSquaresProblem{false}(
    cf1,
    m0,
    (; model = m_test, r_obs = r_obs, W = sqrt.(diag(W)), ω = ω, μ = μ, resp = r2, h = h),
)
sol = solve(prob, GaussNewton(), maxiters = 100)

function cost_function3!(L, m, vars)
    # @show propertynames(vars)
    forward!(vars[:resp], MTModel(m, vars[:h]), vars[:ω])
    # L1 = ([vars[:resp].ρₐ..., vars[:resp].ϕ...] .- [vars[:r_obs].ρₐ..., vars[:r_obs].ϕ...]) .* vars[:W]
    # L2 = MT.∂(length(m)) * m
    # copyto!(L, [L1..., L2...])
    # L[1:end-length(m)] .= [(vars[:resp].ρₐ .- vars[:r_obs].ρₐ)..., (vars[:resp].ϕ .- vars[:r_obs].ϕ)...] .* vars[:W]
    # L[1:end-length(m)] .= ([vars[:resp].ρₐ..., vars[:resp].ϕ...] .- [vars[:r_obs].ρₐ..., vars[:r_obs].ϕ...]) .* vars[:W]
    # L[end - length(m) + 1: end] .= sqrt(first(vars[:μ])) .* MT.∂(length(m)) * m
    # @show L1, L2
    return nothing
end

function cost_function3!(L, m, vars)
    copyto!(getfield(vars[:model], :m), m)
    forward!(vars[:resp], vars[:model], vars[:ω])
    # L1 = ([vars[:resp].ρₐ..., vars[:resp].ϕ...] .- [vars[:r_obs].ρₐ..., vars[:r_obs].ϕ...]) .* vars[:W]
    # L2 = MT.∂(length(m)) * m
    # copyto!(L, [L1..., L2...])
    @show size(L)
    copyto!(
        L[1:end-length(m)],
        [(vars[:resp].ρₐ .- vars[:r_obs].ρₐ)..., (vars[:resp].ϕ .- vars[:r_obs].ϕ)...] .*
        vars[:W],
    )
    # L[1:end-length(m)] .= ([vars[:resp].ρₐ..., vars[:resp].ϕ...] .- [vars[:r_obs].ρₐ..., vars[:r_obs].ϕ...]) .* vars[:W]
    L[end-length(m)+1:end] .= sqrt(first(vars[:μ])) .* MT.∂(length(m)) * m
    # @show L1, L2
    return nothing
end

m_test = copy(model_occam);
# resp = copy(r_obs)

L = zeros(2length(ω) + length(m0));
@time cost_function3!(
    L,
    m0,
    (; model = m_test, r_obs = r_obs, W = sqrt.(diag(W)), ω = ω, μ = μ, resp = r2, h = h),
)

prob = NonlinearLeastSquaresProblem{true}(
    cost_function3!,
    m0,
    (; model = m_test, r_obs = r_obs, W = sqrt.(diag(W)), ω = ω, μ = μ, resp = r2, h = h),
)
sol = solve(prob, GaussNewton(), maxiters = 100)

function sigmoid2(m, bounds)
    σ(x) = 1 / (1 + exp(-x / (bounds[2] - bounds[1])))
    # return 10^(σ(m)*log10(bounds[2]/bounds[1])+ log10(bounds[1]))
    return σ(m) * (bounds[2] - bounds[1]) + bounds[1]
end
sg(x) = sigmoid2(x, [-2, 6])
resp_ = copy(r_obs);
@time forward!(resp_, m_test, ω)

μ = [1e-1]
function cf5(m, p)
    m2 = sg.(m)
    resp = forward(MTModel(m2, h), ω)
    L1 =
        ([resp.ρₐ..., resp.ϕ...] .- [r_obs.ρₐ..., r_obs.ϕ...])' *
        W *
        ([resp.ρₐ..., resp.ϕ...] .- [r_obs.ρₐ..., r_obs.ϕ...])
    L2 = norm(MT.∂(length(m)) * sg.(m))
    @show L1, L2
    # @show sqrt(L1/50)
    return (L1 .+ μ * L2)
end

m0 = zero(model_occam.m) .+ 2.0; # .+ randn(length(model_occam.m)).*10;
@time cost_function3(m0, (; r_obs = r_obs, W = W, ω = ω, μ = μ, resp = r2, h = h))


@time cf5(m0, []) #, (;r_obs = r_obs, W = W, ω = ω, μ = [1.], resp = r2, h = h))

# prob = NonlinearProblem(cost_function3, m0, (;r_obs = r_obs, W = W, ω = ω, μ = 1., resp = r2, h = h))
# prob = NonlinearProblem(cost_function5, m0, [])
# sol = solve(prob, GaussNewton())
# mutable struct params

prob = SciMLBase.NonlinearLeastSquaresProblem(
    cf5,
    m0,
    (; r_obs = r_obs, W = W, ω = ω, μ = 1.0, resp = r2, h = h),
)
nlcache = init(prob, GaussNewton(); store_trace = Val(true))

# @time step!(nlcache)
@time for i = 1:30
    # check condition
    step!(nlcache)
end

@time sol = solve!(nlcache; maxiters = 30) #, store_trace = Val(true))

@time sol = solve(prob, NewtonRaphson(); maxiters = 30, store_trace = Val(true))

@time cf5(sol.u, [])


optfn = OptimizationFunction(cf5, Optimization.AutoForwardDiff())
prob = OptimizationProblem(optfn, m0) #, lcons = m0./m0 .*1e-2, ucons = m0./m0 .*1e6) #, (;r_obs = r_obs, W = W, ω = ω, μ = 1., resp = r2, h = h))
@time sol = solve(prob, BFGS(), maxiters = 100)

# cost_function5(sol.u, []) #, (;r_obs = r_obs, W = W, ω = ω, μ = [1.], resp = r2, h = h))


m2 = sg.(sol.u)
norm(diff(sol.u)) * first(μ)
r2 = forward(MTModel(m2, h), ω)
plt = prepare_plot(r2, ω, label = "pred")
prepare_plot!(plt, r_obs, ω, label = false)

# sol.u[sol.u .<= 0.] .= 0.1

plot_response(plt)

plt_mod = plot_model(MTModel(m2, h))

inverse!(
    model_occam,
    r_obs,
    ω,
    Occam(μgrid = [1e-2, 1e6]),
    W = W,
    χ2 = 1.0,
    max_iters = 50,
    verbose = false,
)
plot_model!(plt_mod, model_occam)


# Define the problem to solve
using Optimization, ForwardDiff, Zygote

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
_p = [1.0, 100.0]
rosenbrock(x0, _p)

ff = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff())
l1 = rosenbrock(x0, _p)
prob = OptimizationProblem(ff, x0, _p)
