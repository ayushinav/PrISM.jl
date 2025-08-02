using Pkg

versioninfo()
Pkg.activate("wip/.")
using Revise
using MT

using Profile
using LinearAlgebra
using BenchmarkTools
# using LinearSolve
using Turing
using Distributions
# using NonlinearSolve
using JET
using Test
using Enzyme
using DifferentiationInterface
using ForwardDiff

freqs = [55.9910   22.0020   15.4990   10.0000    5.9999    3.8751    2.5000    1.5000    0.9688    0.6250][:];
ω = 2π .* freqs

r_obs = MTResponse(
    10 .^[0.7838    0.7813    0.7378    0.7213    0.7135    0.5780    0.4802    0.3223    0.2344    0.1300][:],
    [37.3922   44.1861   48.4856   54.6212   59.8403   60.6105   61.0912   60.3664   59.8032   54.7307][:]
)

err_resp = MTResponse(
    log(10) .* r_obs.ρₐ .* [0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217    0.0217][:],
    [1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324    1.4324][:]
)

h = 20. .*ones(60);
z = [0; cumsum(h)];

respD = MTResponseDistribution(normal_dist, normal_dist)

modelD = MTModelDistribution(
        # product_distribution([truncated(Uniform(-1.0, 5.0), 0, nothing) for i in eachindex(z)]),
        product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
        vec(h),
    )


n_samples = 200
mcache = mcmc_cache(modelD, respD, n_samples, NUTS(100, 0.65)) #; adtype = AutoEnzyme(mode = set_runtime_activity(Reverse))))

# sampler_= NUTS(1000, 0.65)
# mcache = mcmc_cache(modelD, respD, n_samples, sampler_) #; adtype = AutoEnzyme(mode = set_runtime_activity(Reverse))))

import MT:*
# LogDensityFunction(mcmc_chain)
@time mcmc_chain1 = stochastic_inverse(
    r_obs,
    err_resp,
    ω,
    mcache,
    model_trans_utils = (; m = MT.lin_tf), 
);

# 100 samples => 120 sec

#=
 33.329893 seconds (19.29 M allocations: 5.475 GiB, 2.27% gc time)

julia> mcmc_chain
Chains MCMC chain (30×73×1 Array{Float64, 3}):

Iterations        = 16:1:45
Number of chains  = 1
Samples per chain = 30
Wall duration     = 33.33 seconds
Compute duration  = 33.33 seconds
=#
function fwd_wrapper(r, m, vars, response_trans_utils)
    forward!(r, m, vars, response_trans_utils)

    return typeof(r).name.wrapper(
        [getfield(r, k) for k in propertynames(r)]...
    )
end

resp_sample = MTResponse(
    [copy((getproperty(r_obs, k))) .|> Real
    for k in propertynames(r_obs)]...
)


MTResponse(
    Real[1.0, 1.0],
    Real[2f0, 1f0]
)

DiffCache(r_obs.ϕ).dual_du

# need cache ?
@time mcmc_chain2 = stochastic_inverse2(
    r_obs,
    err_resp,
    ω,
    mcache,
    model_trans_utils = (; m = MT.lin_tf),
);


@profview_allocs mcmc_chain2 = stochastic_inverse2(
    r_obs,
    err_resp,
    ω,
    mcache,
    model_trans_utils = (; m = MT.lin_tf),
);
#=
┌ Info: Found initial step size
└   ϵ = 0.05
 52.655692 seconds (30.72 M allocations: 8.719 GiB, 1.86% gc time)

  = >

┌ Info: Found initial step size
└   ϵ = 0.025
 44.474864 seconds (33.13 M allocations: 7.987 GiB, 1.61% gc time)
=#

m1 = m_rto.m
m2 = convert.(ForwardDiff.Dual, copy(m1))

resp_sample = MTResponse(
        [DiffCache(copy((getproperty(r_obs, k))), length(getproperty(r_obs, k)))
        for k in propertynames(r_obs)]...
    )

r_sample_ = get_tmp2(resp_sample, m1);


m_tt = MTModel(m1, h);
@time forward!(r_sample_, m_tt, ω)

r_sample_ = get_tmp2(resp_sample, m2);

m_tt = MTModel(m2, h);
@time forward!(r_sample_, m_tt, ω)
@btime forward!(r_sample_, m_tt, ω)
@code_warntype forward!(r_sample_, m_tt, ω)

r_sample_ = MTResponse(
        [convert.(ForwardDiff.Dual, copy((getproperty(r_obs, k))))
        for k in propertynames(r_obs)]...
    )

m_tt = MTModel(m2, h);
@btime forward!(r_sample_, m_tt, ω)
@code_warntype forward!(r_sample_, m_tt, ω)

# f = 10 .^ range(-4, stop = 1, length = 25)
# ω = vec(2π .* f)

# r_obs = forward(m_test, ω)

# err_phi = asin(0.01) * 180 / π .* ones(length(ω))
# err_appres = 0.02 * r_obs.ρₐ
# err_resp = MTResponse(err_appres, err_phi)

# r_obs.ρₐ .= r_obs.ρₐ .+ err_appres
# r_obs.ϕ .= r_obs.ϕ .+ err_phi

# respD = MTResponseDistribution(normal_dist, normal_dist)

# z = 10 .^ collect(range(1, stop = 4, length = 150))
# h = diff(z)

h = 20. .*ones(60);
z = [0; cumsum(h)];

m_rto = MTModel(fill(2.0, length(z)), h);

W = Diagonal(
    vcat([inv.(getfield(err_resp, k)) for k in propertynames(err_resp)]...)
) .^ 2;

retcode = inverse!(m_rto, r_obs, ω, 
    Occam(μgrid = [1e-2, 1e2]), W = W, χ2 = 1.0,
    max_iters = 50,
    verbose = true,
    # response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)
)
# rto_c = MT.rto_cache(m_rto, [1e-6, 1e2], Occam(), 50, 10, 1.0, [:ρₐ, :ϕ], true)

m_rto = MTModel(fill(2.0, length(z)), h);

retcode = inverse2!(m_rto, r_obs, ω, 
    Occam(μgrid = [1e-2, 1e2]), W = W, χ2 = 1.0,
    max_iters = 50,
    verbose = true,
    # response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)
)

function occam_check(m_occam, r_obs, ω)
    fill!(m_occam.m, 2.)
    inverse!(m_rto, r_obs, ω, 
        Occam(μgrid = [1e-2, 1e2]), W = W, χ2 = 1.0,
        max_iters = 50,
        verbose = false,
        # response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)
    )
end

@profview for i in 1:200 occam_check(m_rto, r_obs, ω) end

@report_call inverse!(m_rto, r_obs, ω, 
Occam(μgrid = [1e-2, 1e2]), W = W, χ2 = 1.0,
max_iters = 50,
verbose = false,
# response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)
)

@btime occam_check(m_rto, r_obs, ω)
#  268.646 ms (35459 allocations: 69.45 MiB)

function occam_check2(m_occam, r_obs, ω)
    fill!(m_occam.m, 2.)
    inverse2!(m_rto, r_obs, ω, 
        Occam(μgrid = [1e-2, 1e2]), W = W, χ2 = 1.0,
        max_iters = 50,
        verbose = false,
        # response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)
    )
end
occam_check2(m_rto, r_obs, ω)

@btime occam_check2(m_rto, r_obs, ω)
# 168.972 ms (174804 allocations: 72.66 MiB)
# 103.664 ms (57624 allocations: 82.61 MiB)

@profview for i in 1:1000 occam_check2(m_rto, r_obs, ω) end

rto_c = MT.rto_cache(m_rto, [1e-6, 1e2], Occam(), 50, 1000, 1.0, [:ρₐ, :ϕ], false)
stochastic_inverse(r_obs, err_resp, ω, rto_c)

jc = MT.jacobian_cache([:ρₐ, :ϕ], r_obs, m_rto, [:m])

@time jacobian!(jc, m_rto, ω, [:m], [:ρₐ, :ϕ])
@inferred jacobian!(jc, m_rto, ω, [:m], [:ρₐ, :ϕ])
@report_call jacobian!(jc, m_rto, ω, [:m], [:ρₐ, :ϕ])
@report_opt jacobian!(jc, m_rto, ω, [:m], [:ρₐ, :ϕ])

function fwd(m, ω)
    resp = forward(m, ω)
    return resp.ρₐ
end

# set_runtime_activity(Reverse)

autodiff(Reverse, fwd, Active, MTModel(Active(m_rto.m), Const(m_rto.h)), Const(ω))


# === PreallocationTools

using ForwardDiff

DiffCache(u)

function fwd_m!(r, m; h = h, ω = ω)
    resp = forward(
        MTModel(m, h), ω
    )

    copyto!(r, resp.ρₐ)
end
import MT:*

fj = ForwardDiff.JacobianConfig(fwd_m!, r_obs.ρₐ, m_rto.m,)

r1 = copy(r_obs.ρₐ)
m1 = copy(m_rto.m)

@time ForwardDiff.jacobian(fwd_m!, r1, m1)
@time ForwardDiff.jacobian(fwd_m!, r1, m1, fj)

fwd_m!(fj, m1)

function checkk(r1, m1)
    ForwardDiff.jacobian(fwd_m!, r1, m1)
end

function checkk2(r1, m1, fjj)
    ForwardDiff.jacobian(fwd_m!, r1, m1, fjj)
end


@btime checkk(r1, m1)
@btime checkk2(r1, m1, fj)

using ForwardDiff #, PreallocationTools
randmat = rand(5, 3)
sto = similar(randmat)
stod = DiffCache(sto)

function claytonsample!(sto, τ, α; randmat = randmat)
    sto = get_tmp(sto, τ)
    sto .= randmat
    τ == 0 && return sto

    n = size(sto, 1)
    for i in 1:n
        v = sto[i, 2]
        u = sto[i, 1]
        sto[i, 1] = (1 - u^(-τ) + u^(-τ) * v^(-(τ / (1 + τ))))^(-1 / τ) * α
        sto[i, 2] = (1 - u^(-τ) + u^(-τ) * v^(-(τ / (1 + τ))))^(-1 / τ)
    end
    return sto
end

ForwardDiff.derivative(τ -> claytonsample!(stod, τ, 0.0), 0.3)
ForwardDiff.jacobian(x -> claytonsample!(stod, x[1], x[2]), [0.3; 0.0])

Array{Float32, ndims(randmat)}(randmat)
Array{ForwardDiff.Dual, ndims(randmat)}(randmat)

Array(undef, size(randmat)...)

@model function my_model1(y, err_y, var_x)

    x ~ MvNormal(ones(2), var_x)
    y_hat = f(x)

    y ~ MvNormal(y_hat, err_y)
end



function f(x)
    return [sum(x), prod(x)]
end

xs = [2., 4.]
ys = f(xs)

err_ys = diagm([1., 1.])

s1 = my_model1(ys, err_ys, 0.5 .*ones(2));

@profview sample(s1, NUTS(), 10000);

@time sample(s1, NUTS(), 1000);
@btime sample(s1, NUTS(), 1000);
#   99.489 ms (1056799 allocations: 60.67 MiB)
function f!(y, x)

    # @show eltype(x), eltype(y)
    y[1] = sum(x) 
    y[2] = prod(x)
    return nothing
end

mutable struct D_cache{T1, T2}
    u::T1
    dual_u::T2
end

y_cache2 = D_cache(y_copy)

function get_val(x, ::Type{T}) where T
    println("inside get_val \t ", T)
    x.u
end

function get_val(x, ::Type{T}) where T <: ForwardDiff.Dual
    println("inside get_val DUAL \t ", T)
    x.dual_u
end
function D_cache(u)
    return D_cache(
        zero(u), 
        convert.(ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 2}, zero(u)))
end

y_cache = D_cache(ys)

function f2!(y, x, ::Type{T} = Float64) where T
    @show T
    y_ = get_val(y, T)
    f!(y_, x)
    nothing
end

xs2 = convert.(ForwardDiff.Dual, copy(xs))
f2!(y_cache, xs2, eltype(xs2))

@model function my_model2(y_hat, y, err_y, var_x)

    x ~ MvNormal(ones(2), var_x)
    T = eltype(x)
    y_ = get_val(y_hat, T)
    @show typeof(y_)
    @show typeof(x)
    f!(y_ , x)
 
    y ~ MvNormal(y_, err_y)
end
y_cache3 = FixedSizeDiffCache(y_copy)
y_copy = copy(ys) #.|> Real
s2 = my_model2(y_cache2, ys, err_ys, ones(2))

@time sample(s2, NUTS(), 1);
#   120.111 ms (1305306 allocations: 76.91 MiB)
fhs = ForwardDiff.JacobianConfig(f!, ys, xs)

x = convert.(ForwardDiff.Dual, [2., 3.])
yx = get_tmp(y_cache3, x)

y_c = rand(100);
y_cache3 = FixedSizeDiffCache(y_c)

@time get_tmp(y_cache3, x);

@btime f!($yx, $x)

yx2 = convert.(ForwardDiff.Dual, copy(ys))

@btime f!($yx2, $x)
