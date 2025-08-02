using Pkg
Pkg.activate(".")

using CairoMakie
using MT
using Turing
using Distributions

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


n_samples = 10_000
mcache = mcmc_cache(modelD, respD, n_samples, Prior()) #; adtype = AutoEnzyme(mode = set_runtime_activity(Reverse))))

@time mcmc_chain1 = stochastic_inverse(
    r_obs,
    err_resp,
    ω,
    mcache,
    model_trans_utils = (; m = MT.lin_tf), 
);


set_theme!(theme_light())

gaussian_kernel(u, σ² = 2) = inv(sqrt(2π)) * exp(-u^2 / σ²)

function get_kde(data, xgrid; K = gaussian_kernel)
    σ = std(data)
    h = 1.06 * σ
    px = zeros(size(xgrid))
    for (i, x) in enumerate(xgrid)
        s = 0
        for idata in data
            s = s + K((x - idata) / h)
        end
        px[i] = inv(length(data) * h) * s
    end
    return px
end

## fixed discretization


function get_kde_image!(fig, chain::C, mDist::mdist;
    cb_kwargs = (;),
    hm_kwargs = (;),
    K = gaussian_kernel,
    half_space_thickness = nothing,
    return_kde_mat = false,
    trans_utils = (m = lin_tf, h = lin_tf),
    grid=(m=collect(-1:0.1:5), z= cumsum(mDist.h))) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:AbstractArray}
    }

    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...) 

    h_length = length(mDist.h)
    kde_img = zeros(length(grid.m), h_length)  # nₘ x nₕ

    for i in 1:h_length
        kde_img[:, i] .= get_kde(pred[:, i], grid.m, K = K)
        norm_factor = sum(kde_img[:, i])
        kde_img[:, i] .= kde_img[:, i] ./ norm_factor
    end

    isnothing(half_space_thickness) && (half_space_thickness = sum(mDist.h) * 1.25)
    zs = [0, cumsum(mDist.h)..., half_space_thickness]

    ax = Axis(fig[1,1])
    cb = fig[1,2]

    ms = broadcast(trans_utils.m.tf, grid.m)
    hm = heatmap!(ax, ms, zs, kde_img; hm_kwargs...)
    ax.yreversed = true

    Colorbar(cb, hm, cb_kwargs...)

    if return_kde_mat
        return kde_img
    else
        return nothing
    end

end

function get_kde_image!(fig, chain::C, mDist::mdist;
    cb_kwargs = (;),
    hm_kwargs = (;),
    K = gaussian_kernel,
    half_space_thickness = nothing,
    return_kde_mat = false,
    trans_utils = (m = lin_tf, h = lin_tf),
    grid=(m=collect(-1:0.1:5), z= cumsum(mean(mDist.h)))) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:Distribution}
    }

    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...) 

    m_length = length(rand(mDist.m))
    h_length = length(rand(mDist.h))

    broadcast!(trans_utils.m.tf,
        view(pred, :, 1:m_length), view(pred, :, 1:m_length))
    broadcast!(trans_utils.h.tf,
        view(pred, :, (m_length + 1):(m_length + h_length)),
        view(pred, :, (m_length + 1):(m_length + h_length))
    )

    isnothing(half_space_thickness) && (half_space_thickness = sum(mean(mDist.h)) * 1.25)
    zs = [0, grid.z..., grid.z[end] .+ range(0., half_space_thickness, length = 10)...]

    m2 = zeros(eltype(pred), size(pred,1), length(zs))
    for i in axes(pred, 1)
        m2[i, :] .= MT.get_ρ_at_z(pred[i, :], zs)
    end

    z_length = length(zs)
    kde_img = zeros(length(grid.m), z_length)  # nₘ x nₕ

    for i in 1:z_length
        kde_img[:, i] .= get_kde(m2[:, i], grid.m, K = K)
        norm_factor = sum(kde_img[:, i])
        kde_img[:, i] .= kde_img[:, i] ./ norm_factor
    end
    
    ax = Axis(fig[1,1])
    cb = fig[1,2]

    ms = broadcast(trans_utils.m.tf, grid.m)
    hm = heatmap!(ax, ms, zs, kde_img; hm_kwargs...)
    ax.yreversed = true

    Colorbar(fig[1,2], hm, cb_kwargs...)

    if return_kde_mat
        return kde_img
    else
        return nothing
    end

end

function get_kde_image2(args...; return_kde_mat = false, kwargs...)
    fig = Figure()
    kde_img = get_kde_image!(fig, args...; kwargs...)

    if return_kde_mat
        return fig, kde_img
    else
        return fig
    end

end

function get_mean_std_image!(ax, chain::C, mDist::mdist;
    confidence_interval = 0.95,
    half_space_thickness = nothing,
    plot_kwargs = nothing,
    trans_utils = (m = lin_tf, h = lin_tf),
    ) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:AbstractArray}
    }

    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...) 

    μ_m = mean(pred, dims = 1)[:]
    σ_m = mean(pred, dims = 1)[:]

    z_score = quantile(Normal(0., 1.), 1 - (1 - confidence_interval)/2) 

    μ₊_m = @. μ_m + z_score * σ_m
    μ₋_m = @. μ_m - z_score * σ_m

    isnothing(half_space_thickness) && (half_space_thickness = sum(mDist.h) * 1.25)
    mean_kwargs = (;)
    std_plus_kwargs = (;)
    std_minus_kwargs = (;)
    if isnothing(plot_kwargs)
        mean_kwargs = (mean_kwargs...,
        label = "mean", color = :blue)

        std_plus_kwargs = (std_plus_kwargs...,
        label = "upper bound", color = :green)

        std_minus_kwargs = (std_minus_kwargs...,
        label = "lower bound", color = :green)
    else
        mean_kwargs = (; plot_kwargs.mean_kwargs...)
        std_plus_kwargs = (; plot_kwargs.std_plus_kwargs...)
        std_minus_kwargs = (; plot_kwargs.std_minus_kwargs...)
    end

    m_type = MT.sample(mDist)

    plot_model!(ax, m_type(μ_m, mDist.h); mean_kwargs...)
    plot_model!(ax, m_type(μ₊_m, mDist.h); std_plus_kwargs...)
    plot_model!(ax, m_type(μ₋_m, mDist.h); std_minus_kwargs...)

    nothing
end

function get_mean_std_image!(ax, chain::C, mDist::mdist;
    confidence_interval = 0.95,
    half_space_thickness = nothing,
    plot_kwargs = nothing,
    trans_utils = (m = lin_tf, h = lin_tf),
    z_points = cumsum(mean(mDist.h))
    ) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:Distribution}
    }

    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...) 

    m_length = length(rand(mDist.m))
    h_length = length(rand(mDist.h))

    broadcast!(getproperty(trans_utils[:m], :tf),
        view(pred, :, 1:m_length), view(pred, :, 1:m_length))
    broadcast!(getproperty(trans_utils[:h], :tf),
        view(pred, :, (m_length + 1):(m_length + h_length)),
        view(pred, :, (m_length + 1):(m_length + h_length))
    )

    isnothing(half_space_thickness) && (half_space_thickness = sum(mean(mDist.h)) * 1.25)
    zs = [0, z_points..., z_points[end] .+ range(0., half_space_thickness, length = 10)...]

    m2 = zeros(eltype(pred), size(pred,1), length(zs))
    for i in axes(pred, 1)
        m2[i, :] .= MT.get_ρ_at_z(pred[i, :], zs)
    end

    μ_m = mean(m2, dims = 1)[:]
    σ_m = mean(m2, dims = 1)[:]

    z_score = quantile(Normal(0., 1.), 1 - (1 - confidence_interval)/2) 

    μ₊_m = @. μ_m + z_score * σ_m
    μ₋_m = @. μ_m - z_score * σ_m

    isnothing(half_space_thickness) && (half_space_thickness = sum(mDist.h) * 1.25)
    mean_kwargs = (;)
    std_plus_kwargs = (;)
    std_minus_kwargs = (;)
    if isnothing(plot_kwargs)
        mean_kwargs = (mean_kwargs...,
        label = "mean", color = :blue)

        std_plus_kwargs = (std_plus_kwargs...,
        label = "upper bound", color = :green)

        std_minus_kwargs = (std_minus_kwargs...,
        label = "lower bound", color = :green)
    else
        mean_kwargs = (; plot_kwargs.mean_kwargs...)
        std_plus_kwargs = (; plot_kwargs.std_plus_kwargs...)
        std_minus_kwargs = (; plot_kwargs.std_minus_kwargs...)
    end

    m_type = MT.sample(mDist)

    # plot_model!(ax, m_type(μ_m, mDist.h); mean_kwargs...)
    lines!(ax, μ_m, zs; mean_kwargs...)
    lines!(ax, μ₊_m, zs; mean_kwargs...)
    lines!(ax, μ₋_m, zs; mean_kwargs...)
    # plot_model!(ax, m_type(μ₊_m, mDist.h); std_plus_kwargs...)
    # plot_model!(ax, m_type(μ₋_m, mDist.h); std_minus_kwargs...)

    nothing
end

function get_mean_std_image(args...; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    get_mean_std_image!(ax, args...; kwargs...)
    fig
end

f = Figure()
hm_kws = (; colormap = :jet, colorrange = (0.001, 0.03))
cb_kws = (;)
get_kde_image!(f, mcmc_chain1, modelD; trans_utils = (m = MT.pow_tf, h = lin_tf), hm_kwargs = hm_kws, cb_kwargs = cb_kws)
ax, cb = f.content
ax.xscale = log10
f

ff = get_kde_image(mcmc_chain1, modelD, trans_utils = (m = MT.pow_tf, h = lin_tf),
hm_kwargs = hm_kws, cb_kwargs = cb_kws)
ax, cb = ff.content
ax.xscale = log10

f = Figure()
hm_kws = (; colormap = :jet) #, colorrange = (0.001, 0.03))
cb_kws = (;)
get_kde_image!(f, mcmc_chain2, modelD2; trans_utils = (m = MT.pow_tf, h = lin_tf), hm_kwargs = hm_kws, cb_kwargs = cb_kws)
ax, cb = f.content
ax.xscale = log10
f

ff = get_kde_image(mcmc_chain2, modelD2, trans_utils = (m = MT.pow_tf, h = lin_tf),
hm_kwargs = hm_kws, cb_kwargs = cb_kws)
ax, cb = ff.content
ax.xscale = log10
ff

f = get_mean_std_image(mcmc_chain1, modelD)
ax = f.content[1]
ax.xscale = log10
f

f = get_mean_std_image(mcmc_chain2, modelD2, trans_utils = (m = MT.pow_tf, h = lin_tf))
ax = f.content[1]
ax.xscale = log10
f

ax, cb = f.content

ax.xscale = log10

get_mean_std_image!(ax, mcmc_chain1, modelD, confidence_interval = 0.5)

f2 = f[2,1]
ax2 = Axis(f2)

get_mean_std_image!(ax2, mcmc_chain1, modelD, confidence_interval = 0.5)
ax2.xscale = log10
xlims!(ax2, (0.1, 1e5))
xlims!(ax, (0.1, 1e5))

f, ax, cb = get_kde_image(mcmc_chain1, modelD, trans_utils = (m = MT.pow_tf, h = lin_tf))

ax = f.content[1]
cb = f.content[2]

ax.xscale = log10


xs = range(0, 2π, length=100)
ys = range(0, 2π, length=100)
zs1 = [sin(x*y) for x in xs, y in ys]
zs2 = [2sin(x*y) for x in xs, y in ys]

joint_limits = (-2, 2)  # here we pick the limits manually for simplicity instead of computing them

fig = Figure()
ax1 = Axis(fig[1,1])

# fig, ax1, hm1 = heatmap(xs, ys, zs1,  colorrange = joint_limits)
hm1 = heatmap!(ax1, xs, ys, zs1,  colorrange = joint_limits)
ax2, hm2 = heatmap(fig[1, end+1], xs, ys, zs2, colorrange = joint_limits)

Colorbar(fig[:, end+1], hm1)                     # These three
Colorbar(fig[:, end+1], hm2)                     # colorbars are
Colorbar(fig[:, end+1], colorrange = (-3,3))  # equivalent

fig


get_mean_std_image!(ax, mcmc_chain1, modelD, confidence_interval = 0.5)

cb.axis.attributes[:limits][] = (1, 2)

f

Colorbar(cb, colorrange = (1,2))

ax = f.content[1]
ax.xscale = log10
cb.axis.attributes[:colorrange][] = (1e-2, 5.)

cb = f.content[2]
cb.axis.attributes[:limits][] = (1e-2, 1.)
cb.axis.attributes[:scale][] = log10
f

function get_mean_std_image!(ax, chain::C, mDist::mdist;
    confidence_interval = 0.95,
    half_space_thickness = nothing,
    plot_kwargs = nothing,
    ) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:AbstractArray}
    }

    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...) 

    μ_m = mean(pred, dims = 1)[:]
    σ_m = mean(pred, dims = 1)[:]


    z_score = quantile(Normal(0., 1.), 1 - (1 - confidence_interval)/2) 

    μ₊_m = @. μ_m + z_score * σ_m
    μ₋_m = @. μ_m - z_score * σ_m

    isnothing(half_space_thickness) && (half_space_thickness = sum(mDist.h) * 1.25)
    mean_kwargs = (;)
    std_plus_kwargs = (;)
    std_minus_kwargs = (;)
    if isnothing(plot_kwargs)
        mean_kwargs = (mean_kwargs...,
        label = "mean", color = :blue)

        std_plus_kwargs = (std_plus_kwargs...,
        label = "upper bound", color = :green)

        std_minus_kwargs = (std_minus_kwargs...,
        label = "lower bound", color = :green)
    else
        mean_kwargs = (; plot_kwargs.mean_kwargs...)
        std_plus_kwargs = (; plot_kwargs.std_plus_kwargs...)
        std_minus_kwargs = (; plot_kwargs.std_minus_kwargs...)
    end

    m_type = MT.sample(mDist)

    plot_model!(ax, m_type(μ_m, mDist.h); mean_kwargs...)
    plot_model!(ax, m_type(μ₊_m, mDist.h); std_plus_kwargs...)
    plot_model!(ax, m_type(μ₋_m, mDist.h); std_minus_kwargs...)

    nothing
end

function get_mean_std_image(args...; kwargs...)

    fig = Figure()
    ax = Axis(fig[1,1])

    get_mean_std_image!(ax, args...; kwargs...)

    return fig, ax
end

F = Figure()
f1 = F[1,1] = GridLayout()
f2 = F[1,2] = GridLayout()

ax = Axis(f1[1,1])
get_mean_std_image!(f1, mcmc_chain1, modelD, confidence_interval = 0.9)

f2 = get_mean_std_image(mcmc_chain1, modelD, confidence_interval = 0.9)
ax2 = f2.content[1]
ax2.xscale = log10

get_mean_std_image!(f, mcmc_chain1, modelD)
ax.xscale = log10
f

ρ= log10.([500., 100., 400., 1000.]);
h= [100., 900., 9000.];
m= MTModel(ρ, h)

plot_model(m)
ax_ = f1.content[1]
stairs!(ax, rand(5), rand(5))


modelD2 = MTModelDistribution(
        # product_distribution([truncated(Uniform(-1.0, 5.0), 0, nothing) for i in eachindex(z)]),
        product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
        product_distribution([Uniform(10.0, 30.0) for i in eachindex(h)]),
)


n_samples = 10_000
mcache = mcmc_cache(modelD2, respD, n_samples, Prior()) #; adtype = AutoEnzyme(mode = set_runtime_activity(Reverse))))

@time mcmc_chain2 = stochastic_inverse(
    r_obs,
    err_resp,
    ω,
    mcache,
    model_trans_utils = (; m = MT.lin_tf), 
);

## fixed discretization
function get_kde_image(chain::C, mDist::mdist;
    K = gaussian_kernel,
    half_space_thickness = nothing,
    return_kde_mat = false,
    trans_utils = (m = lin_tf, h = lin_tf),
    grid=(m=collect(-1:0.1:5), z= cumsum(mean(mDist.h)))) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:Distribution}
    }

    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...) 

    m_length = length(rand(mDist.m))
    h_length = length(rand(mDist.h))

    broadcast!(getproperty(trans_utils[:m], :tf),
        view(pred, :, 1:m_length), view(pred, :, 1:m_length))
    broadcast!(getproperty(trans_utils[:h], :tf),
        view(pred, :, (m_length + 1):(m_length + h_length)),
        view(pred, :, (m_length + 1):(m_length + h_length))
    )

    isnothing(half_space_thickness) && (half_space_thickness = sum(mean(mDist.h)) * 1.25)
    zs = [0, grid.z..., grid.z[end] .+ range(0., half_space_thickness, length = 10)...]

    m2 = zeros(eltype(pred), size(pred,1), length(zs))
    for i in axes(pred, 1)
        m2[i, :] .= MT.get_ρ_at_z(pred[i, :], zs)
    end

    z_length = length(zs)
    kde_img = zeros(length(grid.m), z_length)  # nₘ x nₕ

    for i in 1:z_length
        kde_img[:, i] .= get_kde(m2[:, i], grid.m, K = K)
        norm_factor = sum(kde_img[:, i])
        kde_img[:, i] .= kde_img[:, i] ./ norm_factor
    end
    
    fig = Figure()
    ax = Axis(fig[1,1])

    ms = broadcast(trans_utils.m.tf, grid.m)
    hm = heatmap!(ax, ms, zs, kde_img)
    ax.yreversed = true

    cb = fig[1,2]
    Colorbar(cb, hm)

    if return_kde_mat
        return fig, kde_img
    else
        return fig
    end

end

grid_=(m=collect(-1:0.1:5), z= 1:10:1500)

f2 = get_kde_image(mcmc_chain2, modelD2, trans_utils = (m = MT.pow_tf, h = lin_tf), grid= grid_)
ax2 = f2.content[1]
ax2.xscale = log10
f2

cb2 = f2.content[2]
cb2.axis.attributes[:limits][] = (1e-2, 1.)
cb2.axis.attributes[:scale][] = log10
f2

f = get_kde_image(mcmc_chain1, modelD, trans_utils = (m = MT.pow_tf, h = lin_tf))
ax = f.content[1]
ax.xscale = log10
f

cb = f.content[2]
cb.axis.attributes[:limits][] = (1e-2, 1.)
cb.axis.attributes[:scale][] = log10
f

mt_list = MT.get_model_list(mcmc_chain2, modelD2);

fig = Figure()
ax = Axis(fig[1,1])

for i in 1:1000
    plot_model!(fig, mt_list[i], color = :gray, alpha = 0.4)
end

ax.xscale = log10

mt_list = MT.get_model_list(mcmc_chain1, modelD);

fig = Figure()
ax = Axis(fig[1,1])

for i in 1:500
    plot_model!(fig, mt_list[i], color = :gray, alpha = 0.4)
end

ax.xscale = log10


fig3 = Figure()
ax31 = Axis(fig3[1,1])
ax32 = Axis(fig3[2,1])
ax33 = Axis(fig3[1,2])
ax34 = Axis(fig3[2,2])


ax31, ax32 = f.content[1:2]
ax33, ax34 = f2.content[1:2]

fig3

f2 = Figure()
ax = Axis(f2[1,1])

hm = heatmap!(ax, rand(10,10))
cb = f2[1, 2]
Colorbar(cb, hm)


ax.xscale = log10

ax2 = f2.content[1]
ax2.xscale = log10
f2

cb2 = f2.content[2]
cb2.axis.attributes[:limits][] = (1e-2, 1.)
cb2.axis.attributes[:scale][] = log10
f2



