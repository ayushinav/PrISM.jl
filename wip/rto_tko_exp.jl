#=
 TODO
reg term in occam and inverse
NonlinearSolve.jl
=#

using Pkg
Pkg.activate("wip/.")

using Revise
using MT
using LinearAlgebra
using BenchmarkTools
using Distributions
using Turing


# using JLD2

# sts = ["2_yx", "3_yx", "6_yx", "4_yx"]
# # for st in sts 
# path_to_kivu = "../../kivu/"
# st = "2_yx"
# kk = jldopen(path_to_kivu * "filter_data_err_floor/"*st*".jld2")
# r_obs = kk["r_obs4"];
# err_resp = kk["resp_err4"];
# ω = kk["ω"];


# idxs = sortperm(ω);
# r_obs.ρₐ .= r_obs.ρₐ[idxs];
# r_obs.ϕ .= r_obs.ϕ[idxs];
# err_resp.ρₐ .= err_resp.ρₐ[idxs];
# err_resp.ϕ .= err_resp.ϕ[idxs];
# ω .= ω[idxs];

# z = 10 .^collect(range(1, stop = 6, length = 100));
# h = diff(z);


# @show log10.(r_obs.ρₐ)
# @show r_obs.ϕ
# @show 1/log(10) * 1 ./ r_obs.ρₐ .* err_resp.ρₐ
# @show err_resp.ϕ

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

W = Diagonal(
    vcat([inv.(getfield(err_resp, k)) for k in propertynames(err_resp)]...)
) .^ 2;

model_occam = MTModel(fill(2., length(z)), vec(h));

@time retcode = inverse!(model_occam, r_obs, ω, 
    Occam(μgrid = [1e-2, 1e2]), W = W, χ2 = 1.0,
    max_iters = 50,
    verbose = false,
    # response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)
)

plot_model(model_occam, yscale = :identity, ylims = (0, 1200), label = false, xlims = (1e-2, 1e6))
plt_occ_resp = prepare_plot(forward(model_occam, ω), ω, plot_type = :scatter, color = "green")
# plt_occ_resp = prepare_plot(forward(model_occam, ω), ω, plot_type = :scatter, label = false, color = "green")
prepare_plot!(plt_occ_resp, r_obs, ω, err_resp, label = "r_obs", color = "blue")
prepare_plot!(plt_occ_resp, r_obs, ω, plot_type = :plot, label = false, color = "blue")
plot_response(plt_occ_resp)

m_rto = MTModel(fill(2.0, length(z)), h);

# sigmoid_tf2 = MT.transform_utils([-2.0, 6.0], sigmoid, inverse_sigmoid, d_sigmoid);

#=
v1 : -2,6
v3 : -4, -2
=#


m_rto = MTModel(fill(2.0, length(z)), h);

rto_c = MT.rto_cache(m_rto, [1e-6, 1e2], Occam(), 50, 400, 1.0, [:ρₐ, :ϕ], false)
# rto_c2 = MT.rto_cache(m_rto, [1e-2, 1e6], NonlinearAlg(;alg = LevenbergMarquardt), 50, 50, 1.0, [:ρₐ, :ϕ], false)

@time rto_chains = stochastic_inverse(r_obs, err_resp, ω, rto_c, verbose = 50); #, response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)) #, trans_utils = (;m = sigmoid_tf2, h = MT.lin_tf));
@profview rto_chains = stochastic_inverse(r_obs, err_resp, ω, rto_c) #, response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)) #, trans_utils = (;m = sigmoid_tf2, h = MT.lin_tf));
# @time rto_chains2 = stochastic_inverse(r_obs, err_resp, ω, rto_c2) #, trans_utils = (;m = sigmoid_tf2, h = MT.lin_tf));
# @profview rto_chains = stochastic_inverse(r_obs, err_resp, ω, rto_c, trans_utils = (;m = MT.sigmoid_tf, h = MT.lin_tf));
# rto_chains = stochastic_inverse(r_obs, err_resp, ω, rto_c)

# m_gn = MTModel(1. *ones(61) .+ (0.00 .* randn(61)), 20. *ones(60))
# inverse!(m_gn, r_obs, ω, NonlinearAlg(;alg = LevenbergMarquardt, μ = 1e-2), W = W, max_iters = 100, χ2 = 1.) #, response_trans_utils = (ρₐ = log_tf,ϕ = MT.lin_tf))

# plot_model!(plt_mods2, m_gn)
# ch_test = rto_chains;

# idcs = broadcast(!isnan, view(ch_test.value.data, :, 1, 1))
# rto_chains2 = Turing.Chains(ch_test.value.data[idcs, :,1], [Symbol("m[$i]") for i in 1:(61 + 1)])

m_dist = MTModelDistribution(Product([Uniform(-1, 5) for i = 1:length(z)]), h)

# rto_chains.value.data[111,:,1] .= rto_chains.value.data[1,:,1];
mt_chain = Turing.Chains(
        (rto_chains.value.data[:,1:length(z),:]), [Symbol("ρ[$i]") for i in 1:length(z)]);

# mt_chain = Turing.Chains(
#     (rto_chains2.value.data[:,1:length(z),:]), [Symbol("ρ[$i]") for i in 1:length(z)]);

mt_model_list = get_model_list(mt_chain, m_dist);
# mt_model_list[111] = copy(mt_model_list[1]);
n_samples = length(mt_model_list)

# using JLD2

# JLD2.jldsave("rto_tko_test_narrow.jld2"; mt_model_list, rto_chains);

using Plots

plt_mods1 = Plots.plot()
for i in 1:(length(mt_model_list) > 500 ? 500 : length(mt_model_list))
    plot_model!(plt_mods1, mt_model_list[i], label = false, color = :green, alpha = 0.4);
    # prepare_plot!(plt_resps, resp_models, ω, alpha = 0.4, label = false, color = :gray); 
end

plt_mods1

plot_model!(plt_mods1, model_occam, yscale = :identity, ylim = (0, 1200), yticks = collect(0:100:1200), 
        label = false, xlim = (0.01, 1e6), xticks = 10 .^ (-2:1.:6), linewidth = 3, color = :black, alpha = 1.)

# plot!(plt_mods, xlims = (1e-1, 1e3))

# plt_mods2 = Plots.plot()
# for i in 1:(length(mt_model_list) > 500 ? 500 : length(mt_model_list))
#     plot_model!(plt_mods2, mt_model_list[i], label = false, color = :green, alpha = 0.4);
#     # prepare_plot!(plt_resps, resp_models, ω, alpha = 0.4, label = false, color = :gray); 
# end

# plt_mods2

# plot_model!(plt_mods2, model_occam, yscale = :identity, ylim = (0, 1200), yticks = collect(0:100:1200), 
#         label = false, xlim = (0.01, 1e6), xticks = 10 .^ (-2:1.:6), linewidth = 3, color = :black, alpha = 1.)

r2 = copy(r_obs);
rms_vals = zeros(length(mt_model_list))
for i in eachindex(mt_model_list)
    forward!(r2, mt_model_list[i], ω)
    # rms_vals[i] = norm(([r2.ρₐ..., r2.ϕ...] .- [r_obs.ρₐ..., r_obs.ϕ...]) .* inv.([err_resp.ρₐ..., err_resp.ϕ...]))
    # rms_vals[i] = sqrt(rms_vals[i]/(2*length(ω)))
    rms_vals[i] = MT.χ²([r2.ρₐ..., r2.ϕ...], [r_obs.ρₐ..., r_obs.ϕ...], W=W)
    # @show rms_vals2, rms_vals[i]
end
# rms_vals .= sqrt.(rms_vals ./(2*length(ω)))
# rms_vals .= sqrt.(rms_vals ./(20))

histogram(rms_vals, xlims = (0.5, 2.5), ylabel = "pdf", xlabel = "rms", label = false, title = "rms histogram", normalize = :pdf)
# savefig("rms_vals.png")

histogram(log10.(rto_chains.value.data[:,end,1]), normalize =:pdf, label = false, xlabel = "μ", ylabel = "pdf")
# savefig("μ_vals.png")

plt_resps = prepare_plot(r_obs, ω, alpha = 0., label = false)
resp_models = forward(mt_model_list[1], ω);

for i in 1:(length(mt_model_list) > 200 ? 200 : length(mt_model_list))
    forward!(resp_models, mt_model_list[i], ω);
    prepare_plot!(plt_resps, resp_models, ω, alpha = 0.1, label = false, color = :violet); 
end

prepare_plot!(plt_resps, r_obs, ω, err_resp, markersize = 4, 
    msc = :orange, color =:red, label = false, falpha = 0.);
plt_resp_ = plot_response(plt_resps, size = (1000, 700), margin = 0Plots.mm)

# trans_utils = (m = MT.lin_tf, h = MT.lin_tf);
# pre_img = pre_image(m_dist, mt_chain, trans_utils = trans_utils);

# grid_ = (;m = collect(-2:0.1:6), h = z)
h2 = 3. .*ones(200);
z2 = [0; cumsum(h2)];

# z2 = 10 .^ range(1, 3.5, length = 200)
# h2 = diff(z2)
grid_ = (;m = collect(-2:0.1:6), h = h)
pow_tf2 = transform_utils(
    Vector{Float32}([]), (x, p) -> 10^x, (x, p) -> log10(x), (x, p) -> (10^x * log(10)));
trans_utils_ = (;m = pow_tf2, h = MT.lin_tf);
pre_img = pre_image(m_dist, mt_chain , grid = grid_, trans_utils = trans_utils_);
# kde_vals, kde_img = get_kde_image(pre_img..., true, xscale = :identity, yscale = :identity, yflip = true, clims = (0., 0.1))# , ylim = (0, 1200))
kde_img = get_kde_image(
    pre_img..., 
    return_vals = false,
    xscale = :log10, yscale = :identity, #yticks = collect(100:200:1200),
    yflip = true, 
    # xlim = (1f-2, 1f5)
    # clims = (0., 0.4), color = :jet,
    # size = (500, 800)
)
# kde_img = get_kde_image(pre_img..., false, xscale = :identity, yscale = :identity, clims = (0., 0.1) , ylim = (1f-1, 1200))
# kde_vals, kde_img = get_

m2_ = [-2.05, -1.95, -1.85, -1.75, -1.65, -1.55, -1.45, -1.35, -1.25, -1.15, -1.05, -0.95, -0.8500000000000001, -0.75, -0.6499999999999999, -0.55, -0.45, -0.35, -0.25, -0.15000000000000002, -0.05, 0.05, 0.15000000000000002, 0.25, 0.35, 0.45, 0.55, 0.6499999999999999, 0.75, 0.8500000000000001, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.1500000000000004, 2.25, 2.3499999999999996, 2.45, 2.55, 2.6500000000000004, 2.75, 2.8499999999999996, 2.95, 3.05, 3.1500000000000004, 3.25, 3.3499999999999996, 3.45, 3.55, 3.6500000000000004, 3.75, 3.8499999999999996, 3.95, 4.05, 4.15, 4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05, 5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 6.05]
h2_ = [0.10000000149011612, 20.100000001490116, 40.100000001490116, 60.100000001490116, 80.10000000149012, 100.10000000149012, 120.10000000149012, 140.10000000149012, 160.10000000149012, 180.10000000149012, 200.10000000149012, 220.10000000149012, 240.10000000149012, 260.1000000014901, 280.1000000014901, 300.1000000014901, 320.1000000014901, 340.1000000014901, 360.1000000014901, 380.1000000014901, 400.1000000014901, 420.1000000014901, 440.1000000014901, 460.1000000014901, 480.1000000014901, 500.1000000014901, 520.1000000014901, 540.1000000014901, 560.1000000014901, 580.1000000014901, 600.1000000014901, 620.1000000014901, 640.1000000014901, 660.1000000014901, 680.1000000014901, 700.1000000014901, 720.1000000014901, 740.1000000014901, 760.1000000014901, 780.1000000014901, 800.1000000014901, 820.1000000014901, 840.1000000014901, 860.1000000014901, 880.1000000014901, 900.1000000014901, 920.1000000014901, 940.1000000014901, 960.1000000014901, 980.1000000014901, 1000.1000000014901, 1020.1000000014901, 1040.1000000014901, 1060.1000000014901, 1080.1000000014901, 1100.1000000014901, 1120.1000000014901, 1140.1000000014901, 1160.1000000014901, 1180.1000000014901, 1200.1000000014901, 1240.1000000014901]
heatmap(m2_, h2_, kde_vals, yflip = true, xlims = (-1, 5))


get_kde_image(pre_img..., return_vals = false, xscale = :log10, yflip = true, clims = (0., 0.2), ylim = (0, 1200))
# plt_dist = heatmap(pre_img[2][:m], pre_img[1][:h], log10.(kde_vals .+ 1e-9), yflip = true, clims = (-3, 0), yscale = :identity, ylim = (0, 1200))
# plt_dist = heatmap(pre_img[2][:m], pre_img[1][:h], log10.(kde_vals .+ 1e-9), clims = (-3, 0))

plt_dist = heatmap(10 .^collect(-2:0.1:6), z2, log10.(kde_vals .+ 1e-9), clims = (-2.5, -0.5), ylim = (0, 1200), xscale =:log10, xlims = (1e-2, 1e6), color = :jet)

mean_m = mean(rto_chains.value.data[:,:,1], dims = 1)[1:end-1]
std_m = std(rto_chains.value.data[:,:,1], dims = 1)[1:end-1]

mean_model = MTModel(mean_m, h);
model_95 = MTModel(mean_m .+ std_m .* 1.96, h);
model_5 = MTModel(mean_m .- std_m .* 1.96, h);

plot_model!(plt_dist, mean_model, yscale =:identity, ylims = (0,1200), xlims = (1e-2, 1e6), color = "black", label = false)
plot_model!(plt_dist, model_95, yscale =:identity, ylims = (0,1200), xlims = (1e-2, 1e6), color = "red", label = false)
plot_model!(plt_dist, model_5, yscale =:identity, ylims = (0,1200), xlims = (1e-2, 1e6), color = "red", label = false)
plot!(plt_dist, size = (500, 700), xlabel = "ρ (Ωm)", ylabel = "depth (m)", margin = 3Plots.mm, xticks = 10. .^collect(-2:6), yticks = 0:100:1200)
# savefig("posteriors.png")
plt_dist = kde_img

grid_ = (m = collect(-1:0.1:5), h = h) #diff(z))
scale_tf = transform_utils([], (x) -> x/1e3, (x) -> x*1e3, (x) -> 1e-3);
trans_utils = (m = lin_tf, h = MT.lin_tf,);
pre_img = pre_image(m_dist, mt_chain, grid = grid_, trans_utils = trans_utils);
mean_std_plt_lin = get_mean_std_image(
    pre_img..., trans_utils = trans_utils, 
    yscale = :log10, ylim = (1, 1.2e3), yticks = collect(100:200:1200),
    xlim = (10. ^ -1.5, 1e5), xticks = 10.0 .^ collect(-1:5), 
    ylabel = "depth (km)"
)
plot(mean_std_plt_lin, yscale = :identity)


lay_out = @layout [Plots.grid(2,1) b{0.6w}]; 
plot(plt_resp_, mean_std_plt_lin, kde_img, layout = lay_out, size = (1000, 1000), 
xtickfontsize=15, xguidefontsize = 18, ytickfontsize=15, yguidefontsize = 18,
leftmargin = 7Plots.mm, rightmargin = 5Plots.mm, topmargin = 3Plots.mm)

savefig(path_to_kivu * "rto_tko_1k_$(st).png")
# end

kde_vals, kde_img = get_kde_image(pre_img..., true, xscale = :log10, yscale = :log10, yflip = true, clims = (0., 0.3))
# kde_img = get_kde_image(pre_img..., xscale = :log10, yscale = :log10, yflip = true, clims = (0., 0.2))

plt_dist = heatmap(pre_img[2][:m], pre_img[1][:h], log10.(kde_vals .+ 1e-9), yflip = true, clims = (-3, 0), yscale = :identity, ylim = (0, 1200))
plt_dist = heatmap(pre_img[2][:m], pre_img[1][:h], log10.(kde_vals .+ 1e-9), clims = (-3, 0))

# plt_dist = heatmap(pre_img[2][:m][1:30], pre_img[1][:h][1:30], log10.(kde_vals[1:30, 1:30] .+ 1e-9), yflip = true, clims = (-3, 0), yscale = :identity, ylim = (0, 1200))

lay_out = Plots.@layout [Plots.grid(2,1) b{0.6w}]; 

plot(plt_resp_, plt_mods, plt_dist, layout = lay_out, 
    size = (1600, 1000)
)

savefig("narrow_bounds.png")

# arr = zeros(100,100);
# for i in 1:100
#     arr[i,:] .= i
# end

# heatmap(1:100, 1:100, arr, yflip = true)

# heatmap(1:100, 1:100, arr, yflip = true, ylim = (1,50))
h = 20. .*ones(60);
z = [0; cumsum(h)];

z = 10 .^ range(-1, 3.5, 61)
h = diff(z)

W = Diagonal(
    vcat([inv.(getfield(err_resp, k)) for k in propertynames(err_resp)]...)
) .^ 2;

model_occam = MTModel(fill(2., length(z)), vec(h));

@time retcode = inverse!(model_occam, r_obs, ω, 
    Occam(μgrid = [1e-2, 1e2]), W = W, χ2 = 1.0,
    max_iters = 50,
    verbose = true,
    # response_trans_utils = (;ρₐ = log_tf, ϕ = MT.lin_tf)
)

plot_model(model_occam)

plot_model(MTModel(
    [2., 1.], [100.]
), ylim = (1, 200), yscale = :identity)

plot_model(model_occam, yscale = :identity, ylims = (0, 1200), label = false, xlims = (1e-2, 1e6))
plot_model(model_occam, yscale = :log10, ylims = (1, 1200), label = false, xlims = (1e-2, 1e6))


arr = [10, 20, 30, 40]

cumsum(arr)- arr/2

arr2 = zeros(eltype(arr), 2*length(arr)+1)

arr2[3:2:2length(arr)+1] .= cumsum(arr)
arr2[2:2:2length(arr)] .= cumsum(arr)- arr/2


arr = xs

# xs = Float64[1,2,3,5,9,10]
xs = [0, 0.05, 0.1, 9.9, 10.]
xs2 = [0, 0.05, 0.1, 9.9, 10., 10.1]
# xs2 = Float64[0,1,2,3,4]
ys = Float64[1,2,3,4,5,6,7,8]
ys2 = Float64[0,1,2,3,4,5,6,7,8]
zs = xs' .+ ys
heatmap(xs, ys, zs)
heatmap(xs2, ys2, zs)

xs = Float64[1,2,3,4]
ys = Float64[1,2,3,4,5,6,7,8]
zs = xs'./10 .+ ys
p1 = heatmap(xs, ys, zs, title = "original")
p2 = heatmap(xs, ys, zs, ylims = (2,4), title = "ylim")
p3 = heatmap(xs, ys, zs, yflip = true, title = "yflip")
p4 = heatmap(xs, ys, zs, ylims = (2,4), yflip = true, title = "yflip + ylim")

plot(p1, p2, p3, p4)
savefig("pltos_issue.png")