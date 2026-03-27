## RTO-TKO

RTO-TKO is a stochastic framework introduced in the electromagnetic geophysical community by Blatter et al., 2022 ([a](https://doi.org/10.1093/gji/ggac241) and [b](https://doi.org/10.1093/gji/ggac242))

Similar to fixed discretization scheme, the grid sizes do not change. RTO-TKO tries to obtain the samples of the posterior distribution by unrolling the optimization for the misfit function, that is, it perturbs the model space obtained after one iteration before optimizing for the next. It does this iteratively for both the values of the model space and the regularization coefficient in alternative steps.

Therefore, if we have

```math
m = [m_1, m_2, m_3, ... , m_N]
```

and

```math
d = \mathcal{F}(m)
```

Then for $C_d$, the inverse of data covariance matrix, and $C_m$ the inverse of model covariance matrix, we want to explore uncertainty for the following misfit function:

```math
J(m) = [\mathcal{F}(m) - d]^T C_d [\mathcal{F}(m) - d] + \mu \; m^T C_m m
```

where $\mu$ is the regularization weight and $L$ is the derivative matrix. Probabilistically, the above equation implies that the *a priori* distribution of $m$ is $\mathcal{N}(0, C_m)$. RTO-TKO explores the uncertainty in $m$, as well as the $\mu$-space instead of fixing it, and gives a family of models that fit the data.

The algorithm was proposed for $C_m$ constructed with $L'L$, where $L$ is the discrete derivative matrix. The algorithm works as:

```math
\begin{align*}
\text{Solving for} & \\
& J(m) = [\mathcal{F}(m) - d]^T C_d [\mathcal{F}(m) - d] + \mu (Lm)^T(Lm) \\

& 1) \quad \text{Solve for } \; m_{i+1} \\
& \text{Sample  } \; \tilde{d} \sim \mathcal{N}(d, C_d) \text{ and } \tilde{m} \sim \mathcal{N}(0, \frac{1}{\mu}(L^T L))\\

& \text{Solve} \\
& m_{i+1} = \argmin_{m_{i+1}} \quad [\mathcal{F}(m_{i+1}) - \tilde{d}]^T C_d [\mathcal{F}(m_{i+1}) - \tilde{d}] + \mu_i [L(m_{i+1} - \tilde{m})]^T[L(m_{i+1} - \tilde{m})] \\

& 2) \quad \text{Solve for } \; \mu_{i+1} \\
& \text{Sample  } \; \tilde{d} \sim \mathcal{N}(d, C_d) \\

& \text{Solve} \\
& \mu_{i+1} = \argmin_{\mu_{i+1}} \quad [\mathcal{F}(m_{i+1}) - \tilde{d}]^T C_d [\mathcal{F}(m_{i+1}) - \tilde{d}] - log(p(\mu_{i+1})) \\
\end{align*}
```

!!! note
    
      - Usually, the prior of $\mu$ is a uniform distribution and we do not have to compute the corresponding log pdf term
      - Implementing the above from scratch might not be trivial because of ``L'L`` being non-invertible, and we do the optimization in the domain defined by $\xi = \sqrt{\mu}Lm$.
        Such a variable will then have a standard normal distribution $\mathcal{N}(0, I)$ when $m \sim \mathcal{N}(0, \frac{1}{\mu}(L^T L))$

In the following section, we demonstrate RTO-TKO for a 6-layered earth, including the half-space being imaged using the MT method. The prior distribution is composed of 26 layers, with the initial model being a half-space.

## Demo

```@setup rto_tko
using ProEM
using Distributions
using Turing
using LinearAlgebra
using CairoMakie
```

Let's create a synthetic dataset first, with 10% error floors:

```@example rto_tko
m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3])
f = 10 .^ range(-4; stop=1, length=25)
ω = vec(2π .* f)

r_obs = forward(m_test, ω)

err_phi = asin(0.02) * 180 / π .* ones(length(ω))
err_appres = 0.1 * r_obs.ρₐ
err_resp = MTResponse(err_appres, err_phi)
nothing # hide
```

We don't need to construct the *a priori* the same way as before. Instead, we need to define the optimiser we'd use. Here, we use Occam, and an initial model composed of half-space of 100 $\Omega m$.

```@example rto_tko
z = collect(0:100:2.5e3)
h = diff(z)
m_rto = MTModel(2 .* ones(length(z)), vec(h))

n_samples = 100

r_cache = rto_cache(
    m_rto, [1e-2, 1e4], Occam(), n_samples, n_samples, 1.0, [:ρₐ, :ϕ], false)

rto_chain = stochastic_inverse(
    r_obs, err_resp, ω, r_cache; model_trans_utils=(; m=sigmoid_tf))
```

Since RTO-TKO also samples the regularization coefficient along with model parameters, we exclude it to obtain another chain as:

```@example rto_tko
mt_chain = Chains((rto_chain.value.data[:, 1:(end - 1), :]), [Symbol("m[$i]")
                                                              for i in 1:length(z)])
```

Note that the chain contains fewer samples than we had asked for. This is because a few unstable samples were filtered out. The obtained `mt_chain` contains the *a posteriori* distributions that can be saved using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

```julia
using JLD2
JLD2.@save "file_path.jld2" mt_chain
```

Since RTO-TKO is closely related to Occam, we also include occam results in our figures below. Also, note that we did not require an *a priori* distribution to obtain samples. However, we would need to create the same to plot our posterior samples. This can just be something that encapsulates the whole prior space (or an envelope of the same), e.g., in the case of magnetotelluric imaging, we know that the resistivity values will always be in [-2, 5] on the log-scale.

```@example rto_tko
modelD = MTModelDistribution(Product([Uniform(-1.0, 5.0) for i in eachindex(z)]), vec(h))
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rto_tko
fig = Figure()
ax = Axis(fig[1, 1])
hm = get_kde_image!(ax, mt_chain, modelD; kde_transformation_fn=log10,
    colormap=:binary, colorrange=(-3.0, 0.0))
Colorbar(fig[1, 2], hm; label="log pdf")

mean_kws = (; color=:seagreen3, linewidth=2)
std_kws = (; color=:red, linewidth=1.5)
get_mean_std_image!(ax, mt_chain, modelD; confidence_interval=0.99, mean_kwargs=mean_kws,
    std_plus_kwargs=std_kws, std_minus_kwargs=std_kws)
ylims!(ax, [2500, 0])

plot_model!(ax, m_test; color=:black, linestyle=:dash, label="true", linewidth=2)
Legend(fig[2, :], ax; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example rto_tko
fig # hide
```

The list of models can then be obtained from chains using

```@example rto_tko
model_list = get_model_list(mt_chain, modelD)
nothing # hide
```

We can then easily check the fit of the response curves

```@example rto_tko
fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])

resp_post = forward(model_list[1], ω);
for i in 1:(length(model_list) > 100 ? 100 : length(model_list))
    forward!(resp_post, model_list[i], ω)
    plot_response!([ax1, ax2], ω, resp_post; alpha=0.4, color=:gray)
end
plot_response!([ax1, ax2], ω, r_obs; errs=err_resp, plt_type=:errors, whiskerwidth=10)
plot_response!([ax1, ax2], ω, r_obs; plt_type=:scatter, label="true")
ylims!(ax1, exp10.([1.2, 3.1]))
ylims!(ax2, [12, 70])

fig
```

!!! warning
    
    It is recommended that the samples from RTO-TKO can be filtered out to reject the samples that have a poor fit on the data.
