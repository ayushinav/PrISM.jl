## Variable discretization

There are a few cases when we can simplify the subsurface structure by expressing it using only a few layers. In those cases, it may be desirable to vary the cell size (layer thickness) along with the model parameters. This allows us to move the layer interfaces up and down while parameterizing the model space in a different way wherein there is more flexibility.

Let's denote the model parameters, eg., conductivity, by `m`, and the layer thickness by `h`. Therefore, in a N-layer case, we will have

```math
m = [m_1, m_2, m_3, ... , m_N] \\
h = [h_1, h_2, h_3, ... , h_{N_1}]
```

such that

```math
m_i  \in \mathcal{D}_{m_i} \text{ ; where } \mathcal{D}_{m_i} = \textit{a priori} \text{ distribution for } m_i \text{ ; and} \\
 h_i  \in \mathcal{D}_{h_i} \text{ ; where } \mathcal{D}_{h_i} = \textit{a priori} \text{ distribution for } h_i 
```

In the example that follows, we demonstrate MCMC inversion for a 6-layered earth, including the half-space being imaged using Rayleigh waves. THe prior distribution assumes all layers have uncorrelated shear wave velocities bounded between $3.5$ and $5 \; km/s$, defined using a uniform distribution, and layered thickness values vary between 15 and 25 km for all the layers above the half-space.

In the following example, we demonstrate MCMC inversion for a 3-layered earth, including the half-space being imaged using DC resistivity method. The prior distribution assumes all layers have uncorrelated resistivities bounded between $10^{-1}$ and $10^5$, defined using a uniform distribution.

!!! tip Important

    The most important thing to be noted here is the specification of the prior distribution, done via:

    ```julia
    modelD = RWModelDistribution(
        Product(
            [Uniform(-1.0, 5.0) for i in eachindex(z)]
        ),
        Product(
            [Uniform(h_lb[i], h_ub[i]) for i in eachindex(h)]
        )
        ...
    );
    ```

    in the example.

## Copy-Pasteable code

```@setup variable_mcmc
using PrISM
using Distributions
using Turing
using LinearAlgebra
using CairoMakie
```

Let's create a synthetic dataset first, with 1% error floors:

```@example variable_mcmc
vs = [4.1, 4.4, 4.8]
vp = [7.0, 7.5, 7.8]
ρ = [3.0, 3.2, 3.3]
h = [20.0, 20.0] .* 1e3
m_test = RWModel(vs, h, ρ, vp)

T = exp10.(0:0.1:3)

r_obs = forward(m_test, T);

err_c = 0.01 * r_obs.c;
err_resp = SurfaceWaveResponse(err_c)
nothing # hide
```

Now, let's define the *a priori* distribution with the variable grid points. Note that we only invert for one physical property, i.e., shear wave velocity.

```@example variable_mcmc
z = collect(0e3:20e3:40e3)
h = diff(z)
h_ub = h .+ 2.5e3
h_lb = h .- 2.5e3

# variable discretization
modelD = RWModelDistribution(Product([Uniform(4.0, 5.0) for i in eachindex(z)]),
    Product([Uniform(h_lb[i], h_ub[i]) for i in eachindex(h)]),
    fill(3.3, length(z)), fill(7.5, length(z)))
nothing # hide
```

then define the likelihood

```@example variable_mcmc
respD = SurfaceWaveResponseDistribution(normal_dist)
nothing # hide
```

Put everything together for MCMC

```@example variable_mcmc
n_samples = 10_000
mcache = mcmc_cache(modelD, respD, n_samples, MH())

rw_chain = stochastic_inverse(r_obs, err_resp, T, mcache; progress=true)
```

The obtained `rw_chain` contains the *a posteriori* distributions that can be saved using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

```julia
using JLD2
JLD2.@save "file_path.jld2" rw_chain
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example variable_mcmc
fig = Figure()
ax = Axis(fig[1, 1])
hm = get_kde_image!(ax, rw_chain, modelD; kde_transformation_fn=log10,
    grid=(m=collect(4.0:0.01:5.0), z=collect(1:1e3:60e3)),
    colormap=:binary, colorrange=(-4, 0.0))
Colorbar(fig[1, 2], hm; label="log pdf")

mean_kws = (; color=:seagreen3, linewidth=2)
std_kws = (; color=:red, linewidth=1.5)
get_mean_std_image!(ax, rw_chain, modelD; confidence_interval=0.9, mean_kwargs=mean_kws,
    std_plus_kwargs=std_kws, std_minus_kwargs=std_kws, z_points=collect(1:1e3:60e3))
ylims!(ax, [6e4, 0])

plot_model!(ax, m_test; color=:black, linestyle=:dash, linewidth=2, label="true")
Legend(fig[2, :], ax; orientation=:horizontal)
```

```@raw html
</details>
```

```@example variable_mcmc
fig # hide
```

The list of models can then be obtained from chains using, which can then be used to check the fits and perform other diagnostics:

```@example variable_mcmc
model_list = get_model_list(rw_chain, modelD)
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example variable_mcmc
fig = Figure()
ax1 = Axis(fig[1, 1])

resp_post = forward(model_list[1], T);
for i in 1:(length(model_list) > 1000 ? 1000 : length(model_list))
    forward!(resp_post, model_list[i], T)
    plot_response!([ax1], T, resp_post; alpha=0.4, color=:gray)
end

plot_response!([ax1], T, r_obs; errs=err_resp, plt_type=:errors, whiskerwidth=10)
plot_response!([ax1], T, r_obs; plt_type=:scatter, label="true")
```

```@raw html
</details>
```

```@example variable_mcmc
fig # hide
```
