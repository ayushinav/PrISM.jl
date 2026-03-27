# Solvers from `NonlinearSolve.jl`

```@setup nl_demo
using ProEM, LinearAlgebra, CairoMakie, NonlinearSolve, DifferentiationInterface
```

## Brief introduction

We can solve the inverse problem using a [suite of solvers](https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_system_solvers/#Full-List-of-Methods) provided by `NonlinearSolve.jl` such as :

  - Gauss Newton
  - Levenberg Marquardt
  - Newton Raphson

## Demo

!!! warning
    
    Owing to the non-uniqueness of the geophysical methods, it might take considerable effort to obtain a model that converges and/or fits the data.

We demonstrate using Levenberg-Marquadt on a couple of geophysical models below:

### Rayleigh waves

We first begin with an example for Rayleigh waves

Let's define a synthetic model:

```@example nl_demo
vs = [3.2, 4.39731, 4.40192, 4.40653, 4.41113, 4.41574]
vp = [5.8, 8.01571, 7.99305, 7.97039, 7.94773, 7.92508]
ρ = [2.6, 3.38014, 3.37797, 3.37579, 3.37362, 3.37145]
h = [20.0, 20.0, 20.0, 20.0, 20.0] .* 1e3
m = RWModel(vs, h, ρ, vp)
```

and then get the forward response, with 1% error floors :

```@example nl_demo
freq = exp10.(-2:0.1:1)
t = inv.(freq)

resp = forward(m, t)

err_resp = SurfaceWaveResponse(0.01 .* resp.c,)
```

We need to define a data covariance matrix $C_d$ and an initial model. Let's assume gaussian noise in this case, and choose an initial model a half-space of 100 $\Omega m$.

```@example nl_demo
h_test = fill(2500.0, 50)
m_lm = RWModel(4.0 .+ zeros(length(h_test) + 1), h_test, fill(3e3, 51), fill(7e3, 51))

C_d = diagm(inv.(err_resp.c)) .^ 2
nothing # hide
```

The final result will be stored in the same variable `m_occam`. All we need to do now is specify using Occam and then calling `inverse!`.

!!! note
    
    Using Occam, you also have the option to perform a smoothing step. Once the model has achieved the threshold misfit, it is smoothened until it fits the data just about the threshold misfit.

```@example nl_demo
alg_cache = NonlinearAlg(; alg=LevenbergMarquardt(; autodiff=AutoFiniteDiff()), μ=100.0)
rw_bounds = transform_utils(
    x -> SubsurfaceCore.sigmoid(x, 3.0, 5.0), x -> SubsurfaceCore.inverse_sigmoid(x, 3.0, 5.0))
C_m = Float64.(I(51))
retcode = inverse!(
    m_lm, resp, t, alg_cache; W=C_d, max_iters=100, verbose=10, model_trans_utils=rw_bounds);
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example nl_demo
fig = Figure()
ax_m = Axis(fig[1, 1]; xlabel="Vs (km/s)", ylabel="depth (m)")

plot_model!(ax_m, m; label="true", color=:steelblue3)
plot_model!(ax_m, m_lm; label="occam", color=:tomato)
fig

ax1 = Axis(fig[1, 2]; xscale=log10)

plot_response!([ax1], t, resp; plt_type=:scatter, color=:steelblue3)
plot_response!(
    [ax1], t, resp; errs=err_resp, plt_type=:errors, whiskerwidth=10, color=:steelblue3)

resp_lm = forward(m_lm, t)
plot_response!([ax1], t, resp_lm; color=:tomato)
Legend(fig[2, :], ax_m; orientation=:horizontal)
fig
nothing # hide
```

```@raw html
</details>
```

```@example nl_demo
fig # hide
```

This one did not converge, and that's geophysical inversion 90% of the time. :)

### DC Resistivity

Let's try and see if we can get things to converge for a simple DC resistivity model. Let's make another synthetic model with 5 % error floors for Wenner array.

```@example nl_demo
ρ = log10.([1000.0, 100.0])
h = [2000.0]
m = DCModel(ρ, h)

locs = get_wenner_array(400:200:5000)
resp = forward(m, locs)

err_resp = DCResponse(0.05 .* resp.ρₐ)
```

Now, defining an initial model and covariance matrix:

```@example nl_demo
h_test = fill(60.0, 50)
m_lm = DCModel(2.0 .+ zeros(length(h_test) + 1), h_test)

C_d = diagm(inv.(err_resp.ρₐ)) .^ 2
nothing # hide
```

and then

```@example nl_demo
alg_cache = NonlinearAlg(; alg=LevenbergMarquardt(; autodiff=AutoFiniteDiff()), μ=10.0)
retcode = inverse!(m_lm, resp, locs, alg_cache; W=C_d, max_iters=100);
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example nl_demo
fig = Figure()
ax_m = Axis(fig[1, 1])

plot_model!(ax_m, m; label="true", color=:steelblue3)
plot_model!(ax_m, m_lm; label="inverse", color=:tomato)
fig

ax1 = Axis(fig[1, 2]; xscale=log10)

ab_2 = abs.(locs.srcs[:, 2] .- locs.srcs[:, 1]) ./ 2
plot_response!([ax1], ab_2, resp; plt_type=:scatter, color=:steelblue3)
plot_response!(
    [ax1], ab_2, resp; errs=err_resp, plt_type=:errors, whiskerwidth=10, color=:steelblue3)

resp_lm = forward(m_lm, locs)
plot_response!([ax1], ab_2, resp_lm; color=:tomato)
Legend(fig[2, :], ax_m; orientation=:horizontal)
fig
nothing # hide
```

```@raw html
</details>
```

```@example nl_demo
fig # hide
```
