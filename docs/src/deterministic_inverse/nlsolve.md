# Solvers from `NonlinearSolve.jl`

```@setup occam_demo
using ProEM, LinearAlgebra, CairoMakie, NonlinearSolve
```
## Brief introduction

We can solve the inverse problem using a suite of solvers provided by `NonlinearSolve.jl` such as :
* Gauss Newton
* Levenberg Marquardt
* Newton Raphson

## Demo

We demonstrate using Levenberg-Marquadt on a couple of geophysical models below:

### Love waves 

Occam1D was initially designed and compiled for 1D MT, and this is the first example we run here.

Let's define a synthetic model:

```@example occam_demo
ρ = log10.([100.0, 10.0, 400.0, 1000.0])
h = [100.0, 1000.0, 10000.0]
m = MTModel(ρ, h)
```

and then get data, and assume 5% error floors:

```@example occam_demo
T = 10 .^ (range(-1, 5; length=19))
ω = 2π ./ T
resp = forward(m, ω)

err_resp = MTResponse(
    0.1 .* resp.ρₐ,
    180 / π .* asin(0.1) .+ zero(ω)
)
```

We need to define a data covariance matrix $C_d$ and an initial model. Let's assume gaussian noise in this case, and choose an initial model a half-space of 100 $\Omega m$.

```@example occam_demo
C_d = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2; # weight matrix

# initial model
h_test = 10 .^ range(0.0, 5.0; length=50)
ρ_test = 2.0 .+ randn(length(h_test) + 1)

m_lm = MTModel(ρ_test, h_test);
```

The final result will be stored in the same variable `m_occam`. All we need to do now is specify using Occam and then calling `inverse!`.

!!! note
    Using Occam, you also have the option to perform a smoothing step. Once the model has achieved the threshold misfit, it is smoothened until it fits the data just about the threshold misfit.

```@example occam_demo
alg_cache = NonlinearAlg(; alg=LevenbergMarquardt(), μ=500.0)
retcode = inverse!(m_lm, resp, ω, alg_cache; W=C_d, max_iters=100, verbose=true);
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example occam_demo
fig = Figure()
ax_m = Axis(fig[1:2,1], yscale = log10)

plot_model!(ax_m, m, label = "true", color = :steelblue3)
plot_model!(ax_m, m_lm, label = "LM", color = :tomato)
fig

ax1 = Axis(fig[1,2], xscale = log10)
ax2 = Axis(fig[2,2], xscale = log10)

plot_response!([ax1, ax2], T, resp; plt_type = :scatter, color = :steelblue3)
plot_response!([ax1, ax2], T, resp; errs = err_resp, plt_type = :errors, whiskerwidth=10, color = :steelblue3)

resp_lm = forward(m_lm, ω)
plot_response!([ax1, ax2], T, resp_lm; color = :tomato)
Legend(fig[3,:], ax_m, orientation = :horizontal)

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
fig # hide
```