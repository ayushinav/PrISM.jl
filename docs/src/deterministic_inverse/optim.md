# Solvers from `Optimisation.jl`

```@setup occam_demo
using ProEM, LinearAlgebra, CairoMakie, Optimization, OptimizationOptimJL, DifferentiationInterface
```

## Brief introduction

`Optimisation.jl` further provides an extensive [list of solvers](https://docs.sciml.ai/Optimization/stable/optimization_packages/optim/) to solve the inverse problem such as:
* Conjugate Gradient
* Gradient Descent
* LBFGS
* Simulated Annealing
* Particle Swarm

Following, we demonstrate using Conjugate Gradient on a couple of models:

## Demo

### DC Resistivity

We begin by defining a simple synthetic model:
with 5 % error floors for Wenner array.

```@example occam_demo
ρ = log10.([1000.0, 100.0])
h = [2000.0]
m = DCModel(ρ, h)

locs = get_wenner_array(400:200:5000)
resp = forward(m, locs)

err_resp = DCResponse(0.05 .* resp.ρₐ)
```

Now, defining an initial model and covariance matrix:
```@example occam_demo
h_test = fill(60., 50)
m_lbfgs = DCModel(2.0 .+ zeros(length(h_test) + 1), h_test)

C_d = diagm(inv.(err_resp.ρₐ)) .^ 2
nothing # hide
```

and then 

```@example occam_demo
alg_cache = OptAlg(; alg=LBFGS, μ=10.0)
retcode = inverse!(m_lbfgs, resp, locs, alg_cache; W=C_d, max_iters=100);
nothing # hide
```


```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example occam_demo
fig = Figure()
ax_m = Axis(fig[1,1])

plot_model!(ax_m, m, label = "true", color = :steelblue3)
plot_model!(ax_m, m_lm, label = "inverse", color = :tomato)
fig

ax1 = Axis(fig[1,2], xscale = log10)

ab_2 = abs.(locs.srcs[:,2] .- locs.srcs[:,1])./2
plot_response!([ax1], ab_2, resp; plt_type = :scatter, color = :steelblue3)
plot_response!([ax1], ab_2, resp; errs = err_resp, plt_type = :errors, whiskerwidth=10, color = :steelblue3)

resp_lbfgs = forward(m_lbfgs, locs)
plot_response!([ax1], ab_2, resp_lbfgs; color = :tomato)
Legend(fig[2,:], ax_m, orientation = :horizontal)
fig
nothing # hide
```

```@raw html
</details>
```

```@example occam_demo
fig # hide
```

### CSEM : TODO